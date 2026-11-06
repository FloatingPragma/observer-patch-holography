"""Identify supplied strain responses from finite, error-bounded records.

The native E instrument supplies a charge census. Volume interventions and
their energy readings are separate supplies; they are never native records.
All certificate arithmetic, including logarithm enclosures, is rational.
"""
from __future__ import annotations

import argparse
import ast
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path

from ..native_equilibration import history_entropy as census

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "results.json"


def exact(value):
    if type(value) not in (int, F):
        raise ValueError("exact integer or rational required")
    return F(value)


@dataclass(frozen=True)
class Interval:
    lo: F
    hi: F

    def __post_init__(self):
        object.__setattr__(self, "lo", exact(self.lo))
        object.__setattr__(self, "hi", exact(self.hi))
        if self.lo > self.hi:
            raise ValueError("reversed interval")

    def __add__(self, other):
        return Interval(self.lo + other.lo, self.hi + other.hi)

    def __sub__(self, other):
        return Interval(self.lo - other.hi, self.hi - other.lo)

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __mul__(self, other):
        values = [a*b for a in (self.lo, self.hi) for b in (other.lo, other.hi)]
        return Interval(min(values), max(values))

    def __truediv__(self, other):
        if other.lo <= 0 <= other.hi:
            raise ValueError("division interval contains zero")
        return self * Interval(1/other.hi, 1/other.lo)

    def contains(self, value):
        return self.lo <= value <= self.hi


def point(value):
    return Interval(value, value)


def log_bounds(value, terms=24):
    """Rational enclosure using range reduction and a positive atanh series."""
    x = exact(value)
    if x <= 0 or type(terms) is not int or not 1 <= terms <= 128:
        raise ValueError("positive logarithm argument and 1..128 terms required")
    if x == 1:
        return point(0)
    if x < 1:
        return -log_bounds(1/x, terms)

    def reduced(y):
        z = (y-1)/(y+1)
        lower = 2*sum((z**(2*k+1)/(2*k+1) for k in range(terms)), F())
        tail = 2*z**(2*terms+1)/((2*terms+1)*(1-z*z))
        return Interval(lower, lower+tail)

    power = 0
    while x >= 2:
        x /= 2
        power += 1
    return point(power)*reduced(F(2)) + reduced(x)


@dataclass(frozen=True)
class EnergyRecord:
    volume: F
    energy: F
    error: F = F()
    state_tag: str = "held-canonical-state"

    def __post_init__(self):
        for name in ("volume", "energy", "error"):
            object.__setattr__(self, name, exact(getattr(self, name)))
        if self.volume <= 0 or self.error < 0:
            raise ValueError("positive volume and nonnegative error required")
        if type(self.state_tag) is not str or not self.state_tag:
            raise ValueError("a declared held-state identifier is required")

    @property
    def interval(self):
        return Interval(self.energy-self.error, self.energy+self.error)


def _records(records, size):
    if (type(records) is not tuple or len(records) != size
            or any(type(row) is not EnergyRecord for row in records)
            or len({row.state_tag for row in records}) != 1):
        raise ValueError("complete ordered records of the same held state required")


def identify_power_response(records):
    """Decode H(V)=C+A(V/V0)^(-alpha), A>0, alpha>=0.

    The response family is a hypothesis. Three settings do not verify it.
    An unresolved signal receives no guessed exponent or derivative.
    """
    _records(records, 3)
    v0, v1, v2 = (row.volume for row in records)
    ratio = v1/v0
    if ratio <= 1 or v2*v0 != v1*v1:
        raise ValueError("three increasing geometrically spaced volumes required")
    energies = tuple(row.interval for row in records)
    first, second = energies[0]-energies[1], energies[1]-energies[2]
    if all(row.error == 0 for row in records) and first == second == point(0):
        return {"status": "zero_work_in_declared_family", "alpha": point(0),
                "amplitude": None, "offset": None, "pressure": point(0),
                "reference_ratio": point(0) if energies[0].lo > 0 else None}
    if max(value.lo for value in energies) <= min(value.hi for value in energies):
        # Even error intervals touching at one common value admit alpha=0.
        return {"status": "unresolved_energy_differences"}
    if first.hi <= 0 or second.hi <= 0:
        raise ValueError("records contradict the decreasing positive-amplitude family")
    if first.lo <= 0 or second.lo <= 0:
        return {"status": "unresolved_energy_differences"}
    contraction = second/first
    if contraction.lo >= 1:
        raise ValueError("records contradict a finite nonnegative power response")
    if contraction.hi >= 1:
        return {"status": "unresolved_contraction"}
    alpha = Interval(log_bounds(1/contraction.hi).lo,
                     log_bounds(1/contraction.lo).hi)/log_bounds(ratio)
    amplitude = first/(point(1)-contraction)
    offset = energies[0]-amplitude
    pressure = alpha*amplitude/point(v0)
    return {"status": "identified_within_declared_family", "contraction": contraction,
            "alpha": alpha, "amplitude": amplitude, "offset": offset,
            "pressure": pressure,
            "reference_ratio": pressure*point(v0)/energies[0]
                if energies[0].lo > 0 else None}


def supplied_power_records(amplitude, offset, contraction, *, volume=F(3),
                           volume_ratio=F(64), error=F(), state_tag="held-canonical-state"):
    """Synthetic calibration cases of the added family, never native output."""
    a, c, t, v, r = map(exact, (amplitude, offset, contraction, volume, volume_ratio))
    if a <= 0 or not 0 < t <= 1 or v <= 0 or r <= 1:
        raise ValueError("invalid positive-amplitude, nonnegative-power family")
    return tuple(EnergyRecord(v*r**j, c+a*t**j, error, state_tag) for j in range(3))


def pressure_stencil(records, *, derivative_order, derivative_bound):
    """Enclose -H'(V0) from V0+-h and a uniform bound on H'' or H'''."""
    _records(records, 2)
    if type(derivative_order) is not int or derivative_order not in (2, 3):
        raise ValueError("a second- or third-derivative bound is required")
    bound = exact(derivative_bound)
    if bound < 0:
        raise ValueError("derivative bound must be nonnegative")
    left, right = records
    h = (right.volume-left.volume)/2
    if h <= 0:
        raise ValueError("ordered distinct volumes required")
    center = (left.volume+right.volume)/2
    pressure = -(right.interval-left.interval)/point(2*h)
    remainder = bound*h/2 if derivative_order == 2 else bound*h*h/6
    return {"reference_volume": center, "half_step": h,
            "pressure": pressure+Interval(-remainder, remainder),
            "record_error_radius": (left.error+right.error)/(2*h),
            "regularity_radius": remainder,
            "energy_zero_required_for_pressure": False}


def two_setting_family(first, second, contraction):
    """Every allowed contraction gives the same two energy settings."""
    e0, e1, t = map(exact, (first, second, contraction))
    if not e0 > e1 or not 0 < t < 1:
        raise ValueError("a decreasing pair and contraction in (0,1) required")
    amplitude = (e0-e1)/(1-t)
    return amplitude, e0-amplitude


def invisible_strain_polynomial(reference_volume, sampled_volumes):
    """u(V0)=0, u'(V0)=1, and u vanishes at every sampled volume."""
    v0 = exact(reference_volume)
    nodes = tuple(exact(x) for x in sampled_volumes)
    if v0 <= 0 or any(x <= 0 for x in nodes):
        raise ValueError("positive volume settings required")
    coefficients = [-v0, F(1)]
    for node in sorted(set(nodes)-{v0}):
        updated = [F()]*(len(coefficients)+1)
        for degree, value in enumerate(coefficients):
            updated[degree] -= value*node/(v0-node)
            updated[degree+1] += value/(v0-node)
        coefficients = updated
    return tuple(coefficients)


def native_census_certificate():
    """Actual two-attempt E tapes realize the initial charge projectors."""
    paths = census.histories("E", prefix=0)
    states = census.native.states(binary=True)
    groups = defaultdict(list)
    for path in paths:
        if path.terminal != "success" or len(path.tape) != 2 or path.chance != F(1, 4):
            raise ArithmeticError("the complete two-attempt E census changed")
        groups[path.tape].append(path)
    totals = [[F()]*66 for _ in range(3)]
    rows = []
    for tape, branches in sorted(groups.items()):
        charge = sum(census.decode_core(tape, 0, "E"))
        # Each branch is a partial permutation. This verifies zero off-diagonal
        # effects, including all source/reference coherences, without sampling.
        if len({p.endpoint for p in branches}) != len(branches):
            raise ArithmeticError("native tape is not an injective branch")
        for path in branches:
            if charge != sum(states[path.initial][:2]):
                raise ArithmeticError("native records misidentify initial core charge")
            totals[charge][path.initial] += path.chance
        rows.append({"tape": tape, "decoded_initial_charge": charge,
                     "input_indices": tuple(p.initial for p in branches),
                     "endpoints": tuple(p.endpoint for p in branches),
                     "branch_weight": F(1, 4)})
    for q in range(3):
        if totals[q] != [F(sum(state[:2]) == q) for state in states]:
            raise ArithmeticError("complete charge effects differ from the source projectors")
    return {"native_attempts": 2, "law": "E", "source_dimension": 66,
            "raw_input_branches": len(paths), "complete_tapes": len(rows),
            "effect_matrix_units_checked": 3*66**2,
            "all_tapes": rows, "success_probability": F(1),
            "native_strain_or_energy_measurements": False,
            "unchanged_quantum_terminal_state_claimed": False}


def sector_work_interval(probabilities, slopes):
    """Average sector force enclosures; the probabilities require their own record law."""
    p = tuple(exact(x) for x in probabilities)
    if (len(p) != 3 or min(p) < 0 or sum(p) != 1 or type(slopes) is not tuple
            or len(slopes) != 3 or any(type(x) is not Interval for x in slopes)):
        raise ValueError("three normalized sector weights and force intervals required")
    return sum((point(weight)*slope for weight, slope in zip(p, slopes, strict=True)), point(0))


def native_scalar_replay():
    path = ROOT/"code/equation_of_state/native/verify_native_eos.py"
    spec = importlib.util.spec_from_file_location("operational_strain_native_verifier", path)
    verifier = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verifier)
    report = json.loads((path.parent/"native_eos_receipt.json").read_text())
    checked = verifier.verify(report)
    if checked["status"] != "PASS":
        raise ArithmeticError(checked["errors"])
    return checked, F(report["cases"][0]["initial_quadratic"])


def source_bindings():
    # The inherited list includes narrative documentation. This independent
    # certificate binds the executable/mathematical sources, never that prose.
    names = {name for name in census.v.SOURCE_PATHS if not name.startswith("docs/")}
    names.update({"Lean/Thermodynamics/EquationOfState.lean",
                  "code/equation_of_state/native/native_eos_receipt.json",
                  "code/equation_of_state/native/verify_native_eos.py",
                  "code/equation_of_state/maxwell/maxwell_eos.py",
                  "model/native_equilibration/history_entropy.py"})
    names.update(p.relative_to(ROOT).as_posix() for p in HERE.glob("*.py"))
    names.update(p.relative_to(ROOT).as_posix() for p in HERE.glob("*.md"))
    # Include function-local imports and eager package initializers. The
    # inherited manifest alone does not close those executable dependencies.
    pending = list(names)
    visited = set()

    def add_module(module):
        if not module.startswith(("model.", "tools.")):
            return
        stem = module.replace(".", "/")
        parts = stem.split("/")
        candidates = [stem+".py", stem+"/__init__.py"]
        candidates.extend("/".join(parts[:end])+"/__init__.py"
                          for end in range(1, len(parts)))
        for name in candidates:
            if (ROOT/name).is_file() and name not in names:
                names.add(name)
                pending.append(name)

    while pending:
        name = pending.pop()
        if name in visited or not name.endswith(".py"):
            continue
        visited.add(name)
        module = name.removesuffix(".py").replace("/", ".")
        package = (module.removesuffix(".__init__") if module.endswith(".__init__")
                   else module.rpartition(".")[0])
        add_module(module)
        for node in ast.walk(ast.parse((ROOT/name).read_text())):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    add_module(alias.name)
            elif isinstance(node, ast.ImportFrom):
                target = "."*node.level+(node.module or "")
                target = importlib.util.resolve_name(target, package) if node.level else target
                add_module(target)
                for alias in node.names:
                    add_module(target+"."+alias.name)
    if any(name.startswith("docs/") for name in names):
        raise ArithmeticError("documentation cannot be a mathematical source binding")
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(names)}


def encode(value):
    if isinstance(value, F):
        return str(value)
    if type(value) is Interval:
        return [str(value.lo), str(value.hi)]
    if type(value) is EnergyRecord:
        return {name: encode(getattr(value, name)) for name in ("volume", "energy", "error", "state_tag")}
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode(item) for item in value]
    return value


def run():
    before = source_bindings()
    native_check, q = native_scalar_replay()
    examples = []
    for contraction, alpha in ((F(1, 4), F(1, 3)), (F(1, 8), F(1, 2))):
        records = supplied_power_records(q, q, contraction, error=F(1, 2**20))
        decoded = identify_power_response(records)
        if not decoded["alpha"].contains(alpha) or not decoded["pressure"].contains(alpha*q/3):
            raise ArithmeticError("power-response enclosure failed")
        examples.append({"explicitly_supplied_alpha": alpha, "records": records,
                         "decoded_only_from_records": decoded})
    if examples[0]["decoded_only_from_records"]["alpha"].hi >= examples[1]["decoded_only_from_records"]["alpha"].lo:
        raise ArithmeticError("the two supplied strain laws are unresolved")
    sector_records = tuple(supplied_power_records(a, c, t, error=F(1, 2**24),
                           state_tag=f"charge-sector-{i}") for i, (a, c, t) in
                           enumerate(((2, 7, F(1, 2)), (3, 11, F(1, 4)), (5, 13, F(1, 8)))))
    sector_slopes = tuple(identify_power_response(rows)["pressure"] for rows in sector_records)
    weights = census.v.common.reference()
    mean, pair = census.v.moments(weights)
    probabilities = (1-mean+pair, mean-2*pair, pair)
    average = sector_work_interval(probabilities, sector_slopes)
    exact_average = sum((p*a*alpha/3 for p, a, alpha in zip(probabilities, (2, 3, 5),
                         (F(1, 6), F(1, 3), F(1, 2)), strict=True)), F())
    if not average.contains(exact_average):
        raise ArithmeticError("native sector force averaging failed")
    polynomial = invisible_strain_polynomial(3, (F(2), F(3), F(4), F(192)))
    result = {"schema": "oph.operational_strain.v1", "source_sha256": before,
        "native_scalar_source": native_check, "native_Q_used_only_as_example_amplitude": q,
        "supplied_strain_examples": examples,
        "zero_work": identify_power_response(supplied_power_records(72, 72, 1)),
        "two_setting_ambiguity": [{"contraction": t,
            "amplitude_offset": two_setting_family(144, 90, t)} for t in (F(1, 4), F(1, 8))],
        "native_charge_readout": native_census_certificate(),
        "supplied_sector_calibration": {"records": sector_records,
            "native_reference_charge_probabilities": probabilities,
            "pressure_interval": average, "example_exact_pressure": exact_average,
            "work_law_identified_as_physical": False},
        "general_regular_response": pressure_stencil((EnergyRecord(F(29, 10), F(841, 100), F(1, 10000)),
            EnergyRecord(F(31, 10), F(961, 100), F(1, 10000))), derivative_order=2, derivative_bound=2),
        "invisible_positive_response": {"u_coefficients_low_to_high": polynomial,
            "family": "H_d(V)=B*exp(d*u(V)/B), B>0; H_d'(3)=d",
            "same_positive_samples_for_every_real_d": True},
        "scope": {"finite_family_derivative_identified_from_added_records": True,
            "absolute_energy_zero_inferred_within_nondegenerate_power_family": True,
            "native_charge_census_used_without_postselection": True,
            "arbitrary_coherent_source_and_reference_effect_identity": True,
            "physical_volume_and_strain_intervention_source_selected": False,
            "energy_meter_or_work_coupling_constructed_from_native_repair": False,
            "generic_response_family_verified_by_three_records": False,
            "equilibrium_or_physical_equation_of_state_derived": False}}
    if source_bindings() != before:
        raise ArithmeticError("source files changed during certification")
    return encode(result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.write:
        OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
    if args.check and json.loads(OUTPUT.read_text()) != result:
        raise SystemExit("operational strain receipt differs from exact replay")
    print("operational strain: exact finite records and native census verified")


if __name__ == "__main__":
    main()
