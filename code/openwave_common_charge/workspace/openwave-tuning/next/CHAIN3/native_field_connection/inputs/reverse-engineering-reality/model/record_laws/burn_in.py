"""Uniform forgetting for every admitted native current-record law.

The comparison state is the uniform binary density with the original,
untouched reference marginal. Warm-up records and controller memory are
discarded; they are not asserted to decouple. See README for the proof.
All retained bounds and structural checks use exact arithmetic.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from pathlib import Path

from tools import strict_json
from model.observer import native, quantum_mixing as mixing
from . import family as f

ROOT = Path(__file__).resolve().parents[2]
D = 66
S = f.S
ZERO = (0,) * 12
R_UPPER = F(19079, 20000)  # (55+sqrt(5))/60 < this rational.
C0_RATE = S(F(43, 60), F(2, 15))
SOURCE_PATHS = (
    "evidence/native_geometric_source_20260925/spec.json",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/__init__.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/__init__.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/claims.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/constants/__init__.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/constants/oph_pixel.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/core/__init__.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/core/array_geometry.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/core/echosahedral_dynamics.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/core/echosahedral_federation.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/core/graph.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/core/icosahedral.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/core/patch_state.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/core/patchnet.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/core/pixel_scale.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/core/screen_microphysics.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/core/screen_ports.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/__init__.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/anneal.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/consensus_certificate.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/kernel_dispatcher.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/modular.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/positive_geometry.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/repair.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/evidence/__init__.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/evidence/bundle.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/evidence/cross_repo_artifacts.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/evidence/hashes.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/evidence/particle_input_policy.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/evidence/production_envelope.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/evidence/verifier.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/finite_groups.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/gauge/__init__.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/gauge/a5_sm_certificate.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/gauge/covariant_overlap.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/gauge/higgs_carrier.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/gauge/mar_sieve.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/gauge/physical_a5_sm_requirements.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/gauge/repair_projection.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/gauge/source_to_eft.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/gauge/yang_mills_gap.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/groups/__init__.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/groups/base.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/groups/clock.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/groups/registry.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/groups/s3.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/groups/z2.py",
    "model/__init__.py",
    "model/common_source/__init__.py",
    "model/common_source/experiment.py",
    "model/common_source/gauge_matter.py",
    "model/common_source/native_join.py",
    "model/common_source/records.py",
    "model/common_source/scalar_coupling.py",
    "model/common_source/slow_geometry.py",
    "model/observer/__init__.py",
    "model/observer/native.py",
    "model/observer/p_response/__init__.py",
    "model/observer/p_response/experiment.py",
    "model/observer/quantum_mixing.py",
    "model/observer/scalar_response/__init__.py",
    "model/observer/scalar_response/experiment.py",
    "model/observer/source_discrimination/certificate.py",
    "model/observer/source_discrimination/circuit.py",
    "model/observer/source_scalar/__init__.py",
    "model/observer/source_scalar/bounds.py",
    "model/observer/source_scalar/compiler.py",
    "model/observer/source_scalar/cutoff.py",
    "model/observer/source_scalar/experiment.py",
    "model/observer/source_scalar/field.py",
    "model/predictions/__init__.py",
    "model/predictions/observer_instrument.py",
    "model/record_laws/__init__.py",
    "model/record_laws/burn_in.py",
    "model/record_laws/family.py",
    "tools/strict_json.py",
)


def _count(value):
    if type(value) is not int or value < 0:
        raise ValueError("nonnegative integer attempt count required")
    return value


def _coverage(value):
    if type(value) not in (int, F) or not 0 < value <= 1:
        raise ValueError("exact scheduler coverage in (0,1] required")
    return F(value)


def validate_schedule(weights):
    """Exact exogenous seam probabilities; no claim to infer independence."""
    weights = tuple(weights)
    if (len(weights) != len(mixing.seams())
            or any(type(w) not in (int, F) or w < 0 for w in weights)
            or sum(weights) != 1):
        raise ValueError("thirty exact nonnegative seam probabilities summing to one required")
    return tuple(map(F, weights))


def coverage_of(weights):
    """Lower envelope alpha such that every seam has probability >=alpha/30."""
    return 30 * min(validate_schedule(weights))


def binary_bound(attempts, coverage=1):
    """Half trace distance, uniformly over untouched references and admitted laws."""
    n, alpha = _count(attempts), _coverage(coverage)
    rate = 1-alpha*(1-R_UPPER)
    return min(F(1), F(D, 2)*rate**n)


def full_bound(attempts, split=None, coverage=1):
    """Full 78-state bound using an even collision-removal split.

    The split is an analytic proof device, not an executed measurement/reset.
    None chooses the best even split of this bound, not a physical mixing time.
    """
    n, alpha = _count(attempts), _coverage(coverage)
    if split is not None and (type(split) is not int or split % 2 or not 0 <= split <= n):
        raise ValueError("an even split between zero and the attempt count required")
    mu = 1-alpha/6
    rate = 1-alpha*(1-R_UPPER)
    splits = range(0, n+1, 2) if split is None else (split,)
    value, selected = min((mu**(s//2)+mu**s+F(D, 2)*rate**(n-s), s)
                          for s in splits)
    return {"error": min(F(1), value), "split": selected}


def sufficient_attempts(tolerance, *, full=False, coverage=1):
    """First count certified by these monotone bounds, not a sharp mixing time."""
    if type(tolerance) not in (int, F) or not 0 < tolerance < 1 or type(full) is not bool:
        raise ValueError("an exact tolerance in (0,1) and a Boolean full-carrier switch required")
    alpha = _coverage(coverage)
    def error(n):
        return full_bound(n, coverage=alpha)["error"] if full else binary_bound(n, alpha)
    lower, upper = 0, 1
    while error(upper) > tolerance:
        lower, upper = upper, 2*upper
    while lower+1 < upper:
        middle = (lower+upper)//2
        if error(middle) <= tolerance:
            upper = middle
        else:
            lower = middle
    return upper


@lru_cache(maxsize=None)
def _images(seam_index):
    """Current probabilities and outputs from actual native coin invocations."""
    seam = mixing.seams()[seam_index]
    result = []
    for state in f.STATES:
        counts = Counter()
        for coin in (0, 1):
            output, (_, current) = native.native(state, seam, coin)
            counts[current, f.INDEX[output]] += 1
        result.append({current: (output, F(count, 2))
                       for (current, output), count in counts.items()})
    return tuple(result)


def _root_probability(value):
    return {F(1): S(1), F(1, 2): f.ROOT_HALF, F(1, 4): S(F(1, 2))}[value]


def dominant_column(x, y):
    """Uniform C0 column, independently assembled from native probabilities."""
    rational, irrational = Counter(), Counter()
    for edge in range(30):
        left, right = _images(edge)[x], _images(edge)[y]
        for current in left.keys() & right.keys():
            a, p = left[current]
            b, q = right[current]
            product = p*q
            if product == F(1, 2):
                irrational[a, b] += 1
            else:
                rational[a, b] += {F(1): 2, F(1, 4): 1}[product]
    return {key: S(F(rational[key], 60), F(irrational[key], 60))
            for key in rational.keys() | irrational.keys()}


def dominant_reduction(x, y):
    """Independent induced-exclusion expression on the binary sector."""
    sx, sy = f.STATES[x], f.STATES[y]
    if max(sx) > 1 or max(sy) > 1:
        raise ValueError("binary configurations required")
    delta = tuple(a-b for a, b in zip(sx, sy, strict=True))
    cut1 = sum(abs(delta[a]-delta[b]) == 1 for a, b in mixing.seams())
    cut2 = sum(abs(delta[a]-delta[b]) == 2 for a, b in mixing.seams())
    survival = 1-((1-f.ROOT_HALF)*cut1+S(F(1, 2))*cut2)*F(1, 30)
    out = defaultdict(S, {(x, y): survival})
    for a, b in mixing.seams():
        if delta[a] or delta[b] or sx[a] == sx[b]:
            continue
        left, right = list(sx), list(sy)
        left[a], left[b] = left[b], left[a]
        right[a], right[b] = right[b], right[a]
        out[x, y] -= F(1, 60)
        out[f.INDEX[tuple(left)], f.INDEX[tuple(right)]] += F(1, 60)
    return {pair: coefficient for pair, coefficient in out.items() if coefficient}


def structural_audit():
    """Every full matrix unit, exact binary reduction and transient structure."""
    if not 2236**2 < 5*1000**2 < 2237**2 or R_UPPER != (55+F(2237, 1000))/60:
        raise ArithmeticError("the rational diagonal-rate upper enclosure failed")
    import sympy as sp
    laplacian = sp.zeros(12)
    for a, b in mixing.seams():
        laplacian[a, a] += 1
        laplacian[b, b] += 1
        laplacian[a, b] -= 1
        laplacian[b, a] -= 1
    variable = sp.Symbol("x")
    polynomial = laplacian.charpoly(variable).as_expr()
    if sp.expand(polynomial-variable*(variable-6)**5*(variable**2-10*variable+20)**3) != 0:
        raise ArithmeticError("the native carrier Laplacian gap changed")
    binary = {i for i, state in enumerate(f.STATES) if max(state) == 1}
    columns, rates, digest = {}, {}, hashlib.sha256()
    transient_max = S()
    for x in range(78):
        for y in range(78):
            column = dominant_column(x, y)
            delta = tuple(a-b for a, b in zip(f.STATES[x], f.STATES[y], strict=True))
            for (u, v), coefficient in column.items():
                if tuple(a-b for a, b in zip(f.STATES[u], f.STATES[v], strict=True)) != delta:
                    raise ArithmeticError("dominant channel changed signed difference")
                if coefficient.sign() < 0:
                    raise ArithmeticError("dominant column is not nonnegative")
            if x in binary and y in binary:
                if column != dominant_reduction(x, y):
                    raise ArithmeticError("native C0 differs from the exclusion reduction")
                columns[x, y] = column
                mass = sum(column.values(), S())
                if delta in rates and rates[delta] != mass:
                    raise ArithmeticError("dominant sector column mass is not constant")
                rates[delta] = mass
            else:
                for (u, v), coefficient in column.items():
                    if (u not in binary or v not in binary) and (u, v) != (x, y):
                        raise ArithmeticError("collision block is not diagonal")
                diagonal = column.get((x, y), S())
                if (diagonal-transient_max).sign() > 0:
                    transient_max = diagonal
            digest.update(json.dumps([x, y, [[u, v, str(c)]
                            for (u, v), c in sorted(column.items())]],
                            separators=(",", ":")).encode())
    for pair, column in columns.items():
        for output, coefficient in column.items():
            if columns[output].get(pair, S()) != coefficient:
                raise ArithmeticError("binary C0 comparison is not symmetric")
    maximum = S()
    for delta, rate in rates.items():
        if delta != ZERO and (rate-maximum).sign() > 0:
            maximum = rate
    if maximum != C0_RATE or any((maximum-r).sign() < 0 for d, r in rates.items() if d != ZERO):
        raise ArithmeticError("sharp dominant coherence rate changed")
    if transient_max != F(5, 6) or rates[ZERO] != 1:
        raise ArithmeticError("collision survival or stochastic sector changed")
    # Strict exact comparison to the known diagonal contraction, using sqrt5>2.
    if not (S(F(57, 60))-maximum).sign() > 0:
        raise ArithmeticError("coherence is not faster than diagonal mixing")
    return {"full_matrix_units": 6084, "binary_matrix_units": 4356,
            "dominant_channel_sha256": digest.hexdigest(),
            "exact_binary_exclusion_reduction": True, "binary_symmetric": True,
            "maximum_offdiagonal_row_and_column_sum": str(maximum),
            "attaining_signed_difference_sectors": sum(r == maximum for d, r in rates.items() if d != ZERO),
            "collision_block_diagonal": True, "collision_survival": str(transient_max),
            "rational_diagonal_rate_upper": str(R_UPPER),
            "root5_squared_integer_bounds": [2236**2, 5*1000**2, 2237**2],
            "carrier_characteristic_polynomial": "x*(x-6)^5*(x^2-10*x+20)^3",
            "fixed_uniform_schedule_gap_all_tables": "(5-sqrt(5))/60"}


def examples():
    """Admissible points including complex phases and non-extremal laws."""
    return {"C0": f.table(1, f.ROOT_HALF, f.ROOT_HALF, F(1, 2), f.ROOT_HALF),
            "C": f.table(1, F(1, 2), F(1, 2), 0, F(1, 2)),
            "E": f.table(1, 0, 0, 0, 0), "T": f.table(0, 0, 0, 0, 0),
            "complex": f.table(S(0, 0, -1), S(0, 0, 0, F(1, 2)),
                               -f.ROOT_HALF, F(1, 2), S(0, 0, 0, F(1, 2))),
            "quarter": f.generic_tables()["C0_quarter_turn"],
            "mixture": f.generic_tables()["C_T_mixture"]}


def coefficient_audit():
    checked = 0
    for table in examples().values():
        if not f.positive_semidefinite(table):
            raise ArithmeticError("an example is outside the CP family")
        for (current, left, right), coefficient in f.multipliers(table).items():
            bound = f.PROBABILITY[left, current]*f.PROBABILITY[right, current]
            if (S(bound)-coefficient*coefficient.conj()).sign() < 0:
                raise ArithmeticError("two-by-two CP domination failed")
            checked += 1
    return {"admissible_examples": len(examples()), "coefficient_checks": checked,
            "universal_proof": "each positive G_j has |G_j(a,b)|^2 <= P(j|a)P(j|b)"}


def branch(table, edge, current, rho):
    """One coarse-current branch on sparse (source,source,reference,reference)."""
    a, b = mixing.seams()[edge]
    multipliers = f.multipliers(table)
    images = _images(edge)
    out = defaultdict(S)
    for (x, y, r, s), value in rho.items():
        if current not in images[x] or current not in images[y]:
            continue
        kx, ky = (f.STATES[x][a], f.STATES[x][b]), (f.STATES[y][a], f.STATES[y][b])
        coefficient = multipliers.get((current, kx, ky), S())
        if coefficient:
            u, v = images[x][current][0], images[y][current][0]
            out[u, v, r, s] += S.of(value)*coefficient
    return {key: value for key, value in out.items() if value}


def execute(rho, policy, schedules):
    """Finite allowed law-only feedback; no scheduler or source feedback.

    policy(prefix, edge) returns a CP table before the current is produced.
    Schedules are supplied exogenously. Every nonzero current branch remains.
    Returned warm-up records must be discarded to use the burn-in theorem.
    """
    histories = {(): {key: S.of(value) for key, value in rho.items() if value}}
    checked = set()
    for schedule in schedules:
        weights = validate_schedule(schedule)
        following = {}
        for prefix, state in histories.items():
            for edge, weight in enumerate(weights):
                if not weight:
                    continue
                table = policy(prefix, edge)
                signature = tuple(table[key] for key in f.PARAMETERS)
                if signature not in checked:
                    if not f.positive_semidefinite(table):
                        raise ValueError("feedback selected a non-CP table")
                    checked.add(signature)
                for current in (-1, 0, 1):
                    output = branch(table, edge, current, state)
                    if output:
                        following[prefix+((edge, current),)] = {k: weight*v for k, v in output.items()}
        histories = following
    return histories


def forget(histories):
    output = defaultdict(S)
    for state in histories.values():
        for key, value in state.items():
            output[key] += value
    return {key: value for key, value in output.items() if value}


def reference_marginal(rho):
    output = defaultdict(S)
    for (x, y, r, s), value in rho.items():
        if x == y:
            output[r, s] += value
    return {key: value for key, value in output.items() if value}


def adaptive_input():
    """An imaginary Bell coherence involving a collision and a binary state."""
    x, y = f.representative((2, 0)), f.representative((1, 1))
    return {(x, x, 0, 0): S(F(1, 2)), (y, y, 1, 1): S(F(1, 2)),
            (x, y, 0, 1): S(0, 0, F(-1, 2)),
            (y, x, 1, 0): S(0, 0, F(1, 2))}


def adaptive_policy(prefix, edge):
    return examples()["complex" if prefix and prefix[-1][1] > 0 else "C0"]


def adaptive_schedules():
    return ((F(1, 30),)*30, (F(31, 60),)+(F(1, 60),)*29)


@lru_cache(maxsize=1)
def adaptive_outputs():
    return execute(adaptive_input(), adaptive_policy, adaptive_schedules())


def adaptive_witness():
    histories = adaptive_outputs()
    reference = reference_marginal(forget(histories))
    if reference != reference_marginal(adaptive_input()):
        raise ArithmeticError("law-only warm-up changed an untouched reference")
    rows = [[list(prefix), [[*key, str(value)] for key, value in sorted(state.items())]]
            for prefix, state in sorted(histories.items())]
    return {"attempts": 2, "complete_record_blocks": len(histories),
            "complete_output_sha256": hashlib.sha256(json.dumps(rows, separators=(",", ":")).encode()).hexdigest(),
            "reference_marginal": [[*key, str(value)] for key, value in sorted(reference.items())],
            "initial_collision_binary_imaginary_entanglement": True,
            "law_selected_from_preceding_current": True,
            "exogenous_schedule_coverage": str(min(map(coverage_of, adaptive_schedules()))),
            "warmup_performed_by_new_source_pulse_words": False}


def controls():
    """Exact exclusions: retained tape and record-dependent seam choice."""
    binary = [i for i, x in enumerate(f.STATES) if max(x) == 1]
    tau = {(i, i, 0, 0): S(F(1, D)) for i in binary}
    edge = mixing.seams().index((0, 1))
    alternative = mixing.seams().index((0, 5))
    c = examples()["C"]
    outputs = {j: branch(c, edge, j, tau) for j in (-1, 0, 1)}
    if forget(outputs) != tau:
        raise ArithmeticError("a seam channel failed to preserve the binary uniform density")
    plus = sum(reference_marginal(outputs[1]).values(), S())
    occupied = sum((v for (x, y, _, _), v in outputs[1].items()
                    if x == y and f.STATES[x][0]), S())
    product = plus * F(1, 6)
    biased = {}
    for first, state in outputs.items():
        selected = alternative if first > 0 else edge
        for second in (-1, 0, 1):
            biased[first, second] = branch(c, selected, second, state)
    final = forget(biased)
    tv = sum((abs((final.get(key, S())-value).a) for key, value in tau.items()), F(0))/2
    if (plus, occupied, product) != (F(5, 66), 0, F(5, 396)) or tv <= 0:
        raise ArithmeticError("a premise-control witness disappeared")
    return {"last_current_plus_probability": str(plus),
            "last_current_plus_and_endpoint_occupied": str(occupied),
            "incorrect_source_tape_product_prediction": str(product),
            "record_selected_seam_changes_stationary_population_TV": str(tv),
            "controller_and_tape_decoupling_claimed": False}


def run():
    structural = structural_audit()
    thresholds = {}
    for target in (F(1, 100), F(1, 1000)):
        # Fixed finite certificates; no data or fitted stopping criterion.
        n = sufficient_attempts(target)
        full_n = sufficient_attempts(target, full=True)
        full = full_bound(full_n)
        if full["error"] > target:
            raise ArithmeticError("declared full-carrier warm-up budget failed")
        def rounded_upper(value):
            scale = 10**12
            scaled = value*scale
            return F(-(-scaled.numerator//scaled.denominator), scale)
        thresholds[str(target)] = {"binary_attempts": n,
                                  "binary_error_upper": str(rounded_upper(binary_bound(n))),
                                  "full_attempts": full_n, "full_even_split": full["split"],
                                  "full_error_upper": str(rounded_upper(full["error"]))}
    return {"schema": "oph.record_laws.burn_in.v1",
            "source_sha256": {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                              for name in SOURCE_PATHS},
            "structural": structural, "coefficient_domination": coefficient_audit(),
            "finite_budgets": thresholds, "controls": controls(),
            "adaptive_reference_execution": adaptive_witness(),
            "contract": {
                "states": "arbitrary total-two source density and untouched quantum reference",
                "law": "any CP table, varying with seam, attempt and preceding coarse currents",
                "schedule": "exogenous seam distributions, each probability >= alpha/30; alpha>0",
                "classical_rule": "unchanged fair native placement probabilities at every call",
                "warmup_readout": "discard the complete warm-up tape and controller at the boundary",
                "clock": "repair attempt count, including zero current; no duration identification",
                "later_protocol": "one common complete CPTP instrument, with every timeout retained",
                "reference_distance": "half trace norm to tau_binary tensor original reference marginal",
                "binary_bound": "min(1,33*(1-alpha*(1-19079/20000))^n)",
                "full_bound": "min_even_s min(1,(1-alpha/6)^(s/2)+(1-alpha/6)^s+33*r_alpha^(n-s))",
                "exclusion_gap_import": "Caputo, Liggett, Richthammer (2010), section 4.1.1, equation (4.1)"},
            "scope": {"exact_initial_stationary_preparation_required": False,
                      "same_native_conserved_sector_required": True,
                      "arbitrary_record_basis_or_coherent_environment_feedback": False,
                      "record_dependent_scheduler_covered": False,
                      "physical_record_law_selected": False,
                      "physical_clock_or_interacting_QFT_derived": False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--check", type=Path)
    modes.add_argument("--write", type=Path)
    args = parser.parse_args()
    result = run()
    payload = json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+"\n"
    if args.check and json.dumps(strict_json.load(args.check), sort_keys=True) != json.dumps(result, sort_keys=True):
        raise ValueError("law-uniform burn-in receipt differs")
    if args.write:
        with args.write.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(payload)
    print(json.dumps({"exact_replay": True, "all_tables_same_gap": result["structural"]["fixed_uniform_schedule_gap_all_tables"],
                      "physical_law_selected": False}, sort_keys=True))


if __name__ == "__main__":
    main()
