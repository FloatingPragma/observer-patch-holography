"""A common reference prepared with two finite, incident-clocked phases.

Warm up on CORE incident ticks, then herald movement at the NEXT specified
CORE tick. Complete on the mixed complement schedule, counting CAP0 incident
ticks. Both waiting loops have a finite private attempt deadline. The public
time is local incident time; the private deadline counter is still a controller
service. Native auxiliary coin instruments supply all seam choices. Their
private current and quantum output, and the preparation's fine tape, are
discarded in the comparison; all four terminal flags and work states remain.

The stopped-channel argument also applies to a conditional expectation E:
each seam fixes and preserves E, and on ker(E), G<= (1-gamma) Id. Writing
G=A+R, p=Tr(A(rho)), M=Id-R gives p Id<=M<=Id. The next-tick map is A M^-1,
similar to B=M^-1/2 A M^-1/2. Since 0<=B<=Id-gamma M^-1<= (1-gamma) Id,
||S^n-E||_2 <= p^-1/2 (1-gamma)^n. Source dimension66 then gives the
arbitrary-reference half trace bound33 p^-1/2 (1-gamma)^n. This does not
assert self-adjointness of S. Here E=Dep in phase1 and E_CORE in phase2.

Actual finite ticks use killed geometric sums, not this infinite-wait map.
If A,R include accepted selector probabilities and tau is selector timeout,
success is A sum_{j<h} R^j; timeout is tau sum_{j<h} R^j+R^h. Mixed-cap
selector acceptance differs between its12- and15-address lists and is not
pulled out as a common scalar. Coupling to the ideal schedule costs at most
h*tau+(1-p)^h per tick. All rational error certificates include that loss.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix

from tools import strict_json

from ..native_clock import regional_preparation as incident
from ..native_equilibration import cut_completion as cut
from ..observer import native, quantum_mixing
from ..repair_scheduler import coins
from . import common_caps

CORE, CAPS = common_caps.CORE, common_caps.CAPS
ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name("local_preparation_results.json")
SOURCE_PATHS = tuple(sorted(set((*incident.SOURCE_PATHS,
    *common_caps.preparation.events.SOURCE_PATHS,
    *common_caps.preparation.relational.SOURCE_PATHS,
    *common_caps.regional.transfer.SOURCE_PATHS,
    "model/native_modular_reference/local_preparation.py",
    "model/native_modular_reference/common_caps.py",
    "model/native_modular_reference/regional.py",
    "model/native_modular_reference/preparation.py",
    "model/repair_scheduler/coins.py"))))
FLAGS = ("success", "herald_failure", "warmup_timeout", "completion_timeout")


def _integer(value, low, high):
    if type(value) is not int or not low <= value <= high:
        raise ValueError(f"an integer in [{low}, {high}] required")
    return value


def _dyadic_upper(value, bits=80):
    scale = 2**bits
    return F((value.numerator*scale+value.denominator-1)//value.denominator, scale)


@lru_cache(maxsize=1)
def _reference_payload():
    """Own immutable scalars after the source derivation has been checked."""
    raw, completed, q = incident.reference(CORE)
    expected = (F(2, 225), F(1, 36), F(2, 45))
    basis = native.states(binary=True)
    if (q != F(5, 33) or sum(raw) != 1 or sum(completed) != 1
            or any(w != expected[sum(x[p] for p in CORE)]
                   for x, w in zip(basis, completed, strict=True))):
        raise ArithmeticError("source-derived common local reference changed")
    return tuple((w.numerator, w.denominator) for w in completed)


def reference():
    """The exact E_CORE reference, with no mutable cached Fraction aliases."""
    return tuple(F(n, d) for n, d in _reference_payload())


@lru_cache(maxsize=1)
def _schedule_payload():
    core, weights, cert = common_caps.cap_schedule()
    if core != CORE:
        raise ArithmeticError("mixed cap schedule lost its protected core")
    p = sum((w for e, w in weights.items() if set(e).intersection(CAPS[0])), F())
    if p != F(5, 24) or F(cert["gap_lower_per_attempt"]) != F(1, 600):
        raise ArithmeticError("source-derived mixed incident fraction or gap changed")
    gamma = F(cert["gap_lower_per_attempt"])
    return (tuple((edge, w.numerator, w.denominator) for edge, w in sorted(weights.items())),
            (p.numerator, p.denominator), (gamma.numerator, gamma.denominator))


def _schedule():
    """Detach cached values and check the actual declared native seam menu."""
    encoded, fraction, gap = _schedule_payload()
    weights = {edge: F(n, d) for edge, n, d in encoded}
    p, gamma = F(*fraction), F(*gap)
    expected = {}
    seams = quantum_mixing.seams()
    for cap in CAPS:
        available = tuple(edge for edge in seams if not set(edge).intersection(cap))
        for edge in available:
            expected[edge] = expected.get(edge, F())+F(1, len(CAPS)*len(available))
    incident_fraction = sum((w for edge, w in weights.items() if set(edge).intersection(CAPS[0])), F())
    if (weights != expected or p != incident_fraction or p != F(5, 24)
            or gamma != F(1, 600)):
        raise ArithmeticError("source-derived mixed schedule premises changed")
    return weights, p, gamma


def selectors(rounds=32):
    """Exact reduced native-coin selector laws and worst-case call caps."""
    _integer(rounds, 1, 32)
    seams = quantum_mixing.seams()
    warm = coins.selector_certificate(30, rounds)
    lists = tuple(tuple(e for e in seams if not set(e).intersection(cap)) for cap in CAPS)
    rows = tuple(coins.selector_certificate(len(es), rounds) for es in lists)
    accepted = {}
    for es, row in zip(lists, rows, strict=True):
        weight = F(row["accepted_probability_per_address"])/2
        for seam in es:
            accepted[seam] = accepted.get(seam, F())+weight
    timeout = sum((F(r["timeout_probability"])/2 for r in rows), F())
    if sum(accepted.values(), F())+timeout != 1:
        raise ArithmeticError("mixed selector omitted a timeout or an address")
    return {"warmup": warm, "completion_cap_coin_attempts": 1,
            "completion_address_selectors": rows,
            "completion_timeout_probability": timeout,
            "completion_accepted_seam_weights": accepted,
            "completion_auxiliary_calls_per_slot_upper":
                1+max(r["maximum_auxiliary_native_attempts"] for r in rows)}


def budget(precision_bits=24, deadline=256, selector_rounds=32):
    """Exact certificates; counts refer to attempts, not seconds.

    B bounds the complete four-flag channel against ideal success q*sigma,
    ideal no-movement (I/66-q*raw), and two zero timeout blocks, tensored with
    the original reference marginal. Consequently success probability is in
    q+-B, conditional success distance is <=2B/(q-B), and the successful
    work marginal has normalized sandwich radius derived from B/(q*minσ).
    """
    _integer(precision_bits, 12, 40)
    _integer(deadline, 1, 4096)
    selection = selectors(selector_rounds)
    _, p_completion, gamma = _schedule()
    p_warmup = F(len(incident.incident_edges(CORE)), 30)
    target = F(1, 2**precision_bits)
    blocks_warm = blocks_complete = 0
    while F(33)/p_warmup/2**blocks_warm > target:
        blocks_warm += 1
    while F(33)/p_completion/2**blocks_complete > target:
        blocks_complete += 1
    block_complete = math.ceil(1/gamma)
    warm_ticks, complete_ticks = 16*blocks_warm, block_complete*blocks_complete
    delta = incident.local_bound(warm_ticks, CORE)
    eta = F(33)/p_completion/2**blocks_complete
    warm_slots, complete_slots = (warm_ticks+1)*deadline, complete_ticks*deadline
    tau_warm = F(selection["warmup"]["timeout_probability"])
    tau_complete = selection["completion_timeout_probability"]
    warm_cutoff = _dyadic_upper((warm_ticks+1)*(deadline*tau_warm+(1-p_warmup)**deadline))
    complete_cutoff = _dyadic_upper(complete_ticks*(deadline*tau_complete+(1-p_completion)**deadline))
    q = F(5, 33)
    bound = delta+warm_cutoff+q*(eta+complete_cutoff)
    radius = bound/(q*min(reference()))
    if bound >= q or radius >= 1:
        raise ValueError("finite deadlines/selectors do not certify faithful preparation")
    auxiliary = warm_slots*selection["warmup"]["maximum_auxiliary_native_attempts"]
    auxiliary += complete_slots*selection["completion_auxiliary_calls_per_slot_upper"]
    # Couple to unbounded rejection and remote waits. A finite abort only
    # shortens that path. Slot costs have the same conditional mean after
    # every past reduced history by the native coarse-coin identity, so the
    # tail-sum argument applies even with an initially correlated auxiliary.
    expected_warm_slots = F(warm_ticks+1)/p_warmup
    expected_complete_slots = (q+bound)*complete_ticks/p_completion
    expected_work = expected_warm_slots+expected_complete_slots
    expected_auxiliary = F(16, 3)*expected_warm_slots+F(29, 5)*expected_complete_slots
    # One explicit reusable counter/register layout. This is a controller
    # subtotal, not a bound on program/lookup tables, interpreter working
    # memory, fine tapes, or the two quantum carrier registers.
    controller_registers = {
        "phase_and_terminal_flag": 3,
        "local_tick_counter": max(warm_ticks, complete_ticks).bit_length(),
        "within_tick_attempt_counter": deadline.bit_length(),
        "selector_round_counter": selector_rounds.bit_length(),
        "selector_bit_position": 3, "selector_word_accumulator": 5,
        "mixed_cap_choice": 1, "native_seam_address": 5,
        "current_or_contrast_reply": 2, "placement_coin_reply": 1,
        "fixed_law_identifier": 1,
    }
    return {
        "core": CORE, "completion_clock_cap": CAPS[0],
        "warmup_incident_fraction": p_warmup, "completion_incident_fraction": p_completion,
        "warmup_ticks": warm_ticks, "specified_herald_ticks": 1,
        "completion_ticks": complete_ticks, "per_tick_attempt_deadline": deadline,
        "selector_rounds": selector_rounds, "warmup_halftrace_upper": delta,
        "completion_halftrace_upper": eta, "warmup_cutoff_probability_upper": warm_cutoff,
        "completion_cutoff_probability_upper": complete_cutoff,
        "warmup_selector_timeout_per_slot": tau_warm,
        "completion_selector_timeout_per_slot": tau_complete,
        "uniform_ideal_herald_probability": q, "success_probability_interval": (q-bound, q+bound),
        "complete_flagged_reference_halftrace_upper": bound,
        "conditional_reference_uniform_halftrace_upper": min(F(1), 2*bound/(q-bound)),
        "operator_sandwich_radius": radius,
        "prepared_modular_generator_operator_error_upper": 2*radius/(1-radius),
        "work_native_attempts_upper": warm_slots+complete_slots,
        "auxiliary_native_attempts_upper": auxiliary,
        "total_native_attempts_upper": warm_slots+complete_slots+auxiliary,
        "expected_work_native_attempts_upper": expected_work,
        "expected_auxiliary_native_attempts_upper": expected_auxiliary,
        "expected_total_native_attempts_upper": expected_work+expected_auxiliary,
        "controller_register_bits": controller_registers,
        "controller_register_bits_subtotal": sum(controller_registers.values()),
        "controller_program_tables_interpreter_and_private_tapes_in_memory_subtotal": False,
        "private_per_tick_attempt_deadline_counter_required": True,
        "public_preparation_time_counts_only_incident_ticks": True,
        "initial_work_sector": "binary two-occupation66",
        "auxiliary_dimension": 78, "auxiliary_reset_or_independence_required": False,
        "fine_preparation_and_auxiliary_tapes_in_comparison": False,
        "terminal_flags_in_comparison": FLAGS,
    }


@lru_cache(maxsize=2)
def _source_channels(law):
    incident._law(law)
    return tuple(incident.seam_counts(e, law) for e in quantum_mixing.seams())


@lru_cache(maxsize=4)
def _herald_channels(law, moving):
    incident._law(law)
    if type(moving) is not bool:
        raise ValueError("Boolean movement flag required")
    states = native.states(binary=True)
    result = {}
    for seam in incident.incident_edges(CORE):
        rows, cols = [], []
        for coin, _, mapping in coins.coin_maps(law, seam):
            selected = [(i, o) for i, o in mapping
                        if (native.native(states[i], seam, coin)[1][1] != 0) == moving]
            for i, a in selected:
                for j, b in selected:
                    rows.append(66*a+b); cols.append(66*i+j)
        result[seam] = coo_matrix((np.ones(len(rows), np.int64), (rows, cols)),
                                  shape=(4356, 4356)).tocsr()
    return result


def exact_audit():
    """Integer source checks on every matrix unit, including quantum modes."""
    source = incident.exact_audit(CORE)
    e, scale = cut.completion_counts(CORE)
    weights, p, gamma = _schedule()
    block_of = np.empty(4356, dtype=int)
    for number, indices in enumerate(difference_blocks()):
        block_of[list(indices)] = number
    laws = {}
    for law in ("C", "E"):
        channels = _source_channels(law)
        for seam, a in zip(quantum_mixing.seams(), channels, strict=True):
            row, column = a.nonzero()
            if np.any(block_of[row] != block_of[column]):
                raise ArithmeticError("native channel changes a signed-difference block")
            if seam in weights and (not cut._zero(a@e-2*e) or not cut._zero(e@a-2*e)):
                raise ArithmeticError("a mixed local seam does not preserve and fix E_CORE")
            if seam in incident.incident_edges(CORE):
                if not cut._zero(_herald_channels(law, True)[seam]
                                 +_herald_channels(law, False)[seam]-a):
                    raise ArithmeticError("the two specified-tick herald branches are incomplete")
        total, denominator = common_caps.mixed_counts(law)
        if (not cut._zero(total-total.T) or not cut._zero(total@e-denominator*e)
                or not cut._zero(e@total-denominator*e)):
            raise ArithmeticError("mixed source contraction/completion premises failed")
        laws[law] = {"matrix_units": 4356, "each_mixed_seam_fixes_and_preserves_completion": True,
                     "both_actual_herald_channels_retained": True,
                     "signed_difference_blocks_preserved": True}
    if not cut._zero(e@e-scale*e):
        raise ArithmeticError("source completion is not a projection")
    # Prove the shared state preserves every core-containing bare algebra.
    sigma = reference()
    for mask in range(1024):
        cap = (*CORE, *(i+2 for i in range(10) if mask & (1 << i)))
        for unit in common_caps.regional.basis_supports(cap, completed=False):
            if len({sigma[i]/sigma[j] for i, j in unit}) != 1:
                raise ArithmeticError("common local state fails a core-containing cap algebra")
    return {"local_source": source, "laws": laws, "completion_tick_fraction": str(p),
            "mixed_attempt_gap_lower": str(gamma), "core_containing_cap_algebras": 1024,
            "native_auxiliary_coin_instrument": coins.effect_certificate()}


def _power_sum(matrix, exponent):
    """Return M^n and sum_{j<n} M^j, without inversion or long iteration."""
    _finite(matrix)
    power, total = np.eye(len(matrix)), np.zeros_like(matrix)
    base, block = matrix.copy(), np.eye(len(matrix))
    while exponent:
        with np.errstate(over="ignore", invalid="ignore"):
            if exponent & 1:
                total += power@block
                power = power@base
            exponent >>= 1
            if exponent:
                block += base@block
                base = base@base
        _finite(power, total, base, block)
    return power, total


def _finite(*arrays):
    if any(not np.isfinite(a).all() for a in arrays):
        raise ArithmeticError("nonfinite finite native channel execution")


@lru_cache(maxsize=1)
def difference_blocks():
    basis = native.states(binary=True)
    groups = {}
    for i, x in enumerate(basis):
        for j, y in enumerate(basis):
            z = tuple(a-b for a, b in zip(x, y, strict=True))
            groups.setdefault(z, []).append(66*i+j)
    return tuple(tuple(ix) for ix in groups.values())


def _block_maps(law, indices, *, warmup_ticks, completion_ticks, deadline, selector_rounds):
    """Actual finite native maps on one invariant signed-difference block."""
    incident._law(law)
    _integer(warmup_ticks, 0, 1000000)
    _integer(completion_ticks, 0, 1000000)
    _integer(deadline, 1, 4096)
    selected = selectors(selector_rounds)
    seams = quantum_mixing.seams()
    core_edges = incident.incident_edges(CORE)
    d = len(indices)
    warm_a, warm_r, move, stay = (np.zeros((d, d)) for _ in range(4))
    complete_a, complete_r = np.zeros((d, d)), np.zeros((d, d))
    warm_weight = float(F(selected["warmup"]["accepted_probability_per_address"]))
    for seam, source in zip(seams, _source_channels(law), strict=True):
        m = source[indices, :][:, indices].toarray()/2
        if seam in core_edges:
            warm_a += warm_weight*m
            move += warm_weight*_herald_channels(law, True)[seam][indices, :][:, indices].toarray()/2
            stay += warm_weight*_herald_channels(law, False)[seam][indices, :][:, indices].toarray()/2
        else:
            warm_r += warm_weight*m
        if seam in selected["completion_accepted_seam_weights"]:
            weighted = float(selected["completion_accepted_seam_weights"][seam])*m
            if set(seam).intersection(CAPS[0]):
                complete_a += weighted
            else:
                complete_r += weighted
    _finite(warm_a, warm_r, move, stay, complete_a, complete_r)
    warm_tail, warm_sum = _power_sum(warm_r, deadline)
    with np.errstate(over="ignore", invalid="ignore"):
        warm_tick = warm_a@warm_sum
        warm_failure = float(F(selected["warmup"]["timeout_probability"]))*warm_sum+warm_tail
    _finite(warm_tick, warm_failure)
    warm_power, warm_accumulated = _power_sum(warm_tick, warmup_ticks)
    with np.errstate(over="ignore", invalid="ignore"):
        heralded = move@warm_sum@warm_power
    _finite(heralded)
    complete_tail, complete_sum = _power_sum(complete_r, deadline)
    with np.errstate(over="ignore", invalid="ignore"):
        complete_tick = complete_a@complete_sum
        complete_failure = float(selected["completion_timeout_probability"])*complete_sum+complete_tail
    _finite(complete_tick, complete_failure)
    complete_power, complete_accumulated = _power_sum(complete_tick, completion_ticks)
    with np.errstate(over="ignore", invalid="ignore"):
        result = {"success": complete_power@heralded,
                  "herald_failure": stay@warm_sum@warm_power,
                  "warmup_timeout": warm_failure@(warm_accumulated+warm_power),
                  "completion_timeout": complete_failure@complete_accumulated@heralded}
    _finite(*result.values())
    return result


def population_maps(precision_bits=24, deadline=256, selector_rounds=32, *, law="C"):
    """Fresh66-by66 matrices, columns=input, with all four output flags.

    These are numerical executions of the finite source protocol; the exact
    rational budget is the certificate. No branch is normalized or discarded.
    """
    b = budget(precision_bits, deadline, selector_rounds)
    return _block_maps(law, tuple(67*i for i in range(66)), warmup_ticks=b["warmup_ticks"],
                       completion_ticks=b["completion_ticks"], deadline=deadline,
                       selector_rounds=selector_rounds)


def execute(density, *, law="C", warmup_ticks=2, completion_ticks=2, deadline=3, selector_rounds=2):
    """Bounded actual quantum execution, with up to three reference levels.

    Explicit small counts support adversarial complete-instrument tests.
    Call budget() separately for certified preparation counts. Subnormalized
    positive inputs, including zero, are preserved without renormalization.
    """
    incident._law(law)
    if (type(density) is not np.ndarray or density.dtype.kind not in "fc" or density.ndim != 2
            or density.shape[0] != density.shape[1] or density.shape[0] % 66
            or not 66 <= density.shape[0] <= 198 or not np.isfinite(density).all()):
        raise ValueError("finite plain66-sector density with up to three reference levels required")
    for value, low, high in ((warmup_ticks, 0, 1000000), (completion_ticks, 0, 1000000), (deadline, 1, 4096)):
        _integer(value, low, high)
    selectors(selector_rounds)
    with np.errstate(over="ignore", invalid="ignore"):
        density = density.astype(complex, copy=False)
        trace = np.trace(density)
        # sqrt(d)*Frobenius dominates the nuclear norm. Screen aggregate
        # anti-Hermitian error, rather than accepting d^2 small defects.
        hermitian_error = math.sqrt(len(density))*np.linalg.norm(density-density.conj().T, "fro")/2
        validation_state = density/2+density.conj().T/2
    if (not np.isfinite(density).all() or not np.isfinite(trace)
            or not np.isfinite(hermitian_error) or not np.isfinite(validation_state).all()
            or hermitian_error > 1e-13):
        raise ValueError("positive subnormalized Hermitian density required")
    with np.errstate(over="ignore", invalid="ignore"):
        eigenvalues = np.linalg.eigvalsh(validation_state)
        negative_mass = -eigenvalues[eigenvalues < 0].sum()
        validation_error = hermitian_error+negative_mass+max(0., trace.real-1., -trace.real)
    if (not np.isfinite(eigenvalues).all() or not np.isfinite(validation_error)
            or validation_error > 1e-13):
        raise ValueError("positive subnormalized Hermitian density required")
    # validation_state is used only for this screen. The original, neither
    # symmetrized nor normalized input is propagated through every branch.
    ref = len(density)//66
    tensor = density.reshape(66, ref, 66, ref).transpose(0, 2, 1, 3).reshape(4356, ref*ref)
    output = {flag: np.zeros_like(tensor, dtype=complex) for flag in FLAGS}
    for ix in difference_blocks():
        if not np.any(tensor[list(ix)]):
            continue
        maps = _block_maps(law, ix, warmup_ticks=warmup_ticks, completion_ticks=completion_ticks,
                           deadline=deadline, selector_rounds=selector_rounds)
        for flag in FLAGS:
            _finite(maps[flag])
            with np.errstate(over="ignore", invalid="ignore"):
                branch = maps[flag]@tensor[list(ix)]
            _finite(branch)
            output[flag][list(ix)] = branch
    return {flag: a.reshape(66, 66, ref, ref).transpose(0, 2, 1, 3).reshape(density.shape)
            for flag, a in output.items()}


def _json(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {str(k): _json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json(v) for v in value]
    return value


def run():
    b = budget()
    maps = population_maps()
    if (type(maps) is not dict or set(maps) != set(FLAGS)
            or any(type(m) is not np.ndarray or m.shape != (66, 66) or m.dtype.kind != "f"
                   or not np.isfinite(m).all() or np.any(m < 0) or np.any(m > 1+5e-10)
                   for m in maps.values())):
        raise ArithmeticError("each retained population channel must be finite and nonnegative")
    with np.errstate(over="ignore", invalid="ignore"):
        masses = {flag: m.sum(axis=0) for flag, m in maps.items()}
        total = sum(maps.values())
    _finite(total, *masses.values())
    if any(np.any(m > 1+5e-10) for m in masses.values()):
        raise ArithmeticError("a retained population channel exceeds unit mass")
    q = float(b["uniform_ideal_herald_probability"])
    raw = np.array(list(map(float, incident.reference(CORE)[0])))
    sigma = np.array(list(map(float, reference())))
    ideal = {"success": q*sigma[:, None], "herald_failure": (1/66-q*raw)[:, None],
             "warmup_timeout": np.zeros((66, 1)), "completion_timeout": np.zeros((66, 1))}
    with np.errstate(over="ignore", invalid="ignore"):
        deviations = sum(np.abs(maps[f]-ideal[f]).sum(axis=0) for f in FLAGS)/2
        mass_residual = np.max(np.abs(total.sum(axis=0)-1))
    _finite(deviations, mass_residual)
    if (mass_residual > 5e-10
            or np.max(deviations) > float(b["complete_flagged_reference_halftrace_upper"])+5e-10):
        raise ArithmeticError("finite source population execution exceeds its certificate")
    return {"schema": "oph.native_modular_reference.local_preparation.v1", "all_pass": True,
            "construction": "finite native-selector incident-clock common preparation",
            "reference_weights_and_multiplicities": {str(k): v for k, v in Counter(reference()).items()},
            "budget": _json(b), "exact_source_audit": exact_audit(),
            "population_execution": {"all_input_columns": 66, "retained_terminal_flags": list(FLAGS),
                "maximum_total_probability_residual_diagnostic": float(mass_residual),
                "maximum_complete_halftrace_diagnostic": float(np.max(deviations)),
                "basis_zero_output_flag_probabilities_diagnostic":
                    {flag: float(m[:, 0].sum()) for flag, m in maps.items()}},
            "source_sha256": {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
                              for path in SOURCE_PATHS}}


def compare_receipt(supplied, actual):
    """Keep replay tolerance, but never extend probability/error domains."""
    common_caps.full.corner.compare_receipt(supplied, actual)
    for receipt in (supplied, actual):
        row = receipt["population_execution"]
        residual = row["maximum_total_probability_residual_diagnostic"]
        distance = row["maximum_complete_halftrace_diagnostic"]
        probabilities = row["basis_zero_output_flag_probabilities_diagnostic"]
        if set(probabilities) != set(FLAGS):
            raise ValueError("all four native preparation flags required")
        values = [residual, distance, *probabilities.values()]
        if any(type(v) not in (int, float) or not math.isfinite(v) or v < 0 for v in values):
            raise ValueError("nonnegative finite preparation diagnostics and probabilities required")
        bound = float(F(receipt["budget"]["complete_flagged_reference_halftrace_upper"]))
        if (residual > 5e-10 or distance > bound+5e-10
                or any(v > 1+5e-10 for v in probabilities.values())):
            raise ValueError("preparation diagnostic or probability exceeds its domain")
        mass = sum(probabilities.values())
        if not math.isfinite(mass) or abs(mass-1) > 5e-10 or abs(mass-1) > residual+1e-13:
            raise ValueError("complete probability mass disagrees with its maximum residual")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write", type=Path)
    group.add_argument("--check", type=Path)
    args = parser.parse_args()
    report = run()
    if args.write is not None:
        with args.write.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(json.dumps(report, indent=2, sort_keys=True)+"\n")
    if args.check:
        compare_receipt(strict_json.load(args.check), report)
    if not args.write and not args.check:
        print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
