"""Common native modular clocks for overlapping caps and their exact boundary.

For a faithful diagonal state on the two-occupation shell over N>=3 ports,
all bare two-port number-preserving algebras are modular invariant iff
w_ij=lambda*r_i*r_j. Equivalently every four distinct ports have equal
matching products w_ij*w_kl=w_ik*w_jl=w_il*w_jk. This is an exact finite
characterization, not a restriction on arbitrary nondiagonal references.
For N=3 the condition is vacuous and every positive weight triple factors.

Necessity follows because |i,k><j,k| must have the same modular ratio for
every spectator k. Sufficiency follows by cancellation of the spectator
factor. Set lambda=w01*w02/w12, r0=1, r1=w12/w02 and ri=wi1/w01 for i>=2;
the four-point identities prove every remaining pair. The generator is
constant minus sum_i log(r_i) N_i, and all regional algebras are invariant.
Native connected symmetric repair preserves such a state only when all r_i
are equal, so its stationary modular action is trivial. Nontrivial fields
are not selected by this factorization test.

There is a constructive common-state domain without demanding all pairs.
Choose caps R=(0,1,5,7), S=(0,1,2,8). Select either cap fairly, then a uniform
complement seam and fair native placement coin; use the unread marginal.
Their union has 20 edges connected on the ten ports outside I=R intersect S.
The native C/E fixed algebra is A_I. The weighted gap is >=1/600 and the
same source preparation theorem gives a reference-uniform convergence bound
to E_I. Acting on the native Palm state selects one common E_I(Palm), with
weights 68/4500,69/4500,60/4500 according to occupation of I.

This density belongs to A_I, hence its one modular action preserves every
A_R with I subset R. It is central in A_I: the overlap's own clock is
trivial. Nontrivial cap probes transfer occupation between the protected
core and another port of their cap. The 2042-pulse common source word and
two complete bit-read Ramsey experiments below implement this same action.
Protected ports, cap/seam scheduling, inaccessible preparation records and
coherent preparation/read services remain specified operational choices.
No arbitrary local Hamiltonian, physical time or spacetime net is supplied.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import combinations
import json
import math
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix

from tools import strict_json
from ..native_equilibration import construction as mixing, cut_completion as cut
from . import regional, completed_clock as previous, full_flow as full, preparation

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name("common_caps_results.json")
CAPS = ((0, 1, 5, 7), (0, 1, 2, 8))
CORE = (0, 1)
WORD_ERROR_GATE = F(1, 10**9)
SOURCE_PATHS = tuple(sorted(set((
    *full.corner.SOURCE_PATHS, *full.corner.report.SOURCE_PATHS,
    *cut.SOURCE_PATHS, *mixing.SOURCE_PATHS, *preparation.burn_in.SOURCE_PATHS,
    "model/native_modular_reference/common_caps.py", "model/native_modular_reference/regional.py",
    "model/native_modular_reference/completed_clock.py", "model/native_modular_reference/preparation.py",
    "model/native_modular_reference/full_flow.py", "model/common_source/exact_output.py",
    "model/common_source/output_certificate.py",
))))


def factor_characterization(vertices, weights):
    if type(vertices) is not int or vertices < 3:
        raise ValueError("at least three integer ports required")
    expected = set(combinations(range(vertices), 2))
    if (type(weights) is not dict
            or any(type(pair) is not tuple or len(pair) != 2
                   or any(type(port) is not int for port in pair) for pair in weights)
            or set(weights) != expected
            or any(type(w) not in (int, F) or w <= 0 for w in weights.values())):
        raise ValueError("one exact positive weight per unordered distinct pair required")
    w = lambda i, j: F(weights[tuple(sorted((i, j)))])
    witness = None
    for i, j, k, l in combinations(range(vertices), 4):
        products = (w(i, j)*w(k, l), w(i, k)*w(j, l), w(i, l)*w(j, k))
        if len(set(products)) != 1:
            witness = {"ports": [i, j, k, l], "matching_products": list(map(str, products))}
            break
    scale = w(0, 1)*w(0, 2)/w(1, 2)
    factors = (F(1), w(1, 2)/w(0, 2), *(w(i, 1)/w(0, 1) for i in range(2, vertices)))
    compatible = all(w(i, j) == scale*factors[i]*factors[j] for i, j in expected)
    if compatible != (witness is None):
        raise ArithmeticError("four-point invariance and constructive factorization disagree")
    return {"all_bare_two_port_algebras_invariant": compatible,
            "factor_scale": str(scale) if compatible else None,
            "relative_factors": list(map(str, factors)) if compatible else None,
            "noninvariance_four_port_witness": witness}


def native_weights(reference):
    states = regional.source()[0].states
    if len(reference) != 66:
        raise ValueError("one weight per native state required")
    return {tuple(p for p, bit in enumerate(x) if bit): value
            for x, value in zip(states, reference, strict=True)}


def cap_schedule(caps=CAPS):
    if type(caps) is not tuple or not caps:
        raise ValueError("a nonempty tuple of declared caps required")
    caps = tuple(regional.region(r) for r in caps)
    if len(set(caps)) != len(caps):
        raise ValueError("distinct caps required")
    core = tuple(sorted(set.intersection(*(set(r) for r in caps))))
    active = tuple(p for p in range(12) if p not in core)
    if len(active) < 3:
        raise ValueError("at least three active ports exclude empty/full coherences")
    weights = {}
    for cap in caps:
        cut.cut_certificate(cap)  # each constituent is an admitted cut protocol
        seams = [e for e in regional.source()[0].seams if not set(e).intersection(cap)]
        for seam in seams:
            weights[seam] = weights.get(seam, F())+F(1, len(caps)*len(seams))
    graph = mixing.graph_certificate(len(active), tuple((active.index(a), active.index(b))
                                                        for a, b in sorted(weights)))
    alpha = len(weights)*min(weights.values())
    gap = alpha*min(F(graph["quantum_gap_lower_per_attempt"]), F(1, 2*len(weights)))
    return core, weights, {"caps": [list(r) for r in caps], "protected_intersection": list(core),
                          "native_seam_union": len(weights), "uniform_component_weight": str(alpha),
                          "graph": graph, "gap_lower_per_attempt": str(gap)}


def mixed_counts(law, caps=CAPS):
    cap_schedule(caps)
    rows = [cut.native_cut_counts(cap, law) for cap in caps]
    scale = math.lcm(*(len(caps)*denominator for _, denominator in rows))
    total = sum((matrix*(scale//(len(caps)*denominator)) for matrix, denominator in rows),
                coo_matrix((4356, 4356), dtype=np.int64).tocsr())
    return total, scale


@lru_cache(maxsize=1)
def reference():
    core, _, _ = cap_schedule()
    if core != CORE:
        raise ArithmeticError("protected common cap intersection changed")
    sigma = cut.completed_palm(CORE)
    if set(sigma) != {F(68, 4500), F(69, 4500), F(60, 4500)}:
        raise ArithmeticError("common native Palm completion weights changed")
    states = regional.source()[0].states
    if any(w != (F(68, 4500), F(69, 4500), F(60, 4500))[sum(x[p] for p in CORE)]
           for w, x in zip(sigma, states, strict=True)):
        raise ArithmeticError("common density is not central in the protected-core algebra")
    return sigma


def exact_audit():
    sigma = reference()
    e, es = cut.completion_counts(CORE)
    density = np.zeros(4356, np.int64)
    for i, w in enumerate(sigma):
        density[67*i] = int(4500*w)
    laws = {}
    for law in ("C", "E"):
        p, ps = mixed_counts(law)
        if (not cut._zero(p@e-ps*e) or not cut._zero(e@p-ps*e)
                or not cut._zero(p-p.T) or not np.array_equal(p@density, ps*density)):
            raise ArithmeticError("actual mixed native channel lost the common completion")
        laws[law] = {"matrix_units_checked": 4356, "native_integer_denominator": ps,
                     "shared_completion_fixed_and_preserved": True}
    # Exact invariance on every cap containing the core, without assuming a
    # common state independently for each observer.
    for mask in range(1 << 10):
        ports = (*CORE, *(p for bit, p in enumerate(range(2, 12)) if mask & (1 << bit)))
        for unit in regional.basis_supports(ports, completed=False):
            if len({sigma[i]/sigma[j] for i, j in unit}) != 1:
                raise ArithmeticError("one common modular action failed a core-containing cap")
    for unit in regional.basis_supports(CORE, completed=False):
        if any(sigma[i] != sigma[j] for i, j in unit):
            raise ArithmeticError("protected intersection unexpectedly acquired a moving clock")
    return {"native_laws": laws, "common_reference_cap_algebras_checked": 1024,
            "protected_intersection_modular_action": "identity",
            "nontrivial_cap_probes_require_transfer_between_core_and_remaining_cap": True}


def unprotected_control():
    """Actual interleaved caps whose seam union has no protected ports."""
    caps = (CAPS[0], (2, 3, 4, 6))
    core, _, schedule = cap_schedule(caps)
    if core:
        raise ArithmeticError("unprotected connected control gained a conserved cap")
    rows, cols = np.meshgrid(np.arange(66)*67, np.arange(66)*67, indexing="ij")
    projection = coo_matrix((np.ones(66*66, dtype=np.int64), (rows.ravel(), cols.ravel())),
                            shape=(4356, 4356)).tocsr()
    laws = {}
    for law in ("C", "E"):
        p, scale = mixed_counts(law, caps)
        if not cut._zero(p@projection-scale*projection) or not cut._zero(projection@p-scale*projection):
            raise ArithmeticError("connected unprotected source lost the tracial projection")
        laws[law] = {"matrix_units_checked": 4356, "fixed_algebra_dimension": 1}
    return {"schedule": schedule, "native_laws": laws,
            "unique_fixed_binary_sector_state": "I66/66",
            "stationary_modular_action": "identity",
            "scope": "unread connected native equilibrium; does not exclude heralded or protected references"}


def preparation_budget(attempts=369, precision_bits=24, window=F(16)):
    palm = preparation.budget(attempts, window)
    _, _, schedule = cap_schedule()
    joined = mixing.preparation_budget(66, F(schedule["gap_lower_per_attempt"]), precision_bits)
    eta = F(joined["trace_distance_upper"])
    t, u = 66*palm["warmup_half_trace_error_upper"], eta/min(reference())
    radius = t+u+t*u
    if radius >= 1:
        raise ValueError("the common prepared reference is not certified faithful")
    return {"warmup_native_attempts": attempts, "global_movement_herald_attempts": 1,
            "mixed_cap_native_attempts_after_success": joined["attempts"],
            "mixed_cap_trace_error_upper": eta,
            "operator_sandwich_radius": radius,
            "actual_common_reference_modular_half_diamond_upper": min(F(1), window*radius/(1-radius)),
            "complete_flagged_reference_half_trace_upper":
                palm["warmup_half_trace_error_upper"]+preparation.Q*eta,
            "minimum_work_eigenvalue_lower":
                (1-radius)*preparation.Q*min(reference())/palm["herald_success_interval"][1],
            "expected_native_preparation_attempts_upper":
                palm["expected_attempts_with_full_warmup_repeated_each_trial_upper"]+joined["attempts"],
            "individual_cap_marginals_preserved_by_other_cap_protocol": False}


@lru_cache(maxsize=1)
def target_ratios():
    values = [F(1)]*512
    for code, w in zip(full.support()[0], reference(), strict=True):
        values[code] = 60*w
    return tuple(values)


def _ratios(values):
    if (type(values) is not tuple or len(values) != 512
            or any(type(value) not in (int, F) or value <= 0 for value in values)):
        raise ValueError("512 original exact positive diagonal ratios in a tuple required")
    return tuple(map(F, values))


def logarithm_interval(ratio):
    """Exact power-of-two range reduction, preserving the original local box."""
    if type(ratio) is not F or ratio <= 0:
        raise ValueError("an original exact positive Fraction ratio required")
    return _logarithm_interval(ratio)


@lru_cache(maxsize=128)
def _logarithm_interval(ratio):
    if F(1, 2) <= ratio <= 2:
        return previous.logarithm_interval(ratio)
    exponent = ratio.numerator.bit_length()-ratio.denominator.bit_length()
    reduced = ratio/F(2)**exponent
    if reduced < 1:
        reduced *= 2
        exponent -= 1
    lo, hi = previous.logarithm_interval(reduced)
    log2_lo, log2_hi = previous.logarithm_interval(F(2))
    shift_lo, shift_hi = sorted((exponent*log2_lo, exponent*log2_hi))
    return lo+shift_lo, hi+shift_hi


def phase_enclosure(s, ratio):
    """Enclose exp(-i*s*log(ratio)); return exact real, imag and error radius."""
    s = full.parameter(s)
    lo, hi = logarithm_interval(ratio)
    if F(1, 2) <= ratio <= 2:
        return previous.phase_enclosure(s, ratio)
    theta, halvings = -s*(lo+hi)/2, 0
    while abs(theta) > 8:
        theta /= 2
        halvings += 1
    real, imag, pr, pi = F(1), F(), F(1), F()
    for k in range(1, 65):
        pr, pi = -theta*pi/k, theta*pr/k
        real += pr
        imag += pi
    radius = abs(theta)**65/math.factorial(65)+abs(s)*(hi-lo)/(2*2**halvings)
    # Dyadic centers keep repeated squaring finite; both component errors
    # enter the Euclidean error bound by their exact L1 majorant.
    scale = 1 << 128
    for step in range(halvings+1):
        re = F((real*scale).__floor__(), scale)
        im = F((imag*scale).__floor__(), scale)
        radius += abs(real-re)+abs(imag-im)
        real, imag = re, im
        if step < halvings:
            real, imag = real*real-imag*imag, 2*real*imag
            radius *= 2+radius
    return real, imag, radius


@lru_cache(maxsize=8)
def _walsh_counts(values):
    groups = {v: tuple(code for code, w in enumerate(values) if w == v)
              for v in sorted(set(values)) if v != 1}
    return tuple((v, tuple(sum((-1)**((code & mask).bit_count())
                              for code in indices) for mask in range(512)))
                 for v, indices in groups.items())


def walsh_counts():
    return _walsh_counts(target_ratios())


def diagonal_word(s, ratios):
    """Compile the declared positive diagonal reference on every code address."""
    s = full.parameter(s)
    ratios = _ratios(ratios)
    angles = [0.]*512
    for ratio, counts in _walsh_counts(ratios):
        lo, hi = logarithm_interval(ratio)
        angle = float(-s*(lo+hi)/2)
        angles = [a+angle*n/512 for a, n in zip(angles, counts, strict=True)]
    tape = [full.compiler.pulse("I", angles[0], [0])]
    for pivot in range(9):
        size = 1 << pivot
        for k in range(size):
            gray = k ^ (k >> 1)
            tape.append(full.compiler.pulse("Z", angles[(1 << pivot)|gray], [pivot]))
            next_k = (k+1) % size
            changed = gray ^ (next_k ^ (next_k >> 1))
            if changed:
                tape.extend(full.compiler.cnot(changed.bit_length()-1, pivot))
    return tape


def word(s):
    return diagonal_word(s, target_ratios())


@lru_cache(maxsize=4)
def _executed(s):
    a = full.compiler.execute(word(s), 9)
    return a.dtype.str, a.tobytes()


def matrix(s):
    dtype, data = _executed(full.parameter(s))
    return np.frombuffer(data, dtype=np.dtype(dtype)).reshape(512, 512)


def certify_diagonal(s, ratios, array):
    """Bound the original array's operator error, including every padding entry."""
    s = full.parameter(s)
    ratios = _ratios(ratios)
    if type(array) is not np.ndarray:
        raise ValueError("an original base array required")
    array = full.exact._array(array, (512, 512))
    scale, phases = 1 << 128, {}
    for ratio in set(ratios):
        re, im, error = phase_enclosure(s, ratio)
        phases[ratio] = (F((re*scale).__floor__(), scale),
                         F((im*scale).__floor__(), scale), error+F(2, scale))
    squared = F()
    for i, row in enumerate(array):
        for j, value in enumerate(row):
            ar, ai = full.exact._component(value)
            re, im = phases[ratios[i]][:2] if i == j else (F(), F())
            squared += (ar-re)**2+(ai-im)**2
    return full.exact._sqrt_interval(squared)[1]+full.exact._sqrt_interval(F(512))[1]*max(x[2] for x in phases.values())


def certify(s, array=None):
    s = full.parameter(s)
    return certify_diagonal(s, target_ratios(), matrix(s) if array is None else array)


@lru_cache(maxsize=4)
def _check(s):
    full.corner.report._gate(certify(s), WORD_ERROR_GATE, "common-cap native modular source")


def checked_bound(s):
    _check(full.parameter(s))
    return {"operator_error_upper": WORD_ERROR_GATE,
            "half_diamond_upper": WORD_ERROR_GATE*(2+WORD_ERROR_GATE)/2,
            "channel_norm_upper": (1+WORD_ERROR_GATE)**2}


def mixer_word(sign, bits, *, compact=False):
    if (type(sign) is not int or sign not in (-1, 1)
            or type(bits) is not tuple or any(type(bit) is not int for bit in bits)
            or bits not in ((0, 8), (0, 3, 6)) or type(compact) is not bool):
        raise ValueError("declared signed native pair mixer required")
    positions = tuple(range(len(bits))) if compact else bits
    fan = [event for target in positions[1:] for event in full.compiler.cnot(positions[0], target)]
    reverse = [event for target in reversed(positions[1:]) for event in full.compiler.cnot(positions[0], target)]
    return [*fan, full.compiler.pulse("Y", sign*math.pi/4, [positions[0]]), *reverse]


def mixer_bound(sign, bits):
    a = full.compiler.execute(mixer_word(sign, bits, compact=True), len(bits))
    d = 1 << len(bits)
    lo, hi = full.exact._sqrt_interval(F(1, 2))
    middle, radius, squared = (lo+hi)/2, (hi-lo)/2, F()
    for i in range(d):
        for j in range(d):
            wanted = middle if i == j else sign*(1 if j & 1 else -1)*middle if i == (j ^ (d-1)) else F()
            re, im = full.exact._component(a[i, j]); squared += (re-wanted)**2+im*im
    eps = full.exact._sqrt_interval(squared)[1]+full.exact._sqrt_interval(F(2*d))[1]*radius
    gate = full.MIXER_ERROR_GATE
    full.corner.report._gate(eps, gate, "common-cap Ramsey mixer")
    return {"half_diamond_upper": gate*(2+gate)/2, "channel_norm_upper": (1+gate)**2}


def witness(cap_index, s=F(16)):
    if type(cap_index) is not int or cap_index not in (0, 1):
        raise ValueError("one of the two declared cap probes required")
    s = full.parameter(s)
    start, end, bits = ((328, 73, (0, 8)), (0, 73, (0, 3, 6)))[cap_index]
    codes = full.support()[0]; i, j = codes.index(start), codes.index(end)
    states = regional.source()[0].states
    if (reference()[j]/reference()[i] != F(69, 68)
            or any(states[i][p] != states[j][p] for p in range(12) if p not in CAPS[cap_index])):
        raise ArithmeticError("common-gap probe changes its cap-exterior context")
    initial = np.zeros(512, complex); initial[start] = 1
    prepared = full.compiler.execute(mixer_word(-1, bits), 9, initial)
    final = full.compiler.execute(mixer_word(1, bits), 9, matrix(s)@prepared)
    probabilities = [float(np.vdot(final[b::2], final[b::2]).real) for b in (0, 1)]
    re, _, error = previous.phase_enclosure(s, F(69, 68))
    match = (1+re)/2
    mixer_bound(-1, bits)
    bound = full.corner.composite_bound(checked_bound(s), mixer_bound(1, bits))
    discrepancy = (abs(probabilities[0]-float(match))+abs(probabilities[1]-float(1-match)))/2
    if discrepancy > float(bound+error/2):
        raise ArithmeticError("complete common-clock Ramsey record exceeds source certificate")
    return {"cap": list(CAPS[cap_index]), "native_codes": [start, end],
            "complete_bit_probabilities": probabilities,
            "ideal_match_law": "(1+cos(s log(69/68)))/2",
            "complete_record_half_l1_diagnostic": discrepancy,
            "whole_sequence_half_diamond_upper": full.compiled._outward(bound),
            "source_pulses": len(word(s))+2*len(mixer_word(1, bits)),
            "coherent_selected_pair_and_complete_bit0_PVM_remain_admitted": True}


def run():
    sigma = reference()
    _, _, schedule = cap_schedule()
    common = preparation_budget()
    gate = checked_bound(F(16))
    common["source_to_actual_common_modular_half_diamond_upper"] = gate["half_diamond_upper"]+common["actual_common_reference_modular_half_diamond_upper"]
    common["source_then_any_core_cap_restrict_vs_actual_cap_modular_half_diamond_upper"] = common["source_to_actual_common_modular_half_diamond_upper"]
    common["complete_flagged_reference_and_source_half_trace_upper"] = (
        common["complete_flagged_reference_half_trace_upper"]
        +preparation.budget(369, F(16))["herald_success_interval"][1]*gate["half_diamond_upper"])
    result = {"schema": "oph.native_modular_reference.common_caps.v1", "all_pass": True,
              "source_sha256": {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
              "schedule": schedule, "exact_source_audit": exact_audit(),
              "unprotected_connected_control": unprotected_control(),
              "reference_weight_census": {str(w): sigma.count(w) for w in sorted(set(sigma))},
              "all_pair_characterization": {"Palm": factor_characterization(12, native_weights(regional.source()[2])),
                  "one_cap_completion": factor_characterization(12, native_weights(cut.completed_palm(CAPS[0]))),
                  "common_cap_completion": factor_characterization(12, native_weights(sigma)),
                  "tracial": factor_characterization(12, native_weights((F(1, 66),)*66))},
              "finite_common_reference_budget": preparation._compact(common),
              "retained_modular_parameter": "16",
              "common_source_pulses": len(word(F(16))), "source_operator_error_gate": str(WORD_ERROR_GATE),
              "cap_readouts": [witness(i) for i in (0, 1)],
              "scope": {"one_common_native_prepared_state_for_1024_core_containing_caps": True,
                        "all_port_regions_or_spacetime_net_supplied": False,
                        "protected_intersection_has_nontrivial_modular_clock": False,
                        "individual_cap_marginals_preserved_by_mixed_preparation": False,
                        "future_access_or_feedback_from_preparation_tape_covered": False,
                        "other_modular_parameters_require_fresh_word_certification": True,
                        "preparation_tape_may_be_stored_if_future_protocol_ignores_it": True,
                        "arbitrary_nondiagonal_reference_classified_by_factor_test": False,
                        "physical_time_temperature_or_passive_modular_evolution_identified": False}}
    return json.loads(preparation.canonical(result))


def _receipt_domains(receipt):
    """Screen raw numerical diagnostics without changing replay tolerances."""
    if (type(receipt) is not dict or receipt.get("schema") != "oph.native_modular_reference.common_caps.v1"
            or type(receipt.get("cap_readouts")) is not list or len(receipt["cap_readouts"]) != 2):
        raise ValueError("a complete common-cap receipt required")
    for row in receipt["cap_readouts"]:
        if type(row) is not dict:
            raise ValueError("a complete common-cap readout required")
        probabilities = row.get("complete_bit_probabilities")
        error = row.get("complete_record_half_l1_diagnostic")
        raw_bound = row.get("whole_sequence_half_diamond_upper")
        try:
            bound = F(raw_bound) if type(raw_bound) is str else None
        except (ValueError, ZeroDivisionError):
            bound = None
        if bound is None or not 0 <= bound < 1 or str(bound) != raw_bound:
            raise ValueError("a canonical nonnegative complete-word bound required")
        if (type(probabilities) is not list or len(probabilities) != 2
                or any(type(p) is not float or not math.isfinite(p)
                       or not 0 <= p <= 1+2*float(bound) for p in probabilities)
                or abs(math.fsum(probabilities)-1) > 2*float(bound)):
            raise ValueError("complete nonnegative bit probabilities within the word mass bound required")
        if (type(error) is not float or not math.isfinite(error)
                or not 0 <= error <= float(bound)):
            raise ValueError("a nonnegative complete-record error within the word bound required")


def compare_receipt(expected, actual):
    full.corner.compare_receipt(expected, actual)
    _receipt_domains(expected)
    _receipt_domains(actual)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write", type=Path); group.add_argument("--check", type=Path)
    args = parser.parse_args(); result = run()
    if args.check:
        compare_receipt(strict_json.load(args.check), result)
        print("Common native cap clocks and exact compatibility boundary verified.")
    elif args.write:
        with args.write.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(full.corner.canonical(result)+"\n")
    else:
        print(full.corner.canonical(result))


if __name__ == "__main__":
    main()
