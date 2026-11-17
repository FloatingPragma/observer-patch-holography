"""One native entropy completion and its source-compiled regional clock.

The actual complement-only C/E protocol converges to E_R. Starting with the
native Palm reference gives bar_sigma=E_R(sigma), with no new Gibbs state
chosen. Its tensor form makes the bare regional algebra modular invariant.
For the four-port cap its charge-one gap is log(38/37), and its charge-two
gap is log(5/4). Both follow from actual native movement degrees.

Heisenberg convention is rho^(is) A rho^(-is). The dual source state word
implements bar_sigma^(-is) up to common native exp(is log60). A second word
implements V_s=bar_sigma^(-is) sigma^(is). Since both references are diagonal,
V_s U_Palm(s)=U_completed(s), V_(s+t)=V_s V_t, and V is the state-action
version of the finite relative modular cocycle. The Heisenberg cocycle is
V_s*. This is a change-of-reference identity, not equality of the two clocks.

Each word reuses the existing 2042-pulse closed Gray Walsh schedule and
existing source I,Z,H,J controls on all 512 addresses. Original binary arrays
are certified against rational logarithm/phase enclosures; every nonnative
address is an identity extension. No numerical array is renormalized.

Let delta be the native Palm warm-up half-trace bound, t=66delta, and eta
the complement-repair bound. Since bar_sigma>=I/75 and the cut map is CP,
r=t+75eta+t*75eta gives an operator sandwich around bar_sigma for the actual
unnormalized successful state, up to the scalar herald normalization.
Its modular half-diamond error is <=|s|r/(1-r). Every actual cap marginal
is preserved throughout cut repair, so the sharper cap-only bound uses t.
This joins the SAME prepared state, entropy completion and clock, while
retaining the cut schedule, discarded records and active coherent services.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import math
from pathlib import Path

import numpy as np

from tools import strict_json
from ..native_equilibration import cut_completion as cut
from . import full_flow as full, preparation, regional

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name("completed_clock_results.json")
WORD_ERROR_GATE = F(1, 10**9)
CAP = regional.CAP


@lru_cache(maxsize=1)
def reference():
    sigma = cut.completed_palm(CAP)
    native_sigma = regional.source()[2]
    independent = [F()] * 66
    bs = regional.blocks(CAP)
    for q, local_weights in regional.bare_reference(CAP):
        group = [b for b in bs if b.charge == q]
        multiplicity = sum(len(b.outside) for b in group)
        for a, w in enumerate(local_weights):
            for block in group:
                for i in block.indices[a]:
                    independent[i] = w/multiplicity
    if tuple(independent) != sigma or min(sigma) != F(1, 75):
        raise ArithmeticError("entropy completion and modular reference do not share a native state")
    if sum(abs(a-b) for a, b in zip(sigma, native_sigma, strict=True))/2 != F(8, 175):
        raise ArithmeticError("Palm completion disturbance changed")
    for unit in regional.basis_supports(CAP, completed=False):
        if len({sigma[i]/sigma[j] for i, j in unit}) != 1:
            raise ArithmeticError("completed reference does not preserve the bare cap algebra")
    return sigma


def ratios(kind):
    if kind not in ("clock", "cocycle"):
        raise ValueError("clock or cocycle word required")
    sigma = regional.source()[2]
    completed = reference()
    result = [F(1)] * 512
    for code, a, b in zip(full.support()[0], completed, sigma, strict=True):
        result[code] = 60*a if kind == "clock" else a/b
    return tuple(result)


def branch_covariance_audit():
    """Every actual labelled complement Kraus branch commutes with bar_sigma."""
    states = regional.source()[0].states
    index = {x: i for i, x in enumerate(states)}
    seams = tuple(e for e in regional.source()[0].seams if not set(e).intersection(CAP))
    sigma = reference()
    counts = {}
    for law in ("C", "E"):
        branches = set()
        checked = 0
        for seam in seams:
            for coin in (0, 1):
                for i, state in enumerate(states):
                    output, (_, current) = cut.native.native(state, seam, coin)
                    j = index[output]
                    label = current if law == "C" else state[seam[0]]-state[seam[1]]
                    branches.add((seam, coin, label))
                    if (sigma[i] != sigma[j]
                            or any(state[p] != output[p] for p in CAP)):
                        raise ArithmeticError("native labelled cut branch violates completed-clock covariance")
                    checked += 1
        counts[law] = {"labelled_Kraus_branches": len(branches),
                       "nonzero_source_Kraus_entries_checked": checked}
    return {"laws": counts, "complete_current_records_and_arbitrary_references_covered": True,
            "ideal_completed_modular_action_commutes_with_each_Kraus_operator": True,
            "scope": "complete modular words between complement attempts; no midword measurement"}
@lru_cache(maxsize=32, typed=True)
def logarithm_interval(ratio):
    if type(ratio) is not F or not F(1, 2) <= ratio <= 2:
        raise ValueError("an exact positive ratio in [1/2,2] required")
    x = (ratio-1)/(ratio+1)
    value = 2*sum((x**(2*j+1)/(2*j+1) for j in range(32)), F())
    remainder = 2*abs(x)**65/(65*(1-x*x))
    return value-remainder, value+remainder


@lru_cache(maxsize=128, typed=True)
def phase_enclosure(s, ratio):
    s = full.parameter(s)
    lo, hi = logarithm_interval(ratio)
    theta = -s*(lo+hi)/2
    real, imag, pr, pi = F(1), F(), F(1), F()
    for k in range(1, 65):
        pr, pi = -theta*pi/k, theta*pr/k
        real += pr; imag += pi
    radius = abs(theta)**65/math.factorial(65)+abs(s)*(hi-lo)/2
    return real, imag, radius


@lru_cache(maxsize=2)
def walsh_counts(kind):
    values = ratios(kind)
    groups = {r: tuple(i for i, v in enumerate(values) if v == r)
              for r in sorted(set(values)) if r != 1}
    return tuple((r, tuple(sum((-1)**((i & parity).bit_count()) for i in indices)
                          for parity in range(512))) for r, indices in groups.items())


def word(kind, s):
    s = full.parameter(s)
    rows = walsh_counts(kind)
    coefficients = [0.] * 512
    for ratio, counts in rows:
        lo, hi = logarithm_interval(ratio)
        angle = float(-s*(lo+hi)/2)
        coefficients = [a+angle*n/512 for a, n in zip(coefficients, counts, strict=True)]
    compiler = full.compiler
    tape = [compiler.pulse("I", coefficients[0], [0])]
    for pivot in range(9):
        count = 1 << pivot
        for index in range(count):
            gray = index ^ (index >> 1)
            tape.append(compiler.pulse("Z", coefficients[(1 << pivot) | gray], [pivot]))
            following = (index+1) % count
            change = gray ^ (following ^ (following >> 1))
            if change:
                tape.extend(compiler.cnot(change.bit_length()-1, pivot))
    if len(tape) != 2042:
        raise ArithmeticError("regional clock source schedule changed")
    return tape


@lru_cache(maxsize=8)
def _executed(kind, s):
    array = full.compiler.execute(word(kind, s), 9)
    return array.dtype.str, array.tobytes()


def matrix(kind, s):
    dtype, payload = _executed(kind, full.parameter(s))
    return np.frombuffer(payload, dtype=np.dtype(dtype)).reshape(512, 512)


def certify(kind, s, array=None):
    s = full.parameter(s)
    if array is None:
        array = matrix(kind, s)
    if type(array) is not np.ndarray:
        raise ValueError("an original base ndarray required")
    array = full.exact._array(array, (512, 512))
    # Dyadic target centers keep the all-entry comparison inexpensive. Their
    # outward enclosure charges BOTH real/imaginary rounding errors; executed
    # source arrays themselves remain their original binary values.
    scale = 1 << 128
    enclosed = {}
    for ratio in set(ratios(kind)):
        re, im, error = phase_enclosure(s, ratio)
        enclosed[ratio] = (F((re*scale).__floor__(), scale),
                           F((im*scale).__floor__(), scale), error+F(2, scale))
    phases = tuple(enclosed[ratio] for ratio in ratios(kind))
    squared = F()
    for i, row in enumerate(array):
        for j, value in enumerate(row):
            ar, ai = full.exact._component(value)
            re, im = phases[i][:2] if i == j else (F(), F())
            squared += (ar-re)**2+(ai-im)**2
    phase_error = full.exact._sqrt_interval(F(512))[1]*max(e for _, _, e in phases)
    eps = full.exact._sqrt_interval(squared)[1]+phase_error
    return {"operator_error_upper": eps, "half_diamond_upper": eps*(2+eps)/2,
            "channel_norm_upper": (1+eps)**2}


@lru_cache(maxsize=8)
def _check_word(kind, s):
    actual = certify(kind, s)
    full.corner.report._gate(actual["operator_error_upper"], WORD_ERROR_GATE,
                            "native completed-reference modular source word")


def checked_bound(kind, s):
    _check_word(kind, full.parameter(s))
    return {"operator_error_upper": WORD_ERROR_GATE,
            "half_diamond_upper": WORD_ERROR_GATE*(2+WORD_ERROR_GATE)/2,
            "channel_norm_upper": (1+WORD_ERROR_GATE)**2}


def budget(s=F(16), attempts=369, precision_bits=24):
    s = full.parameter(s)
    palm = preparation.budget(attempts, abs(s))
    cap = cut.cut_certificate(CAP, precision_bits=precision_bits)
    delta = palm["warmup_half_trace_error_upper"]
    eta = F(cap["preparation"]["trace_distance_upper"])
    t, u = 66*delta, 75*eta
    radius = t+u+t*u
    if radius >= 1:
        raise ValueError("combined preparation does not certify a faithful completed reference")
    source = checked_bound("clock", s)
    return {"warmup_native_attempts": attempts, "global_movement_herald_attempts": 1,
            "complement_native_attempts_after_success": cap["preparation"]["attempts"],
            "complement_precision_bits": precision_bits,
            "combined_reference_operator_sandwich_radius": radius,
            "actual_completed_reference_modular_half_diamond_upper": min(F(1), abs(s)*radius/(1-radius)),
            "source_to_actual_completed_modular_half_diamond_upper":
                source["half_diamond_upper"]+min(F(1), abs(s)*radius/(1-radius)),
            "source_then_cap_restrict_vs_actual_cap_modular_half_diamond_upper":
                source["half_diamond_upper"]+min(F(1), abs(s)*t/(1-t)),
            "complete_flagged_reference_half_trace_upper": delta+preparation.Q*eta,
            "complete_flagged_reference_and_source_half_trace_upper":
                delta+preparation.Q*eta+palm["herald_success_interval"][1]*source["half_diamond_upper"],
            "minimum_work_eigenvalue_lower":
                (1-radius)*preparation.Q/(75*palm["herald_success_interval"][1]),
            "expected_native_preparation_attempts_upper":
                palm["expected_attempts_with_full_warmup_repeated_each_trial_upper"]
                +cap["preparation"]["attempts"],
            "actual_cap_marginal_preserved_by_every_complement_attempt": True,
            "discarded_preparation_tapes_and_controller_required": True,
            "padding_retained_without_projection_or_renormalization": True,
            "complete_flagged_errors_apply_to_one_preparation_trial": True,
            "expected_cost_repeats_failed_Palm_trials_then_runs_one_complement_stage": True}


def mixer_word(sign, *, compact=False):
    if type(sign) is not int or sign not in (-1, 1):
        raise ValueError("an integer mixer sign required")
    target = 1 if compact else 8
    fanout = full.compiler.cnot(0, target)
    return [*fanout, full.compiler.pulse("Y", sign*math.pi/4, [0]), *fanout]


def mixer_bound(sign):
    """Four-dimensional complete-word certificate, lifted with spectators."""
    actual = full.compiler.execute(mixer_word(sign, compact=True), 2)
    lo, hi = full.exact._sqrt_interval(F(1, 2))
    middle, radius = (lo+hi)/2, (hi-lo)/2
    squared = F()
    for i in range(4):
        for j in range(4):
            expected = middle if i == j else sign*(1 if j & 1 else -1)*middle if i == (j ^ 3) else F()
            real, imag = full.exact._component(actual[i, j])
            squared += (real-expected)**2+imag*imag
    eps = full.exact._sqrt_interval(squared)[1]+full.exact._sqrt_interval(F(8))[1]*radius
    gate = full.MIXER_ERROR_GATE
    full.corner.report._gate(eps, gate, "native completed-cap Ramsey mixer")
    return {"half_diamond_upper": gate*(2+gate)/2, "channel_norm_upper": (1+gate)**2}


def witness(s=F(16)):
    """Existing source pulses and complete bit0 PVM read the cap-only gap.

    Codes 328 and 73 have the same outside configuration and cap charge one.
    The admitted two-code configuration preparation and projective bit read
    remain apparatus. A seven-pulse mixer joins the pair using bit0->bit8
    CNOT, Y, CNOT. Intermediate nonnative amplitudes remain in the register.
    No membership measurement or free arbitrary pair matrix is inserted.
    """
    s = full.parameter(s)
    codes = full.support()[0]
    left, right = codes.index(328), codes.index(73)
    states = regional.source()[0].states
    if (reference()[right]/reference()[left] != F(38, 37)
            or any(states[left][p] != states[right][p] for p in range(12) if p not in CAP)
            or sum(states[left][p] for p in CAP) != 1):
        raise ArithmeticError("native cap Ramsey pair no longer reads log(38/37)")
    initial = np.zeros(512, complex); initial[328] = 1
    prepared = full.compiler.execute(mixer_word(-1), 9, initial)
    phased = matrix("clock", s)@prepared
    analyzed = full.compiler.execute(mixer_word(1), 9, phased)
    probabilities = [float(np.vdot(analyzed[bit::2], analyzed[bit::2]).real) for bit in (0, 1)]
    real, _, error = phase_enclosure(s, F(38, 37))
    match = (1+real)/2
    interval = ((1+real-error)/2, (1+real+error)/2)
    mixer_bound(-1)
    delta = full.corner.composite_bound(checked_bound("clock", s), mixer_bound(1))
    half_l1 = (abs(probabilities[0]-float(match))+abs(probabilities[1]-float(1-match)))/2
    if half_l1 > float(delta+error/2):
        raise ArithmeticError("completed-cap Ramsey output exceeds its raw whole-word bound")
    outside = sum(abs(analyzed[i])**2 for i in range(512) if i not in set(codes))
    return {"initial_native_code": 328, "partner_native_code": 73,
            "same_outside_configuration_and_cap_charge_one": True,
            "native_modular_gap": "log(38/37)",
            "ideal_bit_zero_probability": "(1+cos(s log(38/37)))/2",
            "ideal_bit_zero_probability_interval": [full.compiled._outward(interval[0], upper=False),
                                                       full.compiled._outward(interval[1])],
            "complete_executed_bit_probabilities": probabilities,
            "complete_record_half_l1_error_diagnostic": half_l1,
            "whole_sequence_half_diamond_upper": full.compiled._outward(delta),
            "source_pulses": len(word("clock", s))+2*len(mixer_word(1)),
            "outside_native_probability_diagnostic": float(outside),
            "declared_initial_configuration_and_projective_bit_read_services": True,
            "initialization_and_readout_physical_costs_included": False,
            "padding_projected_or_outcomes_renormalized": False}


def run():
    sigma = reference()
    s = F(16)
    rows = []
    for kind in ("clock", "cocycle"):
        bound = checked_bound(kind, s)
        rows.append({"word": kind, "modular_parameter": str(s),
                     "source_pulses": len(word(kind, s)),
                     "source_word_sha256": hashlib.sha256(full.corner.canonical(word(kind, s)).encode()).hexdigest(),
                     "certified_bounds": preparation._compact(bound)})
    # A numerical diagnostic, separately charged from the exact cocycle proof.
    full.checked_bound(s)
    difference = np.linalg.norm(matrix("cocycle", s)@full.matrix(s)-matrix("clock", s), "fro")
    allowance = full.exact._sqrt_interval(F(512))[1]*WORD_ERROR_GATE*(3+WORD_ERROR_GATE)
    if difference > float(allowance):
        raise ArithmeticError("executed diagonal cocycle and completed word disagree")
    paths = ("model/native_modular_reference/completed_clock.py",
             "model/native_modular_reference/regional.py",
             "model/native_modular_reference/full_flow.py",
             "model/native_modular_reference/preparation.py",
             "model/native_equilibration/cut_completion.py",
             "model/native_equilibration/construction.py")
    return {"schema": "oph.native_modular_reference.completed_clock.v1", "all_pass": True,
            "source_sha256": {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
            "cap": list(CAP), "reference": "native complement completion of diag(degree/600)",
            "reference_weight_census": {str(w): sigma.count(w) for w in sorted(set(sigma))},
            "bare_cap_algebra_dimension": 53,
            "derived_cap_modular_gaps": {"charge_zero": "0", "charge_one": "log(38/37)",
                                          "charge_two": "log(5/4)"},
            "state_cocycle": "V_s=bar_sigma^(-is) sigma^(is); U_bar(s)=V_s U_Palm(s)",
            "source_words": rows,
            "cap_Ramsey_readout": witness(s),
            "native_complete_instrument_clock_covariance": branch_covariance_audit(),
            "executed_cocycle_composition_Frobenius_error_upper":
                full.compiled._outward(allowance),
            "finite_preparation_and_clock_budget": preparation._compact(budget(s)),
            "scope": {"same_native_prepared_reference_for_entropy_and_regional_clock": True,
                      "bare_cap_algebra_invariant_under_ideal_completed_reference_flow": True,
                      "finite_preparation_and_arbitrary_probe_references_covered": True,
                      "actual_global_prepared_modular_flow_exact_regional_invariance_claimed": False,
                      "nontrivial_cap_clock_requires_coherent_probe_and_source_controls": True,
                      "complement_schedule_and_discarded_records_remain_protocol_choices": True,
                      "passive_cut_repair_is_modular_evolution": False,
                      "physical_seconds_temperature_or_spacetime_identified": False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write", type=Path)
    group.add_argument("--check", type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        full.corner.compare_receipt(strict_json.load(args.check), result)
        print("Native entropy completion and shared regional modular clock verified.")
    elif args.write:
        with args.write.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(full.corner.canonical(result)+"\n")
    else:
        print(full.corner.canonical(result))


if __name__ == "__main__":
    main()
