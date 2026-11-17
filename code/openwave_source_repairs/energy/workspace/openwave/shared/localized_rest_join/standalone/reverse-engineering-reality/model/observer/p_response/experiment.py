#!/usr/bin/env python3
"""Conditional P -> electromagnetic coupling -> binary detector -> native feedback."""
from __future__ import annotations

import argparse
from collections import defaultdict
from decimal import Decimal as D, localcontext
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import random
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PAPER = "code/P_derivation/paper_math.py"
CERT = "code/P_derivation/runtime/p_interval_contraction_certificate_2026-07-14.json"
NATIVE = "model/observer/native.py"
# Published inputs are pinned independently of this module's output receipt.
PINS = {
    PAPER: "bda031185dca55e54f4c6c0c020e9fac2e2ba5db5193998135f294898af620b8",
    CERT: "20ab213f836fbe95149d4b40ba57aa144038580a0acb899322f4b72d963c2ad6",
    NATIVE: "dd7c87bcda73f55990690effd0633020ff05b9131205d2314793329f8e0c811d",
}
MODE = "thomson_structured_running"
EXPOSURE = D(1)  # closed-loop integral of the normalized external potential
LOOP_CELLS = D(1)  # aperture normalization, not the size of the apparatus
SAMPLE_COUNT = 200000
SAMPLE_SEED = 20261012


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def sources():
    for path, expected in PINS.items():
        if hashlib.sha256((ROOT / path).read_bytes()).hexdigest() != expected:
            raise ValueError(f"Source pin mismatch: {path}")
    from model.observer import native
    return load_module(PAPER, "oph_p_response_paper"), native


def sin_cos(x):
    """Decimal series for our bounded phases; no host-math dependency."""
    with localcontext() as ctx:
        ctx.prec = 60
        s = term_s = x
        c = term_c = D(1)
        for k in range(1, 180):
            term_s *= -x*x / D((2*k)*(2*k+1))
            term_c *= -x*x / D((2*k-1)*(2*k))
            s += term_s
            c += term_c
            if abs(term_s) + abs(term_c) < D("1e-58"):
                return +s, +c
        raise ArithmeticError("Series did not converge")


def show(x):
    return format(x, ".24g")


def adapter_probability(current, outcome, g, exposure=EXPOSURE):
    if current not in (-1, 0, 1) or outcome not in (-1, 1):
        raise ValueError("Invalid current or detector symbol")
    sine, _ = sin_cos(g * D(current) * exposure)
    return (D(1) + D(outcome) * sine) / 2


def decoder_seam(detector_bit):
    """The feedback actuator sees only the new detector bit, never native j."""
    if detector_bit not in (-1, 1):
        raise ValueError("Invalid detector bit")
    return (0, 5) if detector_bit == 1 else (0, 1)


def controlled_native_step(native, loads, detector_bit, coin):
    """Preserve the committed detector bit while it selects a native repair."""
    new_loads, native_record = native.native(loads, decoder_seam(detector_bit), coin)
    return new_loads, detector_bit, native_record


def coupling(ctx, p):
    d10 = ctx.build_d10_from_p(p)
    alpha, inverse, running = ctx.alpha_external_from_d10(d10, MODE)
    g = (D(4) * ctx.pi * alpha).sqrt()
    return {"p": p, "alpha": alpha, "inverse": inverse, "g": g,
            "alpha_u": d10.alpha_u, "mz_anchor_inverse": d10.alpha_em_inv_mz,
            "thomson_shift": running["total_delta_alpha_inv"]}


def histories(native, states, g, exposure=EXPOSURE):
    """Every source state, coin, detector bit, and second coin is retained."""
    binary = [x for x in states if max(x) == 1]
    index = {x: k for k, x in enumerate(states)}
    rows = []
    mass = occupation = correlation = outcome_mean = next_current = D(0)
    for x in binary:
        for c1 in (0, 1):
            after_first, r1 = native.native(x, (0, 1), c1)
            for bit in (-1, 1):
                chance = adapter_probability(r1[1], bit, g, exposure)
                # The stable detector bit is committed before the actuator runs.
                committed_bit = bit
                seam = decoder_seam(committed_bit)
                for c2 in (0, 1):
                    final, retained_bit, r2 = controlled_native_step(native, after_first, committed_bit, c2)
                    weight = chance / D(264)
                    assert retained_bit == committed_bit
                    mass += weight
                    occupation += weight * D(final[0])
                    correlation += weight * D(bit * r1[1])
                    outcome_mean += weight * D(bit)
                    next_current += weight * D(r2[1])
                    rows.append({"initial_source": index[x], "first_native_record": list(r1),
                                 "detector_bit": bit, "second_seam": list(seam),
                                 "second_native_record": list(r2), "final_source": index[final],
                                 "detector_probability": show(chance), "history_probability": show(weight)})
    return rows, {"mass": mass, "occupation": occupation,
                  "correlation": correlation, "outcome_mean": outcome_mean,
                  "next_current": next_current}


def sampled_execution(native, states, g):
    """Execute fixed-seed sequential trials; the actuator receives one bit."""
    rng = random.Random(SAMPLE_SEED)
    binary = [x for x in states if max(x) == 1]
    plus = {j: float(adapter_probability(j, 1, g)) for j in (-1, 0, 1)}
    counts = defaultdict(int)
    total = squares = 0
    for _ in range(SAMPLE_COUNT):
        x = binary[rng.randrange(len(binary))]
        x, r1 = native.native(x, (0, 1), rng.randrange(2))
        # Router j sets the detector law; the later controller only gets bit.
        stable_bit = 1 if rng.random() < plus[r1[1]] else -1
        x, retained_bit, r2 = controlled_native_step(native, x, stable_bit, rng.randrange(2))
        assert retained_bit == stable_bit
        counts[r1, stable_bit, r2] += 1
        total += r2[1]
        squares += r2[1]**2
    mean = total/SAMPLE_COUNT
    sample_variance = (squares - total*total/SAMPLE_COUNT)/(SAMPLE_COUNT-1)
    return {"seed": SAMPLE_SEED, "trials": SAMPLE_COUNT,
            "mean_second_native_current": mean,
            "standard_error_of_mean": math.sqrt(sample_variance/SAMPLE_COUNT),
            "sum_second_native_current": total, "sum_squared_second_native_current": squares,
            "terminal_word_counts": [{"r1": list(r1), "bit": bit, "r2": list(r2), "count": n}
                                     for (r1, bit, r2), n in sorted(counts.items())]}


def run():
    paper, native = sources()
    cert = json.loads((ROOT / CERT).read_text())["modes"][MODE]
    states = native.states()
    checks = {}

    def check(name, assertion):
        checks[name] = bool(assertion)
        if not assertion:
            raise AssertionError(name)

    with localcontext() as work:
        work.prec = 60
        box = cert["certified_enclosure"]["P"]
        p = (D(box["lo"]) + D(box["hi"])) / 2
        # The existing declared source function is actually evaluated here.
        # Moderate cutoffs are enough at this branch; the independent certificate
        # endpoint checks the resulting scalar, without rerunning that certificate.
        ctx = paper.PaperMathContext(precision=40, su2_cutoff=20, su3_cutoff=15)
        point = coupling(ctx, p)
        alpha_box = cert["certified_enclosure"]["alpha"]
        alpha_reference = (D(alpha_box["lo"]) + D(alpha_box["hi"])) / 2
        source_error = abs(point["alpha"] - alpha_reference)
        check("executed_source_map_matches_pinned_closure_branch", source_error < D("1e-13"))
        closure_residual = p - ctx.phi - ctx.sqrt_pi * point["alpha"]
        check("accepted_outer_closure_numerically_recovered", abs(closure_residual) < D("1e-13"))
        g = point["g"]
        sine, cosine = sin_cos(g * EXPOSURE)
        check("nontrivial_unsaturated_detector", D(0) < sine < D(1) and D(0) < g < ctx.pi / 2)

        # Exact conservation/gauge phase exponent identity, for all 78 states.
        chi = [Q(k*k - 4*k + 1, 7) for k in range(12)]
        branch_checks = 0
        for x in states:
            for seam in ((0, 1), (0, 5)):
                for coin in (0, 1):
                    y, r = native.native(x, seam, coin)
                    j = r[1]
                    assert sum(y) == sum(x)
                    assert y[seam[0]] - x[seam[0]] == -j
                    assert y[seam[1]] - x[seam[1]] == j
                    assert sum(chi[k] * (y[k] - x[k]) for k in range(12)) == j*(chi[seam[1]]-chi[seam[0]])
                    branch_checks += 1
        check("all_312_native_branches_have_conserved_charge_and_covariant_transfer", branch_checks == 312)

        # Diagonal effects E_y have eigenvalues p(y|j). The corresponding
        # Ramsey amplitudes m_y(j) obey completeness analytically for all g,chi.
        for j in (-1, 0, 1):
            probs = [adapter_probability(j, y, g) for y in (-1, 1)]
            assert abs(sum(probs) - 1) < D("1e-50")
            assert all(0 <= q <= 1 for q in probs)
            assert adapter_probability(j, 1, g, D(0)) == D("0.5")
            assert adapter_probability(j, 1, D(0)) == D("0.5")
            assert abs(adapter_probability(j, 1, g) - adapter_probability(-j, -1, g)) < D("1e-50")
        check("binary_povm_complete_positive_zero_and_reversal_limits", True)
        eps = D("1e-8")
        finite_slope = (adapter_probability(1, 1, g, eps) - D("0.5")) / eps
        check("weak_exposure_slope_is_charge_over_two", abs(finite_slope - g/2) < D("1e-16"))
        perfect_exposure = ctx.pi / (2*g)
        check("calibrated_pi_over_two_limit_discriminates_signed_current",
              abs(adapter_probability(1, 1, g, perfect_exposure) - 1) < D("1e-48")
              and abs(adapter_probability(-1, 1, g, perfect_exposure)) < D("1e-48"))

        rows, stats = histories(native, states, g)
        expected = D(1)/6 - D(5)*sine/264
        check("all_528_histories_total_one_without_postselection",
              len(rows) == 528 and abs(stats["mass"] - 1) < D("1e-50"))
        check("binary_detector_drives_p_dependent_native_feedback",
              abs(stats["occupation"] - expected) < D("1e-50")
              and abs(stats["next_current"] - D(5)*sine/264) < D("1e-50"))
        check("detector_bit_marginal_half_and_record_correlation_nonzero",
              abs(stats["outcome_mean"]) < D("1e-50")
              and abs(stats["correlation"] - D(5)*sine/33) < D("1e-50"))
        _, disconnected = histories(native, states, g, exposure=D(0))
        check("zero_exposure_disconnect_restores_one_sixth",
              abs(disconnected["occupation"] - D(1)/6) < D("1e-50"))
        # The actual bit marginal is exactly fair. Independent matched shuffling
        # therefore equals the zero-exposure kernel (same half/half seam count).
        check("same_marginal_bit_shuffle_matches_disconnect", abs(stats["outcome_mean"]) < D("1e-50"))

        variants = []
        for ratio in (D("0.999"), D("1.001")):
            other = coupling(ctx, p * ratio)
            other_sine, _ = sin_cos(other["g"] * EXPOSURE)
            _, other_stats = histories(native, states, other["g"])
            assert abs(other_stats["occupation"] - (D(1)/6 - D(5)*other_sine/264)) < D("1e-50")
            variants.append({"P_multiplier": str(ratio), "P": show(other["p"]),
                             "alpha": show(other["alpha"]), "g": show(other["g"]),
                             "detector_plus_given_j_plus": show((1+other_sine)/2),
                             "final_occupation_0": show(other_stats["occupation"]),
                             "scope": "off-closure parameter sensitivity, not a second closure solution"})
        check("fixed_exposure_parameter_controls_change_actual_response",
              variants[0]["detector_plus_given_j_plus"] != variants[1]["detector_plus_given_j_plus"])
        words = {(tuple(r["first_native_record"]), r["detector_bit"], tuple(r["second_native_record"])) for r in rows}
        first_words = {w[0] for w in words}
        detected_words = {w[:2] for w in words}
        check("retained_audit_support_counts", len(words) == 24 and len(first_words) == 4 and len(detected_words) == 8)
        samples = sampled_execution(native, states, g)
        discrepancy = abs(samples["mean_second_native_current"] - float(stats["next_current"]))
        check("fixed_seed_sequential_execution_matches_analytic_mean_within_six_se",
              discrepancy < 6*samples["standard_error_of_mean"] and sum(r["count"] for r in samples["terminal_word_counts"]) == SAMPLE_COUNT)
        return {"schema": "oph.p_response.v1", "all_pass": all(checks.values()), "checks": checks,
                "assumption": "P and N closure outputs are accepted for this conditional hardware design; N is not fitted or used to alter the local coupling",
                "source_pins": PINS, "source_branch": {
                    "mode": MODE, "scale": "declared low-energy Thomson structured-running branch",
                    "P": show(p), "P_imported_interval": box,
                    "alpha": show(point["alpha"]), "inverse_alpha": show(point["inverse"]),
                    "alpha_U": show(point["alpha_u"]), "mz_anchor_inverse": show(point["mz_anchor_inverse"]),
                    "low_energy_inverse_shift": show(point["thomson_shift"]),
                    "g_sqrt_4pi_alpha": show(g), "recomputed_alpha_difference_from_certificate_midpoint": show(source_error),
                    "outer_closure_residual": show(closure_residual),
                    "point_solver": {"precision": 40, "su2_cutoff": 20, "su3_cutoff": 15},
                    "numerical_scope": "deterministic point calculation checked against existing enclosure; this packet supplies no new interval or physical-endpoint certificate"},
                "designed_hardware": {"charge_map": "Q_v=g(P) N_v; transferred charge=g(P) j",
                    "gauge_scope": "branch charge-transfer covariance; coherent record/link transformation must be included for a gauge-invariant enlarged isometry",
                    "probe": "unit-charge path qubit, prepare +X, route loop orientation by j, measure Y",
                    "exposure": str(EXPOSURE), "exposure_definition": "closed-loop integral of normalized externally prepared A; phi=g*j*exposure",
                    "aperture_area_over_ellstar_squared": show(LOOP_CELLS*p),
                    "aperture_canonical_cells": str(LOOP_CELLS),
                    "aperture_scope": "one cell is a chart normalization; arbitrary feasible K_loop scales aperture and inversely scales field at fixed exposure, not the total device size",
                    "normalized_field_density_B_ellstar_squared_for_this_aperture": show(EXPOSURE/(LOOP_CELLS*p)),
                    "feedback": "new detector bit +1 selects (0,5), -1 selects (0,1); it never reads j",
                    "retained_records": "native r1, new binary detector bit, native r2; all separately retained",
                    "actuator_access": "phase and one stable binary detector bit only; two alternative bit values",
                    "rich_audit_reachable_terminal_words": len(words),
                    "rich_audit_prefix_count": 1+len(first_words)+len(detected_words)+len(words),
                    "rich_audit_prefix_encoding_bits": 6,
                    "one_fixed_tape_encoding": "6 x 2 x 6 = 72 terminal record labels, excluding source, phase, blank sectors, probe and control",
                    "physical_readiness": "conditional designed interface; apparatus initialization, routing, measurement and isolation are assumed available, not fabricated"},
                "response": {"phase_for_j_plus": show(g*EXPOSURE), "sine_visibility": show(sine),
                    "plus_probability_by_current": {str(j): show(adapter_probability(j,1,g)) for j in (-1,0,1)},
                    "balanced_signed_current_discrimination_success": show((1+sine)/2),
                    "source_record_detector_bit_correlation": show(stats["correlation"]),
                    "mean_second_native_current": show(stats["next_current"]),
                    "mean_second_transferred_charge": show(g*stats["next_current"]),
                    "charge_weak_exposure_law": "5*g(P)^2*exposure/264 = (5*pi/66)*alpha(P)*exposure, plus O(exposure^3)",
                    "detector_plus_marginal": "0.5", "final_occupation_0": show(stats["occupation"]),
                    "disconnected_or_independent_shuffled_bit_occupation": show(D(1)/6),
                    "feedback_difference": show(-D(5)*sine/264),
                    "feedback_formula": "1/6 - 5*sin(g(P)*exposure)/264",
                    "information_scope": "a noisy electromagnetic transcription of the known native current, not new information about it or a learned predictor"},
                "sensitivity_controls": variants, "sampled_execution": samples, "histories": rows}


def canonical(result):
    return json.dumps(result, indent=2, sort_keys=True) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write", type=Path, help="write a full receipt to a new file")
    group.add_argument("--check", type=Path, help="compare the full receipt with an existing file")
    args = parser.parse_args()
    result = run()
    content = canonical(result)
    if args.write:
        with args.write.open("x") as stream:
            stream.write(content)
    if args.check and args.check.read_text() != content:
        raise SystemExit("Result mismatch")
    print(json.dumps({"all_pass": result["all_pass"], "checks": len(result["checks"]),
                      "result_sha256": hashlib.sha256(content.encode()).hexdigest(),
                      "response": result["response"]}))


if __name__ == "__main__":
    main()
