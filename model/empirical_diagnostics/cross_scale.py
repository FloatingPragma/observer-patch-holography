"""Fixed spectral coupling, its conditional primordial map, and joint scope.

Exact rational enclosures are sealed before reading comparison data. This is
a retrospective test of declared identifications, not a new physical bridge.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
from math import isqrt
from pathlib import Path

from tools import strict_json
from model.predictions import endpoint_audit, shared_parameter
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "cross_scale_results.json"
ALPHA = "model/observer_particles/results/entangled_alpha.json"
CHIRAL = "model/common_source/chiral_candidate/results.json"
QUANTUM = "model/scalar_observations/experimental_readout/executable_results.json"
EXPRESSION = "(79+48 sqrt(2)+42 sqrt(3)+3 sqrt(6))/31104"


def canonical(value):
    return shared_parameter.canonical(value)


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sqrt_bounds(value, digits=60):
    """Enclose sqrt of a nonnegative rational by exact integer squaring."""
    if not isinstance(value, F) or value < 0 or type(digits) is not int or digits < 1:
        raise ValueError("nonnegative rational and positive integer precision required")
    scale = 10**digits
    low = isqrt(value.numerator * scale**2 // value.denominator)
    lower = F(low, scale)
    upper = lower if lower**2 == value else F(low + 1, scale)
    assert lower**2 <= value <= upper**2
    return lower, upper


def arctan_reciprocal_bounds(denominator, terms):
    """Alternating-series enclosure for atan(1/denominator)."""
    if type(denominator) is not int or denominator <= 1 or type(terms) is not int or terms < 1:
        raise ValueError("integer denominator above one and positive term count required")
    partial = sum((F((-1)**k, (2*k + 1)*denominator**(2*k + 1))
                   for k in range(terms)), F())
    next_partial = partial + F((-1)**terms, (2*terms + 1)*denominator**(2*terms + 1))
    return min(partial, next_partial), max(partial, next_partial)


def pi_bounds():
    # Machin's identity pi = 16 atan(1/5) - 4 atan(1/239).
    a, b = arctan_reciprocal_bounds(5, 48), arctan_reciprocal_bounds(239, 16)
    return 16*a[0] - 4*b[1], 16*a[1] - 4*b[0]


def spectral_prediction():
    """No observation or receipt is opened by this forward calculation."""
    roots = {n: sqrt_bounds(F(n)) for n in (2, 3, 5, 6)}
    alpha = tuple((79 + 48*roots[2][i] + 42*roots[3][i] + 3*roots[6][i])/31104
                  for i in (0, 1))
    pi = pi_bounds()
    sqrt_pi = sqrt_bounds(pi[0])[0], sqrt_bounds(pi[1])[1]
    phi = tuple((1 + roots[5][i])/2 for i in (0, 1))
    pixel = tuple(phi[i] + alpha[i]*sqrt_pi[i] for i in (0, 1))
    tilt = shared_parameter.tilt_interval(*pixel)
    return {
        "spectral_expression": EXPRESSION,
        "alpha_interval": list(map(str, alpha)),
        "inverse_alpha_interval": list(map(str, (1/alpha[1], 1/alpha[0]))),
        "P_interval": list(map(str, pixel)),
        "scalar_tilt_interval": list(map(str, tilt)),
        "declared_maps": ["P=(1+sqrt(5))/2+sqrt(pi)*alpha", "n_s=1-P/48"],
        "method": "Exact rational square-root bounds and alternating Machin-series bounds",
        "arithmetic_error_is_physical_theory_error": False,
        "observation_inputs": [],
    }


def _interval(pair):
    return endpoint_audit._interval(*pair)


def comparison(prediction, target):
    centre, uncertainty, _ = endpoint_audit._target(target)
    low, high = map(F, prediction["inverse_alpha_interval"])
    return {
        "target": target,
        "inverse_alpha_decimal_enclosure": _interval((low, high)),
        "signed_difference_over_reported_standard_uncertainty": _interval(
            ((low-centre)/uncertainty, (high-centre)/uncertainty)),
        "within_one_reported_standard_uncertainty": centre-uncertainty <= low <= high <= centre+uncertainty,
        "interpretation": "Retrospective scalar agreement; reported standard uncertainty is not a model-error bound or model probability.",
    }


def joint_status(*, alpha_identified, same_cell_proved, chiral, quantum):
    """Fail closed on actual physical predicates, even if a scalar agrees."""
    if type(alpha_identified) is not bool or type(same_cell_proved) is not bool:
        raise ValueError("explicit physical-identification predicates required")
    verdict = chiral["verdict"]
    scope = quantum["scope"]
    missing = []
    if not alpha_identified:
        missing.append("spectral_scalar_to_Thomson_photon_residue")
    if not same_cell_proved:
        missing.append("same_cell_P_alpha_identification")
    for key in ("same_reconstructed_native_spacetime_field_history_proved",
                "quantum_gauge_and_Higgs_dynamics_executed",
                "observed_standard_model_mass_spectrum_derived"):
        if verdict[key] is not True:
            missing.append(key)
    if scope["preparation_or_action_selected_by_OPH"] is not True:
        missing.append("oscillator_apparatus_action_and_preparation_from_fixed_repair")
    if scope["physical_clock_calibrated"] is not True:
        missing.append("physical_clock_for_quantum_comparison")
    # None of the consumed packets supplies these two contracts. Numerical
    # precision or posterior coverage must not silently stand in for them.
    missing.extend(("native_curvature_to_CMB_transfer_and_background",
                    "common_physical_action_and_complete_comparison_error_contract"))
    comparisons = chiral["separate_empirical_comparison"]["comparisons"]
    if not comparisons or len({row["prescription"] for row in comparisons}) != len(comparisons):
        raise ValueError("distinct nonempty declared mass comparisons required")
    failures = [row["prescription"] for row in comparisons if row["both_ratios_pass"] is False]
    if any(type(row["both_ratios_pass"]) is not bool for row in comparisons):
        raise ValueError("explicit mass-comparison predicates required")
    for row in comparisons:
        intersections = row["ratio_intersects_comparison_windows"]
        if (not isinstance(intersections, list) or len(intersections) != 2
                or any(type(value) is not bool for value in intersections)
                or row["both_ratios_pass"] != all(intersections)):
            raise ValueError("mass verdict disagrees with its two interval comparisons")
    return {
        "status": "not_evaluable_as_one_complete_physical_construction",
        "missing_contracts": missing,
        "declared_chiral_prescriptions": [row["prescription"] for row in comparisons],
        "declared_chiral_prescriptions_passing_both_lepton_ratios": [
            row["prescription"] for row in comparisons if row["both_ratios_pass"] is True],
        "declared_chiral_prescriptions_failing_both_ratio_requirement": failures,
        "joint_likelihood_or_significance": None,
        "all_public_datasets_reproduced": False,
        "source_only_measured_universe_established": False,
    }


def build(*, root=ROOT):
    root = Path(root)
    prediction = spectral_prediction()
    seal = hashlib.sha256(canonical(prediction).encode()).hexdigest()
    alpha = strict_json.loads((root/ALPHA).read_text(encoding="utf-8"))
    if alpha["spectral_expression"] != EXPRESSION or alpha["empirical_inputs"] != []:
        raise ValueError("the retained spectral construction or its input contract changed")
    for key, field in (("alpha_interval", "alpha"), ("inverse_alpha_interval", "inverse_alpha"),
                       ("P_interval", "P_after_declared_same_cell_attachment")):
        lower, upper = map(F, prediction[key])
        if not lower <= F(alpha[field]) <= upper:
            raise ValueError("retained spectral scalar disagrees with independent enclosure")
    target = endpoint_audit.target_fixture()
    coupling = comparison(prediction, target)
    cosmology = shared_parameter.compare_projection(prediction, root=root)
    for chain in cosmology["per_chain"]:
        # These are remote archive-member hashes, not repository source maps.
        chain["archive_member_sha256"] = chain.pop("source_sha256")
    cosmology["physical_comparison_contract"] = {
        "quantity": "primordial scalar power-law spectral index n_s",
        "pivot_Mpc_inverse": "0.05",
        "cosmological_model": "baseline six-parameter flat Lambda-CDM",
        "data_combination": "P-ACT-LB: Planck, full multifrequency ACT DR6, ACT/Planck lensing, DESI 2024 DR1 BAO",
        "archive_product": "p-actbase-l-b_lcdm_camb",
        "ns_prior": {"family": "uniform", "bounds": ["0.85", "1.1"]},
        "supplied_physics": "standard radiation transfer, recombination, foreground, calibration and nuisance models",
        "posterior_marginal_not_joint_model_likelihood": True,
        "documentation": "code/cosmology/tilt_posterior/README.md",
        "primary_parameter_paper": "https://arxiv.org/html/2503.14452v1",
        "native_primordial_embedding_derived": False,
    }
    chiral = strict_json.loads((root/CHIRAL).read_text(encoding="utf-8"))
    quantum = strict_json.loads((root/QUANTUM).read_text(encoding="utf-8"))
    if (chiral["schema"] != "oph.chiral_candidate.v1"
            or quantum["schema"] != "oph.scalar_observations.executable_experimental_readout.v1"):
        raise ValueError("comparison receipt schema changed")
    scope = alpha["physical_scope"]
    status = joint_status(alpha_identified=scope["Thomson_photon_residue_identified"],
                          same_cell_proved=scope["same_cell_P_alpha_bridge_proved"],
                          chiral=chiral, quantum=quantum)
    if hashlib.sha256(canonical(prediction).encode()).hexdigest() != seal:
        raise ArithmeticError("forward prediction changed after opening observations")
    source_paths = (ALPHA, CHIRAL, QUANTUM, shared_parameter.PROJECTION,
                    "code/cosmology/tilt_posterior/README.md",
                    "model/empirical_diagnostics/cross_scale.py",
                    "model/predictions/shared_parameter.py", "model/predictions/endpoint_audit.py")
    return {
        "schema": "oph.empirical_diagnostics.fixed_spectral_cross_scale.v1",
        "classification": "retrospective conditional cross-scale comparison",
        "attribution": alpha["attribution"],
        "prediction": prediction, "prediction_sha256_before_comparison": seal,
        "coupling": coupling, "cosmology": cosmology,
        "chiral_mass_comparison": chiral["separate_empirical_comparison"],
        "quantum_comparison_scope": quantum["scope"],
        "joint_status": status,
        "custody": {"input_sha256": bound_pins({path: _sha(root/path) for path in source_paths}),
                    "upstream_histories_reexecuted_here": False,
                    "exact_spectral_enclosure_recomputed_here": True,
                    "ACT_weighted_projection_recomputed_here": True},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--check", type=Path)
    action.add_argument("--write", type=Path)
    args = parser.parse_args()
    result = build()
    if args.check:
        if strict_json.loads(args.check.read_text(encoding="utf-8")) != result:
            raise SystemExit("fixed spectral cross-scale receipt drift")
    else:
        with args.write.open("x", encoding="utf-8") as stream:
            stream.write(canonical(result))
    print(canonical({"inverse_alpha": result["coupling"]["inverse_alpha_decimal_enclosure"],
                     "ACT_tail_weights": result["cosmology"]["forward_tail_weight_interval"],
                     "joint_status": result["joint_status"]["status"]}), end="")


if __name__ == "__main__":
    main()
