"""Recompute hierarchy obstructions without changing the retained action.

The source spectrum and algebra checks are independent of empirical data.
Benchmark windows enter only the final exclusion and correction-size tests.
This audits retained quantum records; it does not replay the source evolution.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import sympy as s

from ..chiral_candidate import action, history, run as candidate_run
from tools import strict_json


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BASELINE = HERE.parent / "chiral_candidate/results.json"
NOTE = ROOT / "docs/research/OPH_PARTICLE_MASS_RESPONSE.md"
REPLAY_ABSOLUTE_TOLERANCE = 5e-11
REPLAY_RELATIVE_TOLERANCE = 2e-9


def compare_replay(retained, current, path="report"):
    """Compare discrete and exact claims strictly, floating diagnostics portably."""
    if type(retained) is not type(current):
        raise ValueError(f"replay type changed: {path}")
    if isinstance(retained, dict):
        if retained.keys() != current.keys():
            raise ValueError(f"replay keys changed: {path}")
        for key in retained:
            compare_replay(retained[key], current[key], f"{path}.{key}")
    elif isinstance(retained, list):
        if len(retained) != len(current):
            raise ValueError(f"replay length changed: {path}")
        for index, (old, new) in enumerate(zip(retained, current, strict=True)):
            compare_replay(old, new, f"{path}[{index}]")
    elif isinstance(retained, float):
        if not (math.isfinite(retained) and math.isfinite(current)
                and math.isclose(retained, current,
                                 rel_tol=REPLAY_RELATIVE_TOLERANCE,
                                 abs_tol=REPLAY_ABSOLUTE_TOLERANCE)):
            raise ValueError(f"floating replay diagnostic changed: {path}")
    elif retained != current:
        raise ValueError(f"replay claim or exact value changed: {path}")


def _zero(matrix):
    return all(s.simplify(x) == 0 for x in matrix)


def _scalar(value):
    value = s.simplify(value)
    return {"exact": str(value), "float_diagnostic": float(value)}


def spectral_data():
    """Actual family feature and its exact rank-one spectral projectors."""
    feature = action.family_feature()
    if feature != feature.H:
        raise ArithmeticError("a Hermitian family feature is required")
    spectrum = feature.eigenvals()
    if len(spectrum) != 3 or set(spectrum.values()) != {1}:
        raise ArithmeticError("three simple family eigenvalues are required")
    values = tuple(sorted(spectrum, key=float))
    if not all(x.is_positive for x in values):
        raise ArithmeticError("positive source eigenvalues are required")
    projectors = []
    for i, value in enumerate(values):
        projector = s.eye(3)
        for j, other in enumerate(values):
            if i != j:
                projector = projector * (feature-other*s.eye(3))/(value-other)
        projectors.append(s.ImmutableMatrix(projector.applyfunc(s.simplify)))
    if not _zero(sum(projectors, s.zeros(3))-s.eye(3)):
        raise ArithmeticError("incomplete family spectral resolution")
    for i, p in enumerate(projectors):
        if p != p.H or s.trace(p) != 1:
            raise ArithmeticError("nonorthogonal family projector")
        for j, q in enumerate(projectors):
            if not _zero(p*q-(p if i == j else s.zeros(3))):
                raise ArithmeticError("family spectral projector identity failed")
    return feature, values, tuple(projectors)


def source_algebra():
    """Exact affine obstruction and complete single-feature function algebra."""
    feature, values, projectors = spectral_data()
    phi = (1+s.sqrt(5))/2
    spacing = s.simplify((values[2]-values[0])/(values[1]-values[0]))
    if s.simplify(spacing-phi**2) != 0:
        raise ArithmeticError("native affine spacing relation changed")
    # u,v are unconstrained complex symbols. This is not a real-sign test.
    u, v = s.symbols("u v")
    z = [u+v*x for x in values]
    if s.simplify(z[2]-phi**2*z[1]+phi*z[0]) != 0:
        raise ArithmeticError("complex affine eigenvalue relation failed")
    # Every function on this simple finite spectrum is a quadratic polynomial.
    # These are symbolic targets, not measured masses or fitted coefficients.
    targets = s.symbols("target0:3")
    interpolant = sum((t*p for t, p in zip(targets, projectors, strict=True)), s.zeros(3))
    for t, p in zip(targets, projectors, strict=True):
        if not _zero(interpolant*p-t*p):
            raise ArithmeticError("symbolic spectral interpolation failed")
    ys = action.yukawa_matrices()
    for left in ys.values():
        for right in ys.values():
            a, b = left*left.H, right*right.H
            if not _zero(a*b-b*a):
                raise ArithmeticError("the retained trial is no longer aligned")
    return {
        "eigenvalues_ascending": [_scalar(x) for x in values],
        "linear_mass_ratios": [_scalar(x/values[0]) for x in values],
        "source_spread": _scalar(values[-1]/values[0]),
        "affine_eigenvalue_spacing": _scalar(spacing),
        "complex_affine_relation": "z3 - phi^2*z2 + phi*z1 = 0; phi=(1+sqrt(5))/2",
        "sorted_affine_singular_value_bound": "m_max <= phi^2*m_middle + phi*m_min",
        "affine_scope": "Y=u*I+v*F, arbitrary complex u,v; canonical common kinetic form and Higgs factor",
        "all_phases_and_eigenaxis_assignments_covered": True,
        "spectral_projectors_exact": [[[str(x) for x in row] for row in p.tolist()] for p in projectors],
        "single_feature_function_algebra_dimension": 3,
        "arbitrary_symbolic_values_interpolate_at_degree_at_most_two": True,
        "current_Y_Ydagger_commutators_zero_exact": True,
        "all_functions_of_F_share_family_axes": True,
        "functional_maps_alone_produce_nontrivial_CKM_mixing": False,
        "empirical_inputs": [],
    }


def correction_bounds(spread, required_ratio):
    """Necessary, not sufficient, sizes for two specified modifications.

    ZL/ZR are positive definite constant kinetic forms. Additive errors are
    measured in spectral norm relative to the original smallest singular value.
    Momentum-dependent self-energies and extra states are not covered.
    """
    spread, required_ratio = s.sympify(spread), s.sympify(required_ratio)
    if not (spread >= 1 and required_ratio > spread):
        raise ValueError("require target spread greater than original spread >= 1")
    return {
        "required_sqrt_kinetic_condition_product_at_least": _scalar(required_ratio/spread),
        "required_kinetic_condition_product_at_least": _scalar((required_ratio/spread)**2),
        "required_epsilon_if_both_Z_within_epsilon_of_identity_at_least": _scalar(
            (required_ratio-spread)/(required_ratio+spread)),
        "required_additive_norm_over_original_minimum_at_least": _scalar(
            (required_ratio-spread)/(required_ratio+1)),
        "kinetic_inequality": "kappa(Mcanonical) <= kappa(M)*sqrt(kappa(ZL)*kappa(ZR))",
        "additive_inequality": "kappa(M+dM) <= (kappa(M)+delta)/(1-delta), delta=norm(dM)/s_min(M)<1",
        "arbitrary_strong_or_momentum_dependent_dynamics_excluded": False,
        "mere_coordinate_change_alters_poles": False,
    }


def benchmark_exclusion(source, retained):
    # Reuse the original broad benchmark; it cannot affect source_algebra().
    windows = retained["separate_empirical_comparison"]["conservative_ratio_windows"]
    (mu_lo, mu_hi), (tau_lo, tau_hi) = [[s.Rational(x) for x in row] for row in windows]
    phi = (1+s.sqrt(5))/2
    bound = s.simplify(phi**2*mu_hi+phi)
    if not (tau_lo-bound).is_positive:
        raise ArithmeticError("the empirical window no longer excludes affine laws")
    spread = s.sympify(source["source_spread"]["exact"])
    return {
        "scope": "empirical comparison only; no source coefficient changed",
        "charged_lepton_ratio_windows": windows,
        "benchmark_source": retained["separate_empirical_comparison"]["source"],
        "maximum_affine_heavy_ratio_at_admitted_middle_ratio": _scalar(bound),
        "heavy_ratio_lower_bound": str(tau_lo),
        "affine_class_excluded_exact": True,
        "constant_kinetic_or_additive_correction_requirements": correction_bounds(spread, tau_lo),
    }


def record_audit(retained):
    report = retained["continuing_source_history"]
    records = {tuple(map(int, key)): value for key, value in report["source_record_probabilities"].items()}
    decoded = history.infer_ratios(records, report["clock_calibration"], report["phase_window_contract"])
    gaps = sorted(decoded["native_family_order_gap_over_clock"])
    _, values, _ = spectral_data()
    coefficient = action.source_coefficients()["J"]
    target = [float(coefficient*x) for x in values]  # v=Omega in this frozen action.
    residual = max(abs(x-y) for x, y in zip(gaps, target, strict=True))
    if residual > 1e-10:
        raise ArithmeticError("records no longer recover the trial action's rest gaps")
    return {
        "record_count": len(records),
        "decoded_gaps_over_clock_ascending": gaps,
        "decoded_ascending_ratios": [x/gaps[0] for x in gaps],
        "independent_action_gaps_over_clock": target,
        "maximum_absolute_record_action_residual": residual,
        "common_clock_lapse_units_or_Higgs_scaling_change_family_ratios": False,
        "phase_window_is_admitted_prior_contract": True,
        "finite_shot_or_physical_preparation_error_bound": False,
        "dynamics_replayed_in_this_audit": False,
    }


def nonlinear_shapes(retained):
    """Certified logarithmic shape of frozen controls, not new textures."""
    from .nonlinear_shapes import report
    return report(retained)


def _hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _bindings(retained):
    from .native_family_dynamics import SOURCE_FILES
    paths = [BASELINE, BASELINE.with_suffix(".npz"), NOTE]
    paths += sorted(HERE.glob("*.py")) + sorted(HERE.glob("*.md"))
    paths += [ROOT/name for name in SOURCE_FILES]
    paths += [ROOT/"model/observer_time_fields/REFINEMENT.md"]
    result = dict(retained["source_bindings"])
    result.update({path.relative_to(ROOT).as_posix(): _hash(path) for path in paths})
    return dict(sorted(result.items()))


def run():
    candidate_run.verify(BASELINE)
    retained = strict_json.load(BASELINE)
    initial = _bindings(retained)
    candidate_run._check_bindings(initial)
    algebra = source_algebra()
    from .native_covariance import report as covariance_report
    from .quadratic_response import report as quadratic_report
    from .native_family_dynamics import report as dynamics_report
    from .charged_sector import report as charged_report
    from .collective_response import report as collective_report
    from .native_exterior_attachment import report as exterior_report
    dynamics = dynamics_report()
    if any(initial.get(name) != value for name, value in dynamics["source_bindings"].items()):
        raise RuntimeError("native response input missing from diagnosis custody")
    result = {
        "schema": "oph.hierarchy_diagnosis.v1",
        "source_algebra": algebra,
        "retained_observer_readout": record_audit(retained),
        "separate_benchmark": benchmark_exclusion(algebra, retained),
        "fixed_nonlinear_controls": nonlinear_shapes(retained),
        "native_connected_covariance": covariance_report(),
        "quadratic_poles_and_state_independence": quadratic_report(),
        "native_family_dynamics": dynamics,
        "native_charged_sector_extension": charged_report(),
        "native_collective_charged_response": collective_report(),
        "native_exterior_matter_attachment": exterior_report(),
        "conclusion": {
            "linear_trial_and_entire_complex_affine_feature_class_excluded": True,
            "clock_or_units_explain_current_hierarchy_failure": False,
            "new_mass_generating_dynamics_implemented": False,
            "all_OPH_architectures_excluded": False,
            "observed_particle_masses_derived": False,
        },
        "source_bindings": initial,
    }
    if _bindings(retained) != initial:
        raise RuntimeError("diagnosis inputs changed during calculation")
    candidate_run._check_bindings(initial)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--output", type=Path)
    modes.add_argument("--verify", "--check", dest="verify", type=Path)
    args = parser.parse_args()
    if args.output is not None and args.output.exists():
        parser.error("choose a fresh output; retained reports are not overwritten")
    current = run()
    if args.verify is not None:
        compare_replay(strict_json.load(args.verify), current)
        print("Verified by recomputing the diagnosis; original source dynamics were not replayed")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("x") as handle:
            json.dump(current, handle, indent=2, allow_nan=False)
            handle.write("\n")
        print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
