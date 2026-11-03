#!/usr/bin/env python3
"""Regression tests for GitHub issue #539's MaxEnt multiplier-family acceptance test."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import maxent_closure_acceptance as acceptance  # noqa: E402

from maxent_closure_acceptance import (  # noqa: E402
    NEGATIVE_DEFECT_TOL,
    cell_densities,
    closure_defect,
    constrained_hamiltonian,
    decimate,
    duhamel_covariance,
    gibbs_state,
    global_sum_constraints,
    i_projection,
    independent_operator_count,
    projection_diagnostics,
    relative_entropy,
    run_acceptance,
    run_lattice_pair,
    trace_norm,
)


def test_multiplier_count_is_cutoff_independent_for_global_sums() -> None:
    for n_sites in (3, 4, 6, 8):
        constraints = global_sum_constraints(n_sites)
        assert independent_operator_count(constraints) == 2
        assert sum(len(ops) for ops in cell_densities(n_sites)) == 2 * n_sites


def test_i_projection_moment_matches_and_is_unique() -> None:
    fine = global_sum_constraints(6)
    coarse = global_sum_constraints(3)
    omega, _ = gibbs_state(fine, np.array([0.7, 0.4]))
    sigma = decimate(omega, 6)
    lam_star, residual = i_projection(sigma, coarse)
    assert residual < 1e-9
    hessian_floor = np.linalg.eigvalsh(duhamel_covariance(coarse, lam_star)).min()
    assert hessian_floor > 1e-9


def test_pinsker_residual_bound() -> None:
    run = run_lattice_pair(6, np.array([0.7, 0.4]))
    assert run["trace_norm_residual"] <= run["pinsker_residual_bound"] + 1e-9
    assert run["closure_defect_nats"] > 1e-6


def test_transverse_field_subfamily_is_exactly_closed() -> None:
    run = run_lattice_pair(6, np.array([0.0, 0.4]))
    assert run["closure_defect_nats"] < 1e-10
    assert abs(run["induced_map_R_multipliers"][0]) < 1e-8
    assert abs(run["induced_map_R_multipliers"][1] - 0.4) < 1e-8


def test_relative_entropy_is_zero_only_on_the_diagonal() -> None:
    coarse = global_sum_constraints(3)
    rho, _ = gibbs_state(coarse, np.array([0.3, 0.5]))
    tau, _ = gibbs_state(coarse, np.array([0.4, 0.5]))
    assert relative_entropy(rho, rho) < 1e-12
    gap = relative_entropy(rho, tau)
    assert gap > 0
    assert trace_norm(rho - tau) <= (2 * gap) ** 0.5 + 1e-9


def test_full_acceptance_run_passes() -> None:
    results = run_acceptance()
    assert results["all_checks_pass"], results["checks"]
    for run in results["generic_branch"] + [results["closed_subfamily"]]:
        assert run["normalized_moment_matching_residual"] < 1e-11
        assert 0 <= run["projection_optimality_gap_bound_nats"] < 1e-9
        assert 0 <= run["projection_trace_distance_bound"] < 1e-4


@pytest.mark.parametrize("lam", ([.5], [.5,.2,.1], [np.nan,.2], [True,False]))
def test_gibbs_cannot_silently_drop_or_corrupt_constraints(lam):
    constraints=global_sum_constraints(3)
    for operation in (gibbs_state,duhamel_covariance):
        with pytest.raises(ValueError):
            operation(constraints,np.asarray(lam))


def test_projection_rejects_nonunique_parameters_and_invalid_target():
    identity=np.eye(2)
    with pytest.raises(ValueError,match="independent"):
        i_projection(identity/2,[identity,identity])
    with pytest.raises(ValueError):
        i_projection(np.diag([1.1,-.1]),[np.diag([1.,-1.])])


def test_projection_budget_does_not_return_an_unconverged_minimizer():
    constraints=global_sum_constraints(3)
    sigma,_=gibbs_state(constraints,np.array([.7,.4]))
    with pytest.raises(RuntimeError,match="converge"):
        i_projection(sigma,constraints,max_iter=1)
    with pytest.raises(ValueError):
        i_projection(sigma,constraints,max_iter=0)


def test_gibbs_underflow_is_not_silently_regularized():
    for operation in (gibbs_state,duhamel_covariance):
        with pytest.raises(ValueError,match="underflow"):
            operation([np.diag([0.,1.])],np.array([1000.]))


def test_legacy_relative_entropy_keyword_order_is_preserved():
    state,reference=np.diag([.8,.2]),np.diag([.3,.7])
    expected=.8*np.log(.8/.3)+.2*np.log(.2/.7)
    assert relative_entropy(sigma=state,rho=reference) == pytest.approx(expected)
    assert abs(relative_entropy(sigma=reference,rho=state)-expected) > .01


def test_finite_inputs_cannot_overflow_into_a_nan_state():
    with pytest.raises(ValueError,match="numerical range"):
        gibbs_state([np.diag([1e308,0.])],np.array([1e308]))


def test_gibbs_normalization_cannot_erase_a_subnormal_positive_weight():
    with pytest.raises(ValueError,match="underflow"):
        gibbs_state([np.diag([0.,0.,0.,745.])],np.array([1.]))


def test_gibbs_reconstruction_cannot_erase_a_faithful_direction():
    # Neither exp(-60) nor its normalization underflows. Adding it to the
    # dominant rotated eigenspace loses it when the dense state is formed.
    x=np.array([[0.,1.],[1.,0.]])
    with pytest.raises(ValueError,match="faithful"):
        gibbs_state([x],np.array([30.]))
    # The response API no longer materializes that unresolvable state.
    assert duhamel_covariance([x], [30.])[0, 0] == pytest.approx(
        1/np.cosh(30.)**2, rel=2e-14, abs=0)
    rho,_=gibbs_state([np.diag([0.,60.])],np.array([1.]))
    assert rho[1,1] == pytest.approx(np.exp(-60),rel=1e-14,abs=0)


@pytest.mark.parametrize("bad", (1j,1+0j,np.bool_(True),True,np.nan,np.inf,[1.]))
def test_projection_rejects_nonreal_or_nonscalar_tolerance(bad):
    with pytest.raises(ValueError,match="finite real scalar"):
        i_projection(np.eye(2)/2,[np.diag([1.,-1.])],tol=bad)


@pytest.mark.parametrize("scale", [1e-100, 1e-5, 1., 1e5, 1e100, -1e-5])
def test_projection_rescaling_preserves_the_classical_closed_form(scale):
    # p/(1-p) = exp(-2*scale*lambda); this control uses no Gibbs producer.
    z = np.diag([1., -1.])
    target = np.diag([.7, .3])
    constraints = [scale*z]
    assert independent_operator_count(constraints) == 1
    lam, residual = i_projection(target, constraints)
    assert scale*lam[0] == pytest.approx(-.5*np.log(7/3), abs=1e-11)
    actual, _ = gibbs_state(constraints, lam)
    np.testing.assert_allclose(actual, target, rtol=0, atol=1e-11)
    report = projection_diagnostics(target, constraints, lam)
    assert report.normalized_residual < 1e-11
    assert residual == report.raw_moment_residual


def test_small_observable_units_cannot_manufacture_convergence():
    # The former raw gradient is 6e-12: below the default tolerance even
    # though the returned maximally mixed state is wrong by 2e-7 in trace norm.
    target = np.diag([.5000001, .4999999])
    constraints = [3e-5*np.diag([1., -1.])]
    lam, residual = i_projection(target, constraints)
    actual, _ = gibbs_state(constraints, lam)
    assert trace_norm(actual-target) < 1e-11
    assert residual < 1e-11
    assert abs(lam[0]) > .006


def test_identity_energy_does_not_erase_noncommuting_gibbs_terms_or_covariance():
    x = np.array([[0., 1.], [1., 0.]])
    z = np.diag([1., -1.])
    h = np.array([.7, -.2])
    radius = np.linalg.norm(h)
    expected = (np.eye(2)-np.tanh(radius)/radius*(h[0]*x+h[1]*z))/2
    # These observables retain X and Z exactly. Summing their full Hamiltonian
    # before removing the identity would erase the Z diagonal contribution.
    constraints = [x+1e16*np.eye(2), z]
    actual, log_z = gibbs_state(constraints, h)
    np.testing.assert_allclose(actual, expected, rtol=0, atol=2e-15)
    assert log_z == pytest.approx(np.log(2*np.cosh(radius))-h[0]*1e16)
    # Hessian of log(2 cosh ||h||), derived directly in Pauli coordinates.
    tangent = np.tanh(radius)/radius
    expected_cov = tangent*np.eye(2)+(1/np.cosh(radius)**2-tangent)*np.outer(h, h)/radius**2
    np.testing.assert_allclose(duhamel_covariance(constraints, h), expected_cov,
                               rtol=1e-13, atol=1e-14)


@pytest.mark.parametrize("mixing,offsets", [
    (np.diag([1e-8, 1e8]), [0., 0.]),
    (np.array([[1., 0.], [1., 1e-4]]), [0., 0.]),
    (np.array([[2., -1.], [1., 3.]]), [1e8, -1e8]),
])
def test_noncommuting_projection_uses_the_same_family_after_affine_recharting(mixing, offsets):
    x = np.array([[0., 1.], [1., 0.]])
    y = np.array([[0., -1j], [1j, 0.]])
    z = np.diag([1., -1.])
    # The Y moment is unconstrained: copying the target would be incorrect.
    target = (np.eye(2)+.4*x+.2*y+.3*z)/2
    expected = (np.eye(2)+.4*x+.3*z)/2
    constraints = [row[0]*z+row[1]*x+offset*np.eye(2)
                   for row, offset in zip(mixing, offsets)]
    assert independent_operator_count(constraints) == 2
    lam, residual = i_projection(target, constraints)
    actual, _ = gibbs_state(constraints, lam)
    np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-11)
    # For a Bloch vector r, h = -atanh(|r|) r/|r| in exp(-h.sigma).
    expected_h = -np.arctanh(.5)*np.array([.3, .4])/.5
    np.testing.assert_allclose(mixing.T @ lam, expected_h, rtol=0, atol=2e-11)
    report = projection_diagnostics(target, constraints, lam)
    assert report.normalized_residual < 1e-11
    assert residual == report.raw_moment_residual


def test_imaginary_hermitian_direction_is_not_lost_in_real_multiplier_coordinates():
    y = np.array([[0., -1j], [1j, 0.]])
    target = (np.eye(2)+.6*y)/2
    constraints = [1e-8*y+1e16*np.eye(2)]
    lam, _ = i_projection(target, constraints)
    assert lam[0]*1e-8 == pytest.approx(-np.arctanh(.6), abs=1e-11)
    actual, _ = gibbs_state(constraints, lam)
    np.testing.assert_allclose(actual, target, rtol=0, atol=1e-11)


def test_normalization_does_not_promote_dependent_or_unresolved_constraints():
    x = np.array([[0., 1.], [1., 0.]])
    z = np.diag([1., -1.])
    for constraints in ([x, 7*x+3*np.eye(2)], [x, x+1e-16*z]):
        assert independent_operator_count(constraints) == 1
        with pytest.raises(ValueError, match="independent"):
            i_projection(np.eye(2)/2, constraints)
    assert independent_operator_count([np.eye(2)]) == 0


def test_original_multipliers_must_preserve_convergence_after_ill_conditioned_conversion():
    x = np.array([[0., 1.], [1., 0.]])
    z = np.diag([1., -1.])
    target = (np.eye(2)+.4*x+.3*z)/2
    constraints = [z, z+1e-10*x]
    # The normalized span is independent, but translating a solution back
    # requires cancellation of large multipliers. Never return a false pass.
    try:
        lam, residual = i_projection(target, constraints)
    except RuntimeError as error:
        assert "replay" in str(error)
    else:
        actual, _ = gibbs_state(constraints, lam)
        assert residual < 1e-11
        np.testing.assert_allclose(actual, target, rtol=0, atol=1e-11)



def test_negative_roundoff_defect_keeps_its_magnitude():
    # Trace-one states that differ only at roundoff can still give a small
    # negative entropy estimate. The old max(..., 0.0) reported exact closure
    # for exactly the case that signals a failed comparison.
    sigma = np.diag([1.0 - 1e-12, 0.0])
    reference = np.diag([1.0 - 2e-13, 2e-13])
    raw = relative_entropy(sigma, reference)
    assert -NEGATIVE_DEFECT_TOL < raw < 0
    assert closure_defect(sigma, reference) == pytest.approx(-raw, rel=1e-6, abs=0)
    assert closure_defect(sigma, reference) > 0


def test_significantly_negative_defect_is_rejected(monkeypatch):
    # A state that far outside the declared trace tolerance never reaches this
    # guard: the state checks reject it first. Inject the value the guard has
    # to refuse rather than read as zero.
    monkeypatch.setattr(acceptance, "relative_entropy", lambda sigma, rho: -1e-6)
    with pytest.raises(ValueError, match="negative relative entropy"):
        acceptance.closure_defect(np.eye(2) / 2, np.eye(2) / 2)


@pytest.mark.parametrize("energies", [[0., 2.], [0., 2., 4.]])
def test_partition_retains_small_energy_after_large_identity_offsets_cancel(energies):
    energies = np.array(energies)
    constraints = [np.diag(1e16+energies), -1e16*np.eye(len(energies))]
    # Their sum is exactly diag(energies). Form the expected probabilities and
    # partition directly, without the centered Gibbs producer.
    actual, log_z = gibbs_state(constraints, np.ones(2))
    expected_z = np.exp(-energies).sum()
    np.testing.assert_array_equal(constrained_hamiltonian(constraints, np.ones(2)),
                                  np.diag(energies))
    np.testing.assert_allclose(actual, np.diag(np.exp(-energies))/expected_z,
                               rtol=0, atol=2e-15)
    assert log_z == pytest.approx(np.log(expected_z), abs=2e-15)


@pytest.mark.parametrize("scale", [1e-154, 1e154])
def test_representable_covariance_survives_large_or_small_observable_units(scale):
    # At the uniform state, Var(scale Z) = scale**2. Symmetrization must not
    # overflow an otherwise finite result by adding it to itself first.
    covariance = duhamel_covariance([scale*np.diag([1., -1.])], np.zeros(1))
    assert covariance[0, 0] == pytest.approx(scale**2, rel=1e-14, abs=0)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
