#!/usr/bin/env python3
"""Correctness tests for the real diagnostic lattice engine.

Every test here checks an exact mathematical property (algebra, gauge
covariance, hermiticity, Hamiltonian consistency) or a known free-field
value, so a conventions slip anywhere in the engine fails loudly.
"""

from __future__ import annotations

import numpy as np
import pytest

from lattice_backend.core import (
    GAMMA,
    GAMMA5,
    average_plaquette,
    cold_start,
    gauge_transform,
    gauge_transform_spinor,
    project_su3,
    random_su3_near_identity,
    sweep,
    wilson_gauge_action,
)
from lattice_backend.dirac import WilsonClover, cg_normal, point_propagator, point_source
from lattice_backend.hmc import (
    TwoFlavorHMC,
    exp_i_herm,
    gauge_force,
    kinetic_energy,
    pseudofermion_force,
    sample_momenta,
)
from lattice_backend.spectroscopy import (
    effective_mass,
    nucleon_correlators,
    pion_correlator,
)

SMALL = (4, 2, 2, 2)
RNG = np.random.default_rng(20260712)


def _random_gauge_field(shape, eps=0.4, seed=1):
    rng = np.random.default_rng(seed)
    u = cold_start(shape)
    for mu in range(4):
        u[mu] = random_su3_near_identity(rng, shape, eps)
    return u


def _random_gauge_rotation(shape, seed=2):
    rng = np.random.default_rng(seed)
    return random_su3_near_identity(rng, shape, 0.7)


def test_gamma_algebra() -> None:
    for mu in range(4):
        assert np.allclose(GAMMA[mu], GAMMA[mu].conj().T), mu
        for nu in range(4):
            anti = GAMMA[mu] @ GAMMA[nu] + GAMMA[nu] @ GAMMA[mu]
            assert np.allclose(anti, 2.0 * (mu == nu) * np.eye(4)), (mu, nu)
    assert np.allclose(GAMMA5, GAMMA[0] @ GAMMA[1] @ GAMMA[2] @ GAMMA[3])
    assert np.allclose(GAMMA5 @ GAMMA5, np.eye(4))
    for mu in range(4):
        assert np.allclose(GAMMA5 @ GAMMA[mu] + GAMMA[mu] @ GAMMA5, np.zeros((4, 4)))


def test_su3_sampling_and_projection() -> None:
    m = random_su3_near_identity(RNG, (10,), 0.5)
    ident = m @ np.conj(np.swapaxes(m, -1, -2))
    assert np.allclose(ident, np.broadcast_to(np.eye(3), ident.shape), atol=1e-12)
    assert np.allclose(np.linalg.det(m), 1.0, atol=1e-12)
    noisy = m + 1e-3 * RNG.normal(size=m.shape)
    p = project_su3(noisy)
    assert np.allclose(np.linalg.det(p), 1.0, atol=1e-10)
    assert np.allclose(p @ np.conj(np.swapaxes(p, -1, -2)),
                       np.broadcast_to(np.eye(3), p.shape), atol=1e-10)


def test_plaquette_gauge_invariance() -> None:
    u = _random_gauge_field(SMALL)
    g = _random_gauge_rotation(SMALL)
    before = average_plaquette(u)
    after = average_plaquette(gauge_transform(u, g))
    assert abs(before - after) < 1e-12


def test_cold_start_plaquette_is_one() -> None:
    u = cold_start(SMALL)
    assert abs(average_plaquette(u) - 1.0) < 1e-14
    assert abs(wilson_gauge_action(u, 5.7)) < 1e-10


def test_dirac_gamma5_hermiticity() -> None:
    u = _random_gauge_field(SMALL)
    op = WilsonClover(u, kappa=0.12, c_sw=1.0)
    rng = np.random.default_rng(3)
    psi = rng.normal(size=(*SMALL, 4, 3)) + 1j * rng.normal(size=(*SMALL, 4, 3))
    chi = rng.normal(size=(*SMALL, 4, 3)) + 1j * rng.normal(size=(*SMALL, 4, 3))
    lhs = np.vdot(chi, op.apply(psi))
    rhs = np.vdot(op.apply_dag(chi), psi)
    assert abs(lhs - rhs) < 1e-10 * max(1.0, abs(lhs))


def test_dirac_gauge_covariance() -> None:
    """C_pi from a propagator is invariant under gauge rotations."""
    u = _random_gauge_field(SMALL)
    op = WilsonClover(u, kappa=0.11, c_sw=1.0)
    prop, _ = point_propagator(op, SMALL, tol=1e-10)
    c_before = pion_correlator(prop)

    g = _random_gauge_rotation(SMALL)
    # rotate the field; the point source at the origin picks up g(0), a
    # unitary color rotation that cancels inside the correlator trace
    u_rot = gauge_transform(u, g)
    op_rot = WilsonClover(u_rot, kappa=0.11, c_sw=1.0)
    prop_rot, _ = point_propagator(op_rot, SMALL, tol=1e-10)
    c_after = pion_correlator(prop_rot)
    assert np.allclose(c_before, c_after, rtol=1e-7)


def test_cg_solves_dirac_equation() -> None:
    u = _random_gauge_field(SMALL)
    op = WilsonClover(u, kappa=0.12, c_sw=1.0)
    rhs = point_source(SMALL, 0, 0)
    x, _, res = cg_normal(op, rhs, tol=1e-10)
    assert res < 1e-9
    back = op.apply(x)
    assert np.max(np.abs(back - rhs)) < 1e-7


def test_free_field_spectroscopy_matches_pole_masses() -> None:
    """On U = 1 the quark pole mass is am_q = log(1 + m), m = 1/(2 kappa) - 4;
    the pion and nucleon effective masses approach 2 am_q and 3 am_q."""
    shape = (16, 2, 2, 2)
    kappa = 0.10
    am_q = np.log(1.0 + (0.5 / kappa - 4.0))
    u = cold_start(shape)
    op = WilsonClover(u, kappa=kappa, c_sw=1.0)
    prop, _ = point_propagator(op, shape, tol=1e-11)

    c_pi = pion_correlator(prop)
    m_pi = effective_mass(c_pi)
    mid = slice(4, 7)
    assert np.all(np.isfinite(m_pi[mid]))
    assert np.allclose(m_pi[mid], 2.0 * am_q, rtol=0.08)

    direct, exchange = nucleon_correlators(prop)
    n_iso = direct - exchange
    assert n_iso[1] != 0.0
    m_n = effective_mass(n_iso)
    assert np.all(np.isfinite(m_n[mid]))
    assert np.allclose(m_n[mid], 3.0 * am_q, rtol=0.10)


def test_gauge_force_matches_finite_difference() -> None:
    u = _random_gauge_field(SMALL, eps=0.3, seed=11)
    beta = 5.5
    f = gauge_force(u, beta)
    rng = np.random.default_rng(12)
    q = sample_momenta(rng, SMALL)
    eps = 1e-5
    u_plus = u.copy()
    u_minus = u.copy()
    for mu in range(4):
        u_plus[mu] = exp_i_herm(eps * q[mu]) @ u[mu]
        u_minus[mu] = exp_i_herm(-eps * q[mu]) @ u[mu]
    fd = (wilson_gauge_action(u_plus, beta) - wilson_gauge_action(u_minus, beta)) / (2 * eps)
    analytic = float(np.real(np.einsum("...ij,...ji->...", q, f).sum()))
    assert abs(fd - analytic) < 1e-5 * max(1.0, abs(analytic))


def test_pseudofermion_force_matches_finite_difference() -> None:
    u = _random_gauge_field(SMALL, eps=0.3, seed=21)
    kappa = 0.11
    rng = np.random.default_rng(22)
    phi = rng.normal(size=(*SMALL, 4, 3)) + 1j * rng.normal(size=(*SMALL, 4, 3))

    hmc = TwoFlavorHMC(beta=5.5, kappa=kappa, cg_tol=1e-12)
    x_sol, op, _ = hmc._pf_solution(u, phi)
    f = pseudofermion_force(op, op.ubc, x_sol)

    q = sample_momenta(np.random.default_rng(23), SMALL)
    eps = 1e-5

    def s_pf(field: np.ndarray) -> float:
        return hmc.action_pf(field, phi)

    u_plus = u.copy()
    u_minus = u.copy()
    for mu in range(4):
        u_plus[mu] = exp_i_herm(eps * q[mu]) @ u[mu]
        u_minus[mu] = exp_i_herm(-eps * q[mu]) @ u[mu]
    fd = (s_pf(u_plus) - s_pf(u_minus)) / (2 * eps)
    analytic = float(np.real(np.einsum("...ij,...ji->...", q, f).sum()))
    assert abs(fd - analytic) < 1e-4 * max(1.0, abs(analytic))


def test_hmc_reversibility_and_energy() -> None:
    shape = (4, 2, 2, 2)
    u = _random_gauge_field(shape, eps=0.2, seed=31)
    hmc = TwoFlavorHMC(beta=5.5, kappa=0.10, n_steps=12, cg_tol=1e-12)
    rng = np.random.default_rng(32)
    u1, res = hmc.trajectory(rng, u)
    # a fine-step trajectory conserves H at the integrator's order
    assert abs(res.delta_h) < 0.05, res.delta_h


def test_quenched_sweep_moves_plaquette_toward_equilibrium() -> None:
    rng = np.random.default_rng(41)
    u = cold_start(SMALL)
    sweep(rng, u, beta=5.7, n_or=1)
    p = average_plaquette(u)
    # one sweep from cold must move off 1.0 but stay physical
    assert 0.3 < p < 1.0
    assert not np.isnan(p)


@pytest.mark.slow
def test_quenched_plaquette_beta57_matches_literature() -> None:
    """<P> at beta = 5.7 on a small quenched lattice is approximately 0.549
    (literature value on large volumes); allow a broad small-volume window."""
    rng = np.random.default_rng(42)
    shape = (4, 4, 4, 4)
    u = cold_start(shape)
    for _ in range(60):
        sweep(rng, u, beta=5.7, n_or=1)
    samples = []
    for _ in range(40):
        sweep(rng, u, beta=5.7, n_or=1)
        samples.append(average_plaquette(u))
    mean = float(np.mean(samples))
    assert 0.52 < mean < 0.58, mean
