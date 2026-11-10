"""Independent CAR, original-history and full-source scope controls."""
from copy import deepcopy
from fractions import Fraction as F

import numpy as np
import pytest
import sympy as s

from . import model as m


def independent_annihilator(mode, modes):
    dimension = 1 << modes
    matrix = s.zeros(dimension)
    for state in range(dimension):
        if state >> mode & 1:
            matrix[state ^ (1 << mode), state] = (-1)**((state & ((1 << mode)-1)).bit_count())
    return matrix


def test_independent_full_CAR_radial_spectrum_and_sharp_deficit_bounds():
    a = [independent_annihilator(j, 6) for j in range(6)]
    phi = [(a[2*j]*a[5]-a[2*j+1]*a[4])/s.sqrt(2) for j in range(2)]
    S = sum((p.H*p for p in phi), s.zeros(64))
    assert S == m.radial_control()['S']
    holes = s.diag(*(6-k.bit_count() for k in range(64)))
    assert min((holes-2*s.eye(64)+S).eigenvals()) == 0
    assert min((s.Rational(3, 4)*holes-s.eye(64)+S*S/4).eigenvals()) == 0
    # Both constants are optimal: reducing them violates a one-singlet-hole vector.
    v = s.eye(64)[:, 63 ^ (1 << 4)]
    assert (v.H*(s.Rational(99, 100)*holes-2*s.eye(64)+S)*v)[0] < 0
    assert (v.H*(s.Rational(74, 100)*holes-s.eye(64)+S*S/4)*v)[0] < 0


def test_gauge_radial_vs_charged_one_point_and_CAR_parity():
    d = m.charged.dressing_control(); S = m.radial_control()['S']
    Q = -3*d['NL']+6*d['NE']
    assert Q*S == S*Q
    assert any(Q*p != p*Q for p in d['phi'])
    rho = s.eye(64)/64
    assert all(s.trace(rho*p) == 0 for p in d['phi'])
    assert s.trace(rho*S) > 0
    parity = s.diag(*((-1)**k.bit_count() for k in range(64)))
    assert parity*S == S*parity


def test_one_site_contraction_and_complete_multiple_address_failure():
    A = np.array(m.radial_control()['A'], complex)
    value, basis = np.linalg.eigh(A)
    assert value.min() >= -1e-13 and value.max() <= 1+1e-13
    # An independent two-site classical spectral chart retains every radial sector.
    labels = np.array([0, .25, .5, .75, 1])
    left = np.kron(np.diag(labels), np.eye(5))/np.sqrt(2)
    right = np.kron(np.eye(5), np.diag(labels))/np.sqrt(2)
    E = left.conj().T@left+right.conj().T@right
    failure = np.diag(np.sqrt(1-np.diag(E)))
    assert np.linalg.norm(E+failure.conj().T@failure-np.eye(25)) < 1e-13
    assert np.all(np.diag(E) <= 1)


def test_literal_all_native_Wilson_first_actions():
    witness = m.wilson_witness()
    assert witness['all_native_blocks_included'] == 5
    assert witness['actual_two_cell_modes'] == 120
    assert all(x['S_first_commutator_norm_squared'] == F(3, 4) for x in witness['axes'])
    assert all(x['S_double_commutator_norm_squared'] == F(9, 16) for x in witness['axes'])
    assert witness['A_first_norm_squared_coefficient'] == F(9, 16)
    assert witness['A_double_norm_squared_coefficient'] == F(27, 256)


def independent_q6_corner_control(axis):
    H, _, _, _, _ = m.current.q6_two_cell_control(axis)
    states = [(mask, 2-(mask & 15).bit_count()) for mask in range(256)
              if mask.bit_count() == 4 and abs(2-(mask & 15).bit_count()) <= 1]
    negative = np.diag([((mask >> 2) & 3).bit_count() for mask, flux in states])
    omega = np.zeros(len(states))
    omega[states.index((204, 0))] = 1
    return H, negative, omega


@pytest.mark.parametrize('axis', [0, 1, 2])
def test_independent_quantum_link_singlet_hole_matrices(axis):
    H, S, omega = independent_q6_corner_control(axis)
    first = H@S-S@H
    double = S@first-first@S
    assert abs(np.vdot(first@omega, first@omega).real-.5) < 1e-13
    assert abs(np.vdot(double@omega, double@omega).real-.5) < 1e-13
    # With a filled singlet and one spectator weak row, S=1+N_weak/2.
    weak = np.eye(len(omega))+S/2
    first_weak = H@weak-weak@H
    double_weak = weak@first_weak-first_weak@weak
    assert abs(np.vdot(first_weak@omega, first_weak@omega).real-.125) < 1e-13
    assert abs(np.vdot(double_weak@omega, double_weak@omega).real-.03125) < 1e-13


def test_original_spectral_moment_and_Cayley_detector_control():
    H, S, omega = independent_q6_corner_control(1)
    A = S/2
    energies, vectors = np.linalg.eigh(H)
    U = energies.max()+1
    weights = (U-energies)**4; weights /= weights.sum()
    rho = (vectors*weights)@vectors.conj().T
    exact = vectors.conj().T@A@vectors
    transfer = energies[:, None]-energies[None, :]
    mu = abs(exact)**2*weights[None, :]
    C = H@A-A@H
    second = (mu*transfer**2).sum()
    first = (mu*transfer).sum()
    assert first > 0 and second > 0
    assert abs(second-np.trace(rho@C.conj().T@C).real) < 1e-9
    assert abs(first-np.trace(rho@(A@C-C@A)).real/2) < 1e-10
    assert (mu*transfer**2)[transfer < 0].sum() <= (mu*transfer**2)[transfer > 0].sum()+1e-12
    R = max(abs(energies))
    tau = 1/(64*R**2)
    I = np.eye(len(H))
    Q = np.linalg.solve(I+1j*tau*H/2, I-1j*tau*H/2)
    B = Q.conj().T@A@Q
    K = (A@B-B@A)/2
    E = K.conj().T@K
    ev, basis = np.linalg.eigh(E)
    failure = (basis*np.sqrt(np.maximum(0, 1-ev)))@basis.conj().T
    assert np.linalg.norm(E+failure.conj().T@failure-I) < 1e-12
    dilation = np.block([[K, failure], [failure, K]])
    assert np.linalg.norm(dilation.conj().T@dilation-np.eye(2*len(H))) < 2e-12
    assert np.trace(rho@E).real > 0
    assert np.linalg.norm(A@A-A@A) == 0
    # Purified source and its reference give the same marginal probability.
    root = (vectors*np.sqrt(weights))@vectors.conj().T
    assert abs(np.linalg.norm(K@root, 'fro')**2-np.trace(rho@E).real) < 1e-20


def test_large_total_detection_does_not_imply_an_inelastic_pole():
    # Commuting control is deliberately distinct from the interacting primary source.
    A = np.diag([.5, 1]); H = np.diag([0., 2.]); rho = np.diag([.2, .8])
    assert np.trace(rho@A@A) > .76
    assert np.linalg.norm(H@A-A@H) == 0


def test_actual_members_preserve_source_and_capacity_domains():
    for n in (0, 1, 2, 8):
        member = m.contract.member(n); bound = m.bind_member(member); r = bound['result']
        assert all(r[k] == v for k, v in member['ids'].items())
        radial = r['complete_radial_records']; dynamics = r['original_radial_dynamics']
        assert F(radial['original_total_success_lower']) >= F(3481, 160000)
        assert radial['complete_probe_outcomes'] == member['source']['side']**3+1
        if n == 0:
            assert F(radial['original_total_success_lower']) > F(383, 500)
        assert radial['success_includes_elastic_response']
        assert not dynamics['finite_time_success_is_an_energy_transfer_bin']
        assert F(dynamics['complete_clock_duration']) > 0
        assert r['physical_mass_table_GeV'] == r['physical_species_confirmed'] == []
        assert not r['physical_Higgs_field_identified']
    bad = deepcopy(bound)
    bad['result']['physical_mass_table_GeV'] = [{'Higgs': 125}]
    bad['result_sha256'] = m.contract.digest(bad['result'])
    with pytest.raises(ValueError):
        m.verify_bound(member, bad)


def test_the_deleted_hopping_bound_does_not_certify_a_negative_mean():
    result = m.deleted_hopping_comparison()
    assert result['resulting_weak_plus_e_hopping_expectation_upper'] > 0
    assert not result['strict_negative_hopping_from_this_estimate_proved']
    assert not result['actual_scalar_dynamics_excluded']
