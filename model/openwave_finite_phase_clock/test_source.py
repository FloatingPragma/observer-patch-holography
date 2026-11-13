from fractions import Fraction as F
import numpy as np
import pytest
import sympy as s
from scipy.linalg import expm
from .source import (exchange, prepared, free_source_law,
                     complete_terminal_records)


def test_native_transfer_dense_exponential_and_full_complements():
    for c, u in ((F(3, 5), F(4, 5)), (F(5, 13), F(12, 13))):
        N, A, R, B, H, U = exchange(c, u)
        independently_exponentiated = expm(-1j*np.arctan2(float(u), float(c))*np.asarray(B, complex))
        assert np.max(abs(independently_exponentiated-np.asarray(U, complex))) < 4e-16
        for i in (1, 2):
            for j in range(2):
                v = s.eye(8)[:, 2*i+j]
                assert U*v == v
        _, rho, record = prepared(c, u)
        assert s.trace(N*rho) == s.Rational(u)**2/2
        assert record['source_plus_reference_energy_balance_exact'] == '0'
    with pytest.raises(ValueError):
        exchange(F(1, 2), F(1, 2))


def test_source_frequency_changes_with_mass_not_preparation_strength():
    # Independently evolve at common source times, without dividing time by m.
    for mass in (1., 2.7, 206.77):
        for c, u in ((F(3, 5), F(4, 5)), (F(5, 13), F(12, 13))):
            N, A, *_ = exchange(c, u)
            _, rho, _ = prepared(c, u)
            for t in (.003, .007, .011):
                U = expm(-1j*mass*t*np.asarray(N, complex))
                value = np.trace(np.asarray(A, complex)@U@np.asarray(rho, complex)@U.conj().T)
                assert abs(value-1j*float(u)*np.exp(1j*mass*t)/2) < 2e-15
                assert abs(value-1j*float(u)*np.exp(2j*mass*t)/2) > 1e-4
    assert free_source_law(twirl=True)['free_source_A_expectation_exact'] == '0'


def test_complete_detector_cancels_coevolving_phase_and_exposes_fresh_resource():
    normal = [complete_terminal_records(k) for k in range(4)]
    fresh = [complete_terminal_records(k, fresh_terminal_phase=True) for k in range(4)]
    assert all(x['complete_sixteen_records'] == normal[0]['complete_sixteen_records'] for x in normal)
    assert [F(x['readout_excited_probability_exact']) for x in fresh] == [F(241,1250), F(481,1250), F(721,1250), F(481,1250)]
    for k in range(4):
        row = complete_terminal_records(k, fresh_terminal_phase=True, twirl=True)
        assert F(row['readout_excited_probability_exact']) == F(481, 1250)
        assert sum(F(r['probability_exact']) for r in row['complete_sixteen_records']) == 1


def test_energy_incompatible_effect_is_a_genuine_negative_control():
    H = np.diag([0., 1., 1., 2.])
    v = np.array([1, 1j, 1, -1j], complex)/2
    rho = np.outer(v, v.conj())
    rng = np.random.default_rng(1702)
    # Energy-compatible unitary rotates the full degenerate energy-one block.
    G = np.zeros((4,4), complex)
    M = rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))
    G[1:3,1:3] = M+M.conj().T
    V = expm(-1j*G)
    E = np.diag([1., 0., 1., 0.])
    bad = np.ones((4,4))/4
    compatible, changed = [], []
    for t in (0., .7, 1.9):
        U = expm(-1j*H*t)
        final = V@U@rho@U.conj().T@V.conj().T
        compatible.append(np.trace(E@final).real)
        changed.append(np.trace(bad@final).real)
    assert np.ptp(compatible) < 3e-16
    assert np.ptp(changed) > .01


def test_static_controller_local_identity_on_every_column():
    from .controller import source, SCALE
    from scipy.sparse import eye, kron, diags, csr_matrix
    H, _, _, A, B1, B2, prep, read = source(2.7, 3.)
    native = kron(eye(len(prep)), csr_matrix(A), format='csr')
    exact_rhs = (2.7/SCALE*native+
        kron(diags(prep*3.), csr_matrix(B1@A-A@B1), format='csr')+
        kron(diags(read*2.), csr_matrix(B2@A-A@B2), format='csr'))
    residual = H@native-native@H-exact_rhs
    assert residual.shape == (4096,4096)
    assert not residual.nnz or max(abs(residual.data)) < 3e-15
    # Omitting the reference/controller terms fails on their actual columns.
    bad = H@native-native@H-2.7/SCALE*native
    assert max(abs(bad.data)) == 3.


def test_static_controller_keeps_all_clock_positions_and_readout_work():
    from .controller import run
    row = run(1., 2., twirl=True)
    assert np.asarray(row['all_native_reference_controller_probability_rows']).shape == (256,16)
    assert row['complete_record_total_diagnostic'] == pytest.approx(1, abs=1e-12)
    assert row['minimum_free_source_coherence_diagnostic'] < 1e-12
    assert row['terminal_complete_occupation_measurement_work_diagnostic'] != 0
