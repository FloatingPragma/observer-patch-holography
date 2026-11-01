"""Independent finite CAR, branch, winding and locality controls."""
from fractions import Fraction as F
import json
from pathlib import Path

import numpy as np
import pytest
import sympy as s

from ..full_fock_field.field import FullFock
from . import momentum as m
from . import packet


@pytest.fixture(scope='module')
def ring():
    return packet.ring_data()


def test_complete_ring_spectral_and_source_certificate(ring):
    report = packet.ring_certificate(ring)
    assert report['exact_spectral_ranks'] == [64]*4
    assert report['contact_pair_matrix_element'] == '1/4'
    assert report['additive_log_commutator_matrix_element_in_pi_units'] == '1/2'
    assert report['ideal_adjacent_FSWAP_count'] == 12


def test_three_axes_cover_the_full_symbolic_algebra():
    result = packet.torus_certificate()
    assert result['CAR_modes'] == 54
    assert result['full_Fock_dimension'] == 2**54
    assert result['exact_translation_covariances'] == 9


@pytest.mark.parametrize('sides,axis', [((3, 1, 1), 0), ((1, 3, 1), 1), ((1, 1, 3), 2), ((1, 1, 4), 2)])
def test_source_action_against_independent_jordan_wigner(sides, axis):
    geometry = m.PeriodicGeometry(sides)
    space = FullFock(geometry.modes)
    a = [space.annihilator(j).toarray() for j in range(geometry.modes)]
    single = np.asarray(m.one_particle(geometry), dtype=complex)
    number = [x.conj().T@x for x in a]
    independent = sum(single[p, q]*(a[p].conj().T@a[q])
                      for p in range(geometry.modes) for q in range(geometry.modes))
    independent += sum(number[2*x]@number[2*x+1] for x in range(len(geometry.sites)))
    h, _, _ = m.action(geometry)
    assert np.array_equal(np.asarray(m.clifford_matrix(h, geometry.modes), dtype=complex), independent)
    shift = np.asarray(m.graded_matrix(geometry.shift(axis)), dtype=complex)
    for j in range(geometry.modes):
        assert np.array_equal(shift@a[j]@shift.conj().T, a[geometry.shift(axis)[j]])


def test_prime_ring_spectral_log_is_conserved():
    geometry = m.PeriodicGeometry((1, 1, 3))
    shift = m.graded_matrix(geometry.shift(2))
    p, projectors = m.spectral_log(shift, 3)
    h = m.clifford_matrix(m.action(geometry)[0], 6)
    assert m.is_zero(h*p-p*h)
    assert m.is_zero(p-p.conjugate().T)
    assert [s.trace(q) for q in projectors] == [24, 20, 20]
    assert m.is_zero(p**3-s.Rational(4, 9)*p)


def test_nonquadratic_even_allowing_a_scalar_identity(ring):
    _, _, p, _, _, _, q, additive, _ = ring
    assert p[0, 0] == 0  # Any scalar in a proposed scalar+dGamma(q) must vanish.
    one = [1 << j for j in range(8)]
    assert p.extract(one, one) == q
    assert not m.is_zero(p-additive)
    # Fully occupied state is invariant under the graded site shift, while
    # the one-particle trace fixes the additive momentum to 2*pi.
    assert p[255, 255] == 0
    assert additive[255, 255] == 2


def test_winding_is_not_itself_conserved(ring):
    _, _, p, _, h, _, _, additive, winding = ring
    assert m.is_zero(h*p-p*h)
    assert not m.is_zero(h*winding-winding*h)
    assert m.is_zero(h*additive-additive*h-2*(h*winding-winding*h))
    assert m.is_zero(winding**3-winding)


def test_both_nyquist_conventions_are_explicit(ring):
    _, _, p, projectors, h, _, _, _, _ = ring
    minus = p-2*projectors[2]
    assert m.is_zero(p*projectors[2]-projectors[2])
    assert m.is_zero(minus*projectors[2]+projectors[2])
    assert m.is_zero(h*minus-minus*h)
    assert not m.is_zero(minus-p)


@pytest.mark.parametrize('permutation', [(2, 3, 0, 1), (3, 1, 0, 2), (5, 0, 2, 4, 3, 1)])
def test_adjacent_fermionic_words_keep_all_occupation_signs(permutation):
    assert m.word_columns(len(permutation), m.adjacent_word(permutation)) == m.graded_columns(permutation)


def test_order_is_not_fermion_parity_on_the_wrapping_bond():
    geometry = m.PeriodicGeometry((1, 1, 3))
    shift = m.graded_matrix(geometry.shift(2))
    assert shift**3 == s.eye(64)
    wrong = m.graded_matrix(geometry.shift(2), grading=False)
    assert shift != wrong
    a = m.clifford_matrix(m.w.annihilator(0), 6)
    moved = m.clifford_matrix(m.w.annihilator(2), 6)
    assert m.is_zero(shift*a*shift.conjugate().T-moved)
    assert not m.is_zero(wrong*a*wrong.conjugate().T-moved)


@pytest.mark.parametrize('length', [3, 4, 5, 6, 7, 8, 11, 16, 31])
def test_kernel_formula_by_exact_geometric_polynomial(length):
    z = s.Symbol('z')
    labels = m.branch_labels(length)
    a = min(labels)
    moment = sum(j*z**(j % length) for j in labels)
    expression = s.rem(s.expand((1-z)*moment+length*z**(a % length)), z**length-1, z)
    assert s.rem(expression, sum(z**j for j in range(length)), z) == 0


def test_even_kernel_matches_the_actual_one_particle_matrix(ring):
    q = ring[6]
    for x in range(4):
        for y in range(4):
            assert s.simplify(q[2*x, 2*y]-m.principal_kernel(4, x-y)) == 0
            assert q[2*x, 2*y+1] == 0
    # Diameter-one operators cannot connect opposite sites: error >= pi/4.
    assert q[0, 4] == F(1, 4)


@pytest.mark.parametrize('length', [3, 5, 7, 17, 97])
def test_all_branch_prime_obstruction_has_nonzero_integer_remainder(length):
    result = m.prime_branch_certificate(length)
    assert result['constant_remainder_mod_period'] == 2
    assert len(set(result['spin_trace_residues'])) == length
    assert result['cyclotomic_coefficients'] == [1]*length


def test_local_higher_body_terms_cannot_fill_a_remote_one_particle_entry():
    # Check every Majorana monomial supported in any single site or nearest
    # pair on a four-site ring. This includes odd and number-changing terms.
    for sites in ((0,), (1,), (2,), (3,), (0, 1), (1, 2), (2, 3), (3, 0)):
        gammas = [4*x+j for x in sites for j in range(4)]
        for subset in range(1 << len(gammas)):
            mask = sum(1 << gamma for i, gamma in enumerate(gammas) if subset >> i & 1)
            # The occupation matrix element can only be nonzero if the XOR
            # flip pattern takes the particle at site two to site zero.
            flipped = 0
            for gamma in gammas:
                if mask >> gamma & 1:
                    flipped ^= 1 << (gamma//2)
            assert (1 << 4)^flipped != 1 << 0


def test_controlled_smooth_initial_sector_is_not_an_invariant_truncation(ring):
    initial, outside = m.fourier_pair(4, 0, 0), m.fourier_pair(4, 1, -1)
    assert s.expand((outside.conjugate().T*ring[5]*initial)[0]) == F(1, 4)
    budget = m.smooth_sector_budget(2, F(1, 4), F(1, 100))
    assert budget['initial_sine_error'] == F(1, 192)
    assert budget['interacting_sector_leakage_upper'] == F(1, 100)
    assert budget['evolved_sine_error_upper'] < F(57, 1000)
    assert m.smooth_sector_budget(1, F(1, 4), 100)['interacting_sector_leakage_upper'] == 0
    assert m.smooth_sector_budget(2, F(1, 4), 100, 0)['interacting_sector_leakage_upper'] == 0
    for args in [(2, F(3, 2), 0), (0, F(1, 4), 0), (2, 0.25, 0), (2, F(1, 4), -1)]:
        with pytest.raises(ValueError):
            m.smooth_sector_budget(*args)


@pytest.mark.parametrize('bad', [(2, 1, 1), (True, 1, 1), (5, 5, 5), (0, 1, 1)])
def test_geometry_rejects_duplicate_or_unbounded_bonds(bad):
    with pytest.raises(ValueError):
        m.PeriodicGeometry(bad)


@pytest.mark.parametrize('bad', [True, -1, 0.125, float('nan')])
def test_exact_coefficients_are_not_silently_rounded(bad):
    with pytest.raises(ValueError):
        m.PeriodicGeometry(speed=bad)


@pytest.mark.parametrize('bad', [2, 4, 9, 99, True])
def test_prime_theorem_rejects_outside_its_domain(bad):
    with pytest.raises(ValueError):
        m.prime_branch_certificate(bad)


def test_bad_period_and_permutations_fail_closed():
    with pytest.raises(ValueError):
        m.spectral_log(s.eye(2), 5)
    with pytest.raises(ValueError):
        m.spectral_log(s.diag(1, 2), 4)
    for permutation in ((0, 0), (1, 2), (True, 0), ()):
        with pytest.raises(ValueError):
            m.graded_columns(permutation)


def test_supplied_receipt_path_is_actually_checked(monkeypatch, tmp_path):
    receipt = tmp_path/'supplied.json'
    receipt.write_text('{"all_pass": false}')
    monkeypatch.setattr(packet, 'run', lambda: {'all_pass': True})
    monkeypatch.setattr('sys.argv', ['periodic_momentum', '--check', str(receipt)])
    with pytest.raises((AssertionError, ValueError, SystemExit)):
        packet.main()
