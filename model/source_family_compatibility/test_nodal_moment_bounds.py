"""Exact local CAR and independent signed-star spectral controls."""
from fractions import Fraction as F
from functools import reduce

import numpy as np
import pytest
import sympy as s

from model.source_family_compatibility import momentum_response as momentum
from model.source_family_compatibility import nodal_moment_bounds as bounds


def test_source_orbitals_use_the_actual_thirteen_site_star_and_eta_form():
    row = bounds.source_orbitals()
    assert len(set(row['physical_sites'])) == 13
    assert row['canonical_orbital_count'] == 8
    assert row['orbital_integer_norm_squares'] == (1, 6, 1, 1, 1, 1, 1, 1)
    assert row['orbital_eta_signs'] == (1, -1, 1, 1, 1, 1, 1, 1)
    assert row['proper_cubic_lifts_checked'] == 24
    assert not row['neighbour_orbital_is_supplied_extra_mode']
    assert row['all_six_axial_correlations_equal_requires_stated_cubic_average']


def test_complete_local_type_decomposition():
    row = bounds.local_doublet_audit()
    assert row['local_Fock_dimension'] == 4
    assert row['spin_doublet_rank']+row['eta_doublet_rank'] == 4
    assert row['both_staggered_eta_signs_checked']
    assert not row['occupation_sectors_discarded']


@pytest.mark.parametrize('sign', [1, -1])
def test_complete_three_orbital_CAR_blocks_with_signed_inactive_leaves(sign):
    row = bounds.full_CAR_three_orbital_audit(sign)
    assert row['full_Fock_dimension'] == 64
    assert row['complete_local_type_blocks_checked'] == 8
    assert row['inactive_leaf_multiplicities_retained']
    assert row['weak_coupling_sign'] == sign


@pytest.mark.parametrize('sign,radical,offset', [(1, 577, 13), (-1, 673, 11)])
def test_full_star_spectrum_against_independent_tensor_Paulis(sign, radical, offset):
    x = np.array([[0, 1], [1, 0]], complex)
    y = np.array([[0, -1j], [1j, 0]], complex)
    z = np.diag([1, -1]).astype(complex)
    identity = np.eye(2)
    direct = np.zeros((256, 256), complex)
    for leaf, coupling in enumerate((12,)+(sign,)*6, 1):
        for pauli in (x, y, z):
            direct += coupling*reduce(np.kron, [pauli if j in (0, leaf) else identity
                                                for j in range(8)])
    # Tensor conventions reverse the bit numbering but preserve this star.
    literal = np.array(bounds.spin_star((12,)+(sign,)*6), float)
    np.testing.assert_allclose(np.linalg.eigvalsh(direct), np.linalg.eigvalsh(literal), atol=2e-13)
    eigenvalues = np.linalg.eigvalsh(direct)/4
    assert abs(eigenvalues[0]-(-offset-np.sqrt(radical))/4) < 2e-13
    assert np.count_nonzero(np.abs(eigenvalues-eigenvalues[0]) < 1e-12) == 7
    row = bounds.eight_spin_star_audit(sign)
    assert row['magnetization_block_dimensions'] == [1, 8, 28, 56, 70, 56, 28, 8, 1]
    assert row['weak_spin_zero_actual_energies'] == [-9, 3]
    factors = row['exact_characteristic_factors_for_4H']
    assert sum((len(r['coefficients'])-1)*r['multiplicity'] for r in factors) == 256
    assert row['full_eight_orbital_CAR_dimension'] == 65536
    assert not row['all_65536_states_diagonalized_as_one_dense_matrix']


@pytest.mark.parametrize('sign', [1, -1])
def test_inactive_leaf_extension_keeps_the_bound_for_signed_couplings(sign):
    # Every possible number of active identical weak leaves, with either
    # presence or absence of the strong leaf, occurs in a full CAR type block.
    radical, offset = (577, 13) if sign == 1 else (673, 11)
    minimum = (-offset-np.sqrt(radical))/4
    for strong in (0, 12):
        for weak_count in range(7):
            matrix = np.array(bounds.spin_star((strong,)+(sign,)*weak_count), float)/4
            assert np.linalg.eigvalsh(matrix)[0] >= minimum-2e-13


def test_wrong_strong_coupling_is_not_silently_the_source_star():
    expected = bounds.spin_star((12,)+(1,)*6, 4).charpoly().all_coeffs()
    assert bounds.spin_star((11,)+(1,)*6, 4).charpoly().all_coeffs() != expected
    assert bounds.spin_star((12,)+(-1,)*6, 4).charpoly().all_coeffs() != expected
    assert bounds.spin_star(iter((12, 1)), 1) == bounds.spin_star((12, 1), 1)


def test_source_observables_preserve_signed_endpoints_and_complete_scores():
    plus, minus = (bounds.source_observable(sign) for sign in (1, -1))
    qstar = momentum.invariant_observables()['qstar']
    axial = momentum.add(*(momentum.spin_eta_dot(0, j, 1, 1) for j in range(7, 13)))
    assert momentum.add(plus, minus) == momentum.scale(qstar, 24)
    assert momentum.add(plus, momentum.scale(minus, -1)) == momentum.scale(axial, 2)
    for operator in (plus, minus):
        row = momentum.score_plan(operator)
        assert len(row['terms']) == 252
        assert row['score_absolute_bound'] == F(117, 2)
        assert sum(row['setting_probabilities']) == 1
        assert row['all_plus_minus_outcomes_retained']


@pytest.mark.parametrize('value', [0, 1, 577, 673, F(7, 13), F(1, 64)])
def test_outward_rational_roots(value):
    for bits in (1, 7, 40):
        lower, upper = bounds.sqrt_interval(value, bits)
        assert lower*lower <= value <= upper*upper
        assert upper-lower <= F(1, 1 << bits)
    assert bounds.sqrt_interval(F(1, 64), 7) == (F(1, 8), F(1, 8))


@pytest.mark.parametrize('bad', [-1, F(-1, 3), .5, True, '577'])
def test_root_input_rejects_inexact_and_outside_domain(bad):
    with pytest.raises(ValueError):
        bounds.sqrt_interval(bad)


@pytest.mark.parametrize('bad', [0, 129, True, F(7)])
def test_root_precision_contract(bad):
    with pytest.raises(ValueError):
        bounds.sqrt_interval(577, bad)


@pytest.mark.parametrize('bad', [0, 2, True, F(1)])
def test_signed_source_contract(bad):
    for function in (bounds.source_observable, bounds.full_CAR_three_orbital_audit,
                     bounds.eight_spin_star_audit):
        with pytest.raises(ValueError):
            function(bad)


def test_exact_paired_endpoint_ceiling_and_source_homogeneity():
    low577, high577 = bounds.sqrt_interval(577)
    low673, high673 = bounds.sqrt_interval(673)
    assert (29+high673)/4 < (31+low577)/4 < (31+high577)/4 < 14
    row = bounds.atom_budget()
    scaled = bounds.atom_budget(F(1, 8), F(1))
    assert row['original_zero_atom_upper'] == scaled['original_zero_atom_upper']
    assert scaled['all_momentum_Y_Gram_upper']*64 == row['all_momentum_Y_Gram_upper']
    assert row['original_zero_atom_upper'] < F(7, 15)
    assert row['default_ballistic_all_packet_escaping_weight_liminf_lower'] > F(8, 15)
    assert bounds.atom_budget(1, 7)['default_ballistic_all_packet_fixed_energy_limsup_upper'] is None
    assert not row['actual_zero_atom_positive_proved']
    assert not row['source_eighth_moment_computed']


def test_degree_two_polynomial_bounds_atoms_without_claiming_their_existence():
    # Same alpha and first4 moments; changing the hidden clock can remove
    # the original zero atom. These are implication controls, not source fits.
    alpha = s.Integer(4)
    for hidden in (s.Integer(0), s.Integer(1)):
        L = s.Matrix([[0, alpha, 0, 0], [alpha, 0, 1, 0],
                      [0, 1, 0, hidden], [0, 0, hidden, 0]])
        moments = [(L**j)[0, 0] for j in range(5)]
        assert moments == [1, 0, 16, 0, 272]
        polynomial = s.eye(4)-L*L/17
        assert (polynomial*polynomial)[0, 0] == s.Rational(1, 17)
        kernel = L.nullspace()
        if hidden:
            assert not kernel
        else:
            basis = s.Matrix.hstack(*kernel)
            projector = basis*(basis.H*basis).inv()*basis.H
            assert projector[0, 0] == s.Rational(1, 17)


def test_certificate_source_scope_and_fresh_payloads():
    row = bounds.certificate()
    assert row['same_centered_cubic_full_CAR_source']
    assert row['source_global_Y_Gram_over_c_squared_strict_upper'] == 14
    assert row['global_default_fixed_window_limsup_strict_upper'] == F(7, 15)
    assert not row['spin_only_ground_or_antiferromagnetic_order_assumed']
    assert not row['universal_local_bound_attained_in_the_actual_infinite_ground']
    assert not row['native_large_tape_executed']
    assert row['physical_species_identified'] == []
    assert row['physical_mass_prediction'] is None
    row['exact_complete_spin_stars'].clear()
    assert len(bounds.certificate()['exact_complete_spin_stars']) == 2
