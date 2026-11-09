from fractions import Fraction as F
import json
from pathlib import Path

import pytest
import sympy as s

from . import state as m


def test_native_character_normalization_and_full_endpoint_twirled_law():
    control = m.exact_one_link_control()
    assert control['single_link_normalization_squared'] == s.Rational(91, 360)
    rows = control['complete_irrep_records']
    assert len(rows) == 11
    assert sum(row['original_probability'] for row in rows) == 1
    assert sum(row['matrix_coefficient_states'] for row in rows) == 119
    assert sum(row['matrix_coefficient_states']*row['density_eigenvalue'] for row in rows) == 1
    assert sum(row['original_probability']*row['electric_energy'] for row in rows) == s.Rational(33, 9100)
    assert sum(row['original_probability'] for row in rows[1:]) == s.Rational(1, 91)
    assert not control['source_only_Gauss_singlet']
    assert not control['periodic_joint_twirl_replaced_by_independent_link_twirl']


def test_one_link_haar_norm_lower_bound_at_actual_L1():
    link = m.link_bound(1, 1)
    lower = link['normalization_squared_lower']
    exact_lower = lower['rational_factor']/lower['integer_base']**lower['inverse_power']
    assert 0 < exact_lower < F(91, 360)
    assert lower['inverse_power'] == 450
    assert not link['normalization_contains_total_volume']
    assert link['bad_link_operator_threshold'] == F(1, 8)


def test_dyadic_exponent_is_integer_and_has_correct_quadratic_threshold():
    for n in (1, 2, 16, 2**20):
        for r in (1, 3, 7):
            link = m.link_bound(n, r)
            cover = (8*n).bit_length()
            assert 2**(cover-1) < 1+8*n <= 2**cover
            assert link['uncapped_binary_exponent'] == 2+450*cover-(n*n)//(480*4**r)
            assert link['electric_energy_per_link_upper'] == F(29, 10)*n**4


@pytest.mark.parametrize('bits,degree,spacing', [(4, 0, '1'), (20, 128, '1'), (28, 4, '3/2')])
def test_uniform_density_planner_certifies_all_allocated_allowances(bits, degree, spacing):
    p = m.density_plan(bits, degree, spacing)
    eta = F(1, 2**bits)
    delta = p['single_link']['plaquette_operator_threshold']
    q = p['chosen_bad_link_probability_upper']
    assert 3*delta*delta/2 <= eta/(4*48000)
    assert 24*q <= eta/(4*48000)
    assert p['single_link']['uncapped_binary_exponent'] <= -(q.denominator.bit_length()-1)
    assert p['initial_magnetic_cost_per_cell_upper'] == 3*delta*delta/2+24*q
    assert p['initial_electric_energy_per_cell_upper'] == F(87, 10)*p['single_link']['per_link_word_degree']**2
    e = p['initial_electric_energy_per_cell_upper']
    b = p['full_Fock_kinetic_norm_per_cell_upper']
    kappa = p['declared_plaquette_coefficient']
    assert 48000*(e+2*b)/kappa == eta/2
    assert 0 < p['all_time_rough_plaquette_probability_upper'] <= eta


def test_every_volume_has_identical_local_state_degree_action_and_probability_bound():
    plan = m.density_plan()
    instances = [m.volume_instance(plan, side, extra_degree=8) for side in (1, 2, 3, 4, 16)]
    assert len({row['same_kappa'] for row in instances}) == 1
    assert len({row['per_link_initial_word_degree'] for row in instances}) == 1
    assert len({row['same_all_time_local_probability_bound'] for row in instances}) == 1
    for row in instances:
        assert row['initial_total_word_degree_upper'] == 3*row['cells']*row['per_link_initial_word_degree']
        assert row['evolution_total_word_cutoff'] == row['initial_total_word_degree_upper']+8
        assert row['full_Fock_dimension']['exponent'] == 60*row['cells']
        assert row['all_Fock_occupations_retained']
        assert row['complete_detector_probability_sum'] == 1
        assert not row['source_only_Gauss_projection_performed']


def test_shared_gauge_twirl_is_not_independent_link_twirling():
    # Exact finite group control for the distinction in the general proof.
    # The actual source group/representation is checked in the native table.
    x = s.Matrix([[0, 1], [1, 0]])
    z = s.diag(1, -1)
    g = s.kronecker_product(x, x)
    seed = s.diag(1, 0, 0, 0)
    twirled = (seed+g*seed*g.H)/2
    invariant = s.kronecker_product(z, z)
    assert g*invariant*g.H == invariant
    assert s.trace(twirled*invariant) == s.trace(seed*invariant) == 1
    assert s.trace((s.eye(4)/4)*invariant) == 0
    # The invariant mixed state contains both irreducible charges.
    assert g*twirled == twirled*g
    assert s.trace(twirled*(s.eye(4)-g)/2) == s.Rational(1, 2)
    psi = s.zeros(16, 1)
    psi[0], psi[15] = 1/s.sqrt(2), 1/s.sqrt(2)
    assert s.kronecker_product(g, g.conjugate())*psi == psi
    assert psi.H*psi == s.ones(1, 1)


def test_translation_symmetry_does_not_supply_orientation_symmetry():
    orientation_means = [F(1, 100), F(2, 100), F(9, 100)]
    per_cell_sum = sum(orientation_means)
    assert all(value <= per_cell_sum for value in orientation_means)
    assert orientation_means[-1] > per_cell_sum/3
    assert m.density_plan()['no_orientation_factor_three_division']


def test_input_domains_are_exact_and_positive():
    for args in ((0, 1), (1, 0), (1.0, 1), (True, 1)):
        with pytest.raises(ValueError):
            m.link_bound(*args)
    with pytest.raises(ValueError):
        m.density_plan(0)
    with pytest.raises(ValueError):
        m.density_plan(degree=-1)
    with pytest.raises(ValueError):
        m.density_plan(spacing=0.1)
    with pytest.raises(ValueError):
        m.density_plan(spacing='0')
    with pytest.raises(ValueError):
        m.volume_instance(m.density_plan(), 0)
    with pytest.raises(ValueError):
        m.volume_instance(m.density_plan(), 2, -1)


def test_complete_certificate_preserves_source_and_reference_boundaries():
    certificate = m.certificate()
    assert certificate['schema'] == 'oph.local_gauge_state.v1'
    root = Path(__file__).resolve().parents[2]
    assert all((root/path).is_file() for path in certificate['source_inputs'])
    plan = certificate['uniform_density_plan']
    for key in ('field_only_Gauss_singlet_support_claimed',
                'global_Wilson_gap_probability_bound_claimed',
                'local_Wilson_low_mode_bound_claimed',
                'physical_mass_or_species_claimed',
                'interacting_relativistic_refinement_claimed'):
        assert not plan[key]
    assert not certificate['comparison_to_occupied_origin_preparation_claimed']
    assert not certificate['numerical_exponential_used_as_proof']
    assert plan['measured_particle_inputs'] == []
    assert len(json.dumps(certificate)) < 40000
