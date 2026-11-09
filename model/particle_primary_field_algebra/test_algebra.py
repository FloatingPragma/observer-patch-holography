from copy import deepcopy
from fractions import Fraction as F
import json

import numpy as np
import pytest
import sympy as s

from model.particle_common_contract import contract
from model.particle_joint_reconstruction import positive_energy
from . import algebra


@pytest.fixture(scope='module')
def member():
    return contract.member(6)


def test_centered_probes_use_actual_bulk_cells_and_source_charge_spectrum(member):
    center = algebra.probe(member)
    assert center['nearest_lower_owner_cell'] == [4096]*3
    assert center['actual_centered_anchor'] == ['0']*3
    assert center['active']
    assert center['sequence_uniform_norm_upper'] == '1'
    assert center['commutes_all_local_Gauss_constraints']
    assert not center['stationary_algebra_includes_idle_reference']
    absent = algebra.probe(member, charge=10000)
    assert not absent['active'] and absent['zero_reason'] == 'absent_charge'
    early = algebra.probe(contract.member(0))
    assert not early['active'] and early['zero_reason'] == 'outside_interior'
    remote = algebra.probe(member, position=('1/3', '-1/7', '2/5'), kind='plaquette')
    assert remote['active'] and len(remote['unsmeared_owner_support']) == 3
    a = F(member['source']['spacing'])
    for target, actual in zip(remote['centered_target'], remote['actual_centered_anchor']):
        assert 0 <= F(target)-F(actual) < a


def test_probe_descriptor_rejects_forged_support_and_clock(member):
    value = algebra.probe(member)
    algebra.verify_probe(member, value)
    forged = deepcopy(value)
    forged['unsmeared_owner_support'][0][0] += 1
    with pytest.raises(ValueError):
        algebra.verify_probe(member, forged)
    forged_member = deepcopy(member)
    forged_member['clock']['Heisenberg'] = 'exp(2itH) A exp(-2itH)'
    with pytest.raises(ValueError):
        algebra.probe(forged_member)


def test_star_word_bounds_include_zero_factors_identity_and_all_translates(member):
    assert algebra.word_bound([]) == {'norm_upper': 1, 'time_Lipschitz_upper': 0, 'degree': 0}
    zero = algebra.word_bound([(0, 1), (3, 2)])
    assert zero['norm_upper'] == zero['time_Lipschitz_upper'] == 0
    word = algebra.word_bound([(2, 3), (3, 5), (7, 11)])
    assert word['norm_upper'] == 42
    assert word['time_Lipschitz_upper'] == 84*(F(1, 3)+F(1, 5)+F(1, 11))
    tail = algebra.word_negative_tail(member, [(2, 3), (3, 5)], '1/4')
    assert F(tail['spectral_bound']['original_negative_spectral_weight_upper']) == F(36, 256)
    assert not tail['spectral_bound']['original_weights_normalized_by_probe_success']


def test_filtered_cone_uses_temporal_halos_in_the_actual_same_clock(member):
    left = algebra.probe(member, half_width='1/100000')
    right = algebra.probe(member, position=('1', '0', '0'), half_width='1/100000')
    bound = algebra.filtered_locality_bound(member, left, right)
    assert bound['graph_owner_distance'] == 64
    assert bound['total_temporal_reach'] == F(1, 50000)
    assert bound['commutator_norm_upper'] < F(1, 2**28)
    widened = algebra.probe(member, position=('1', '0', '0'), half_width='1')
    inside = algebra.filtered_locality_bound(member, left, widened)
    assert inside['rational_cone_margin'] < 0 and inside['commutator_norm_upper'] == 2
    zero = algebra.filtered_locality_bound(member, left, algebra.probe(member, charge=10000))
    assert zero['exact_zero'] and zero['commutator_norm_upper'] == 0


def test_polynomial_support_error_proves_vanishing_without_float_exponentials():
    values = [algebra.polynomial_support_remainder(F(1, 2**n), 1, '1/50000', 3, 3)
              for n in (8, 10, 12)]
    assert values[0]['explicit_vanishing_upper'] == 4*values[1]['explicit_vanishing_upper']
    assert values[1]['explicit_vanishing_upper'] == 4*values[2]['explicit_vanishing_upper']
    with pytest.raises(ValueError):
        algebra.polynomial_support_remainder('1/2', 1, 0)
    with pytest.raises(ValueError):
        algebra.polynomial_support_remainder('1/1024', 1, 1)


def test_exact_filter_amplitude_from_direct_integral():
    x = s.symbols('x', real=True)
    direct = 2*s.integrate((1-x)*s.cos(s.pi*x), (x, 0, 1))
    control = algebra.exact_two_level_control()
    assert s.simplify(direct-control['filter_amplitude']) == 0
    P = control['unfiltered_projection']
    assert P*P == P
    assert s.simplify(s.det(control['filtered'])) > 0
    assert control['filtered']*control['filtered'] != control['filtered']


def test_complete_star_word_gram_and_canonical_GNS_are_exact():
    control = algebra.exact_two_level_control()
    V = control['gram_factor']
    G = control['word_gram']
    assert V.rank() == 2
    assert (G-V.H*V).applyfunc(s.simplify) == s.zeros(6)
    H = control['GNS_generator']; omega = control['canonical_ground_vector']
    assert H*omega == s.zeros(2, 1)
    assert H.eigenvals() == {0: 1, s.pi: 1}
    assert control['positive_frequency_weight'] == 4/s.pi**4
    assert control['negative_spectral_weight'] == 0
    B, Bt = control['filtered'], control['translated']
    assert (B*Bt-Bt*B).applyfunc(s.simplify) != s.zeros(2)
    assert not control['control_is_primary_spatial_source']


def test_GNS_null_space_is_a_left_ideal_not_a_Cstar_quotient_ideal():
    control = algebra.exact_two_level_control()
    assert control['null_left_ideal_example_norm_squared'] == 0
    assert control['same_element_right_product_norm_squared'] == 1
    omega = s.Matrix([1, 0]); a = s.Matrix([[0, 1], [0, 0]])
    for b in (s.eye(2), s.Matrix([[3, -2], [5, 7]])):
        assert b*a*omega == s.zeros(2, 1)


def test_uniform_filter_translation_bound_by_independent_quadrature():
    H = np.array([0., 1., 7.])
    A = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], complex)
    times = np.linspace(-1, 1, 4001)
    f = 1-np.abs(times)
    samples = np.exp(1j*times[:, None, None]*(H[:, None]-H[None, :]))*A
    B = np.trapezoid(f[:, None, None]*samples, times, axis=0)
    assert np.linalg.norm(B, 2) <= 1+1e-12
    for h in (.001, .03, .2):
        moved = np.exp(1j*h*(H[:, None]-H[None, :]))*B
        assert np.linalg.norm(moved-B, 2) <= 2*h+1e-6


def test_unsmoothed_energy_escape_and_filtered_collapse_are_not_hidden():
    escaping = positive_energy.escaping_energy_control(8)
    assert escaping['mean_excitation_energy'] > 4
    # A_n=|1><0|, H_n=diag(0,n), ground|0>. Its response is exp(int).
    # The original response has unit weight but is not equicontinuous.
    for n in (8, 32, 128):
        time = np.pi/n
        assert abs(np.exp(1j*n*time)-1) > 1.99
        filtered_amplitude = np.sinc(n/(2*np.pi))**2
        assert filtered_amplitude <= min(1, 4/n**2)+1e-15


def test_complete_finite_charge_measurement_does_not_force_limit_completeness():
    for n in (3, 8, 16):
        row = algebra.escaping_charge_control(n)
        assert row['complete_finite_probability_sum'] == 1
        assert all(p == 0 for p in row['fixed_charge_probabilities'].values())
        assert not row['limiting_fixed_charge_PVM_completeness_follows']
        assert row['state_charge'] in row['finite_PVM_labels']


def test_complete_two_instrument_pointer_keeps_every_history_and_coherence():
    row = algebra.complete_pointer_control()
    K0, K1 = row['Kraus']
    assert (K0*K0-row['effect']).applyfunc(s.simplify) == s.zeros(2)
    assert (K1*K1-(s.eye(2)-row['effect'])).applyfunc(s.simplify) == s.zeros(2)
    assert row['completeness'] == s.eye(2)
    density = row['pointer_density']; V = row['pointer_Gram_factor']
    assert (density-V.H*V).applyfunc(s.simplify) == s.zeros(4)
    assert s.trace(density) == 1
    assert sum(row['all_four_record_probabilities']) == 1
    assert all(p > 0 for p in row['all_four_record_probabilities'])
    assert any(density[i, j] != 0 for i in range(4) for j in range(i))
    assert not row['postselection'] and row['failure_outcome_retained']


def test_typed_payload_keeps_primary_and_comparison_controls_separate(member):
    receipt = algebra.bind_member(member)
    algebra.verify_bound(member, receipt)
    value = receipt['result']
    assert value['one_subsequence_for_all_finite_star_moments']
    assert value['GNS_source_energy_generator_nonnegative']
    assert not value['faithful_nontrivial_field_or_nonzero_energy_limit_proved']
    assert not value['arbitrary_idle_reference_observables_in_stationary_algebra']
    assert not value['perfect_finite_interlevel_state_embeddings_assumed']
    assert value['all_fixed_finite_complete_effect_history_probabilities_converge']
    assert not value['growing_canonical_reference_state_limit_claimed']
    forged = deepcopy(receipt)
    forged['result']['Lorentz_symmetry_interactions_or_particles_proved'] = True
    forged['result_sha256'] = contract.digest(forged['result'])
    with pytest.raises(ValueError):
        algebra.verify_bound(member, forged)
    json.dumps(receipt, allow_nan=False)


@pytest.mark.parametrize('kwargs', [
    {'position': (0., 0, 0)}, {'half_width': '0'}, {'kind': 'fermion'},
    {'axis': True}, {'charge': F(1, 2)}, {'axes': (2, 1)},
])
def test_invalid_or_undeclared_generator_labels_are_rejected(member, kwargs):
    with pytest.raises(ValueError):
        algebra.probe(member, **kwargs)
