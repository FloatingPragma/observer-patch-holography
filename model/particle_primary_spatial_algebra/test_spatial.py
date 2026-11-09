"""Exact controls for smearing, boundary transfer, and the momentum premise."""
from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path

import pytest
import sympy as s

from model.particle_common_contract import contract
from model.particle_primary_spatial_algebra import spatial


@pytest.fixture(scope='module')
def member():
    return contract.member(4)


def test_exact_hat_cells_and_normalization():
    # Independent antiderivative controls include both linear pieces and tails.
    assert spatial.hat_integral('-2', '-1') == 0
    assert spatial.hat_integral('-1', '0') == F(1, 2)
    assert spatial.hat_integral('0', '1/2') == F(3, 8)
    assert spatial.hat_integral('1/2', '3') == F(1, 8)
    for a, center, radius in [('1/2', ('0','0','0'), '1'),
                              ('2/3', ('1/7','-3/5','2/9'), '4/5')]:
        weights = spatial.hat_weights(a, center, radius)
        assert sum(weights.values()) == 1
        assert all(v > 0 for v in weights.values())
    assert len(spatial.hat_weights()) == 64


def test_cell_weights_exactly_intertwine_integer_lattice_shifts():
    a = F(1, 3)
    shift = (2, -3, 1)
    origin = (F(1,7), F(-2,5), F(0))
    translated = tuple(x+a*k for x,k in zip(origin, shift))
    first = spatial.hat_weights(a, origin, '2/3')
    second = spatial.hat_weights(a, translated, '2/3')
    shifted_first = {tuple(x+k for x,k in zip(cell, shift)): w
                     for cell,w in first.items()}
    assert second == shifted_first


def test_exact_discrete_smear_difference_obeys_L1_derivative_bound():
    # The exact one-dimensional L1 difference of two unit tents at d<=2
    # is 2d-d²/2. Cell integration contracts this distance.
    d = F(1, 5)
    first = spatial.hat_weights('1/4')
    second = spatial.hat_weights('1/4', (d, '0', '0'))
    distance = sum(abs(first.get(k,0)-second.get(k,0)) for k in first.keys() | second.keys())
    assert 0 < distance <= 2*d-d*d/2
    assert distance <= spatial.spatial_translation_bound((d,'0','0'))
    assert spatial.spatial_translation_bound(('100','0','0')) == 2


def test_spatial_allocation_and_bad_exact_inputs_fail():
    with pytest.raises(ValueError, match='allocation cap'):
        spatial.hat_weights('1/128', max_cells=1000)
    for args in [('1','0'), ('0','1','0','0')]:
        with pytest.raises(ValueError):
            spatial.hat_integral(*args)
    with pytest.raises(ValueError):
        spatial.rounded_translation(2, ('0','0'))
    with pytest.raises(ValueError):
        spatial.rounded_translation(2, (0.1,'0','0'))


def test_rounded_translations_have_a_nonzero_finite_group_defect():
    first = ('1/3','-2/5','1/7')
    second = ('2/7','1/9','-1/6')
    old_cap = None
    nonzero = False
    for n in range(1,9):
        row = spatial.rounded_translation(n, first)
        assert all(0 <= e < F(1,2**n) for e in row['coordinate_error'])
        defect = spatial.translation_group_defect(n, first, second)
        assert all(k in (-1,0) for k in defect['lattice_group_defect'])
        nonzero |= any(defect['lattice_group_defect'])
        cap = defect['vanishing_upper']
        assert defect['smeared_operator_norm_difference_upper'] <= cap
        if old_cap is not None: assert 2*cap == old_cap
        old_cap = cap
        assert defect['finite_rounded_maps_are_exact_R3_group'] is False
    assert nonzero


@pytest.mark.parametrize('side',[2,3,4])
def test_all_crossing_interaction_counts(side):
    row = spatial.crossing_incidence(side)
    assert row['crossing_bonds'] == 6*side**2
    assert row['crossing_plaquettes'] == 12*side**2-3*side
    exact_norm_times_a = 30*row['crossing_bonds']+F(row['crossing_plaquettes'],7)
    assert exact_norm_times_a <= row['crossing_norm_times_spacing_upper']


def test_same_source_bulk_bound_is_nontrivial_and_zero_time_is_exact(member):
    row = spatial.bulk_transfer(member)
    assert row['centered_collar_lower'] == 238
    assert row['physical_bulk_support_admitted']
    assert row['operator_difference_norm_upper']['negative_exponent'] > 50
    zero = spatial.bulk_transfer(member, duration='0')
    assert zero['operator_difference_norm_upper']['exact_zero']
    assert row['finite_physical_algebra_embedding_assumed_isometric'] is False
    assert row['infinite_state_assumed_stationary_before_limit'] is False
    assert row['full_spatial_matrices_executed'] is False


def test_small_box_cannot_claim_bulk_support_or_positive_exponent():
    row = spatial.bulk_transfer(contract.member(0))
    assert row['physical_bulk_support_admitted'] is False
    assert row['operator_difference_norm_upper']['negative_exponent'] == -1


def test_bulk_exponential_certificate_is_a_conservative_exact_bound(member):
    row = spatial.bulk_transfer(member)
    margin = F(row['rational_exponential_margin'])
    factor = F(row['Duhamel_prefactor'])
    bits = row['operator_difference_norm_upper']['negative_exponent']
    # exp(-m)<=2^-floor(m); round the prefactor upward independently.
    ceil_log = 0
    while F(2)**ceil_log < factor: ceil_log += 1
    assert bits <= margin.__floor__()-ceil_log
    assert F(2)**ceil_log >= factor


def test_time_positive_state_need_not_implement_spatial_translations():
    row = spatial.noninvariant_state_control()
    assert row['H'] == s.zeros(4)
    assert row['rho']**2 == row['rho']
    assert s.trace(row['rho']) == 1
    assert row['original_expectation'] == 1
    assert row['translated_expectation'] == 0
    assert row['null_vector_norm_squared'] == 0
    assert row['translated_null_norm_squared'] == 1
    assert row['canonical_spatial_GNS_implementation_exists'] is False


def test_finite_cyclic_relabelling_is_not_open_clock_symmetry():
    row = spatial.open_shift_counterexample()
    U = row['cyclic_permutation']; H = row['H_open']
    assert U*U.H == s.eye(3)
    assert U*H != H*U
    assert row['squared_Hilbert_Schmidt_defect'] == 4
    assert row['cyclic_label_is_open_source_symmetry'] is False


def test_bound_rejects_forged_momentum_and_source(member):
    receipt = spatial.bind_member(member)
    spatial.verify_bound(member, receipt)
    forged = deepcopy(receipt)
    forged['result']['canonical_GNS_momentum_or_vacuum_momentum_proved'] = True
    forged['result_sha256'] = contract.digest(forged['result'])
    with pytest.raises(ValueError, match='receipt differs'):
        spatial.verify_bound(member, forged)
    altered = deepcopy(member)
    altered['source']['spacing'] = '1/32'
    with pytest.raises(ValueError):
        spatial.bulk_transfer(altered)
    json.dumps(receipt, allow_nan=False)


def test_sources_exist_and_no_physical_mass_is_exported(member):
    for path in spatial.SOURCE_INPUTS: assert Path(path).is_file(), path
    payload = spatial.common_payload(member)
    assert payload['regular_algebra_has_point_norm_continuous_R3_action'] is True
    assert payload['selected_limiting_state_spatial_invariance_proved'] is False
    assert payload['physical_zero_momentum_rest_energy_or_mass_proved'] is False
    assert payload['extra_spatial_averaging_coin_inserted'] is False
