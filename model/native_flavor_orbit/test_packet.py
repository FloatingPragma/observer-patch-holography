"""Cold source custody and complete scientific-payload rejection controls."""
from copy import deepcopy

import pytest

from tools import strict_json
from . import packet


@pytest.fixture(scope='module')
def fresh():
    return packet.report()


def test_stored_receipt_matches_complete_cold_calculation(fresh):
    packet.validate(strict_json.load(packet.OUTPUT), fresh)
    assert fresh['complete_current_commutant']['family_commutant_complex_dimension'] == 45
    assert len(fresh['current_gated_joint_spectrum']['policies']) == 4


@pytest.mark.parametrize('mutation', ['commutant', 'record', 'zero', 'CP', 'mixing', 'selection', 'gauge',
                                    'history_CP', 'noise', 'logarithm', 'coherent'])
def test_scientific_omissions_and_promotions_are_rejected(fresh, mutation):
    wrong = deepcopy(fresh)
    if mutation == 'commutant':
        wrong['complete_current_commutant']['family_commutant_complex_dimension'] = 54
    elif mutation == 'record':
        wrong['complete_memory_theorem']['complete_Kraus_label_count'] = 12
    elif mutation == 'zero':
        wrong['current_gated_joint_spectrum']['policies'][2]['full_K45_nullity'] = 3
    elif mutation == 'CP':
        wrong['current_gated_joint_spectrum']['policies'][0]['charged_current_CP_commutator_trace_cube_exact'] = '1'
    elif mutation == 'mixing':
        wrong['current_gated_joint_spectrum']['policies'][0]['quark_mixing_modulus_squared'][0][1] = ['0', '0']
    elif mutation == 'selection':
        wrong['scope']['source_selected_physical_couplings_or_particle_map'] = True
    elif mutation == 'gauge':
        wrong['scope']['gauge_covariance_alone_forces_zero_CP'] = True
    elif mutation == 'history_CP':
        wrong['finite_correlated_memory']['symmetry_is_physical_CP_identification'] = True
    elif mutation == 'noise':
        wrong['finite_correlated_memory']['stationary_interval_witnesses'][0]['normalized_Choi_purity_defect_interval'] = ['0', '0']
    elif mutation == 'logarithm':
        wrong['finite_correlated_memory']['channel_logarithms']['complete_CQ_channel_has_any_matrix_logarithm'] = True
    else:
        wrong['finite_correlated_memory']['exact_controls']['two_step_diagonal_A02_exact'] = '0'
    with pytest.raises(ValueError, match='complete fresh calculation'):
        packet.validate(wrong, fresh)


def test_stale_inherited_source_is_rejected_before_scientific_replay(monkeypatch):
    original = strict_json.load

    def corrupt(path):
        result = original(path)
        if str(path).endswith('native_family_protocol_results.json'):
            result['source_bindings']['model/common_source/native_gauge_family_records.py'] = '0'*64
        return result

    monkeypatch.setattr(strict_json, 'load', corrupt)
    with pytest.raises(ValueError, match='stale inherited source'):
        packet.source_bindings()
