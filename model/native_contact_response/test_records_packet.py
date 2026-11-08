from copy import deepcopy
from fractions import Fraction as F
import hashlib
import json

import pytest

from . import contact as c, records, packet


def test_native_complete_records_and_finite_warmup():
    result = records.certificate()
    assert result['expanded_control_status'] == 'RESOURCE_STOP'
    assert result['preparation']['warmup_attempts'] == 480
    assert result['preparation']['attempts_upper'] == 488
    assert result['preparation']['normalized_score_error_upper'] == '33/65536'
    assert [row['normalized_response'] for row in result['exact_scored_records'][:3]] == [F(1), F(4, 5), F(577, 900)]
    assert result['complete_probe_pair_slots'] == 625
    assert result['uniform_attempt_slots'] == 721


def test_wrong_reflection_cut_breaks_actual_response(monkeypatch):
    original = records.scored_step
    def wrong(value, command, density, side=None):
        return original(value, command, density, 'incoming' if side == 'outgoing' else side)
    monkeypatch.setattr(records, 'scored_step', wrong)
    with pytest.raises(ArithmeticError, match='scored records'):
        records.exact_scores(1)


@pytest.mark.parametrize('value', [True, -1, 0.0, '1'])
def test_program_and_sparse_query_count_validation(value):
    with pytest.raises(ValueError):
        records.program(value)
    with pytest.raises(ValueError):
        records.exact_scores(value)


def test_attempt_ceiling_and_finite_query_cap_remain_distinct():
    with pytest.raises(ValueError):
        records.program(records.interpreter.MAX_ATTEMPTS-3)
    with pytest.raises(ValueError):
        records.exact_scores(13)
    assert records.program(100)['declared_resources']['attempts_upper'] == 104


def test_sampling_budget_retains_low_density_normalization():
    small = records.sampling_budget(12, F(1, 100), F(1, 1000))
    large = records.sampling_budget(240, F(1, 100), F(1, 1000))
    assert small['raw_variance'] == F(2, 33)
    assert large['raw_variance'] == F(1, 7170)
    assert large['independent_complete_trials'] > 100000*small['independent_complete_trials']
    assert large['independence_and_repeated_preparation_are_required']
    with pytest.raises(ValueError):
        records.sampling_budget(12, F(1), F(1))


@pytest.fixture(scope='module')
def receipt():
    return packet.produce()


def test_retained_receipt_matches_full_source_replay(receipt):
    assert packet.canonical(json.loads(packet.OUTPUT.read_text())) == packet.canonical(receipt)
    assert receipt['exact_response']['observed_cyclic_dimension'] == 12
    assert receipt['federation_leakage']['distinct_boundary_triplets'] == 220
    assert receipt['scope']['correlated_mode_is_a_localized_physical_particle'] is False


def test_source_bindings_are_fresh_and_exclude_docs(receipt):
    assert receipt['source_bindings'] == packet.source_bindings()
    assert all(not name.startswith('docs/') for name in receipt['source_bindings'])
    for name, digest in receipt['source_bindings'].items():
        assert digest == hashlib.sha256((packet.ROOT/name).read_bytes()).hexdigest()


@pytest.mark.parametrize('mutation', ['float', 'bool', 'extra', 'source', 'pole'])
def test_typed_receipt_mutations_fail(tmp_path, monkeypatch, receipt, mutation):
    value = deepcopy(receipt)
    if mutation == 'float':
        value['source']['ports'] = 12.0
    elif mutation == 'bool':
        value['exact_response']['probe_seams'][0][0] = False
    elif mutation == 'extra':
        value['unproved_particle'] = True
    elif mutation == 'source':
        value['source_bindings'][next(iter(value['source_bindings']))] = '0'*64
    else:
        value['exact_response']['lowest_pole']['pole_multiplicity_in_primitive_sector'] = 1
    path = tmp_path/'bad.json'; path.write_text(packet.canonical(value))
    monkeypatch.setattr(packet, 'produce', lambda: receipt)
    with pytest.raises(ValueError):
        packet.verify(path)


def test_duplicate_json_keys_rejected_before_receipt_comparison(tmp_path):
    path = tmp_path/'duplicate.json'; path.write_text('{"schema":1,"schema":1}')
    with pytest.raises(ValueError):
        packet.verify(path)


def test_outward_receipt_intervals_enclose_rational_proof():
    proof = c.pole_certificate()
    for key in ('pole_interval', 'normalized_residue_interval'):
        lo, hi = packet.outward(proof[key])
        assert lo <= proof[key][0] < proof[key][1] <= hi
