"""Exact export/marking controls; no spatial quantum preparation is supplied."""
from copy import deepcopy
from fractions import Fraction as Q
from pathlib import Path
import sys
from unittest.mock import patch

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import proper_code_attachment_probe as producer
import verify_proper_code_attachment_independent as verifier


@pytest.fixture(scope='module')
def packet():
    return verifier.m.load(producer.REFERENCE)


def reject(packet):
    with pytest.raises(ValueError):
        verifier.verify(packet)


def test_source_producer_reproduces_frozen_probe(packet):
    assert producer.build() == packet


def test_independent_replay_recognizes_declared_internal_positive_control(packet):
    assert verifier.verify(packet)['verified']
    four = packet['services']['d4_coherent']
    assert four['upper_identity_including_failure']
    assert four['success_effect_commuting_action_count'] == 60
    assert four['outcomes'] == ['success', 'failure']
    assert four['kraus_counts'] == {'success': 1, 'failure': 2}
    assert four['full_M6_success_effect_commuting_action_count'] == 4
    assert four['full_M6_fixed_mark_probability_witness']['processor_basis_index'] == 3
    assert packet['mathematical_attachment']['hom_slow_to_upper'] == 1
    assert packet['mathematical_attachment']['hom_slow_to_lower'] == 0


def test_mathematical_inverse_and_source_scale_are_exact(packet):
    math = packet['mathematical_attachment']
    t = [[producer.F(Q(a), Q(b)) for a, b in row] for row in math['T_upper_from_spatial_loads']]
    frozen = verifier.m.load(producer.ROOT/'code/a5_closure/manifests/marked_triplet_attachment_reference.json')
    p = [[producer.F(Q(a), Q(b)) for a, b in row] for row in frozen['projectors'][0]]
    assert math['TT_star_scale'] == ['10', '2']
    inverse = [[x/producer.F(10, 2) for x in row] for row in producer.transpose(t)]
    assert producer.mm(t, inverse) == producer.eye(3)
    assert producer.mm(inverse, t) == p


def test_d2_has_a_complete_rational_fixed_mark_obstruction(packet):
    row = packet['services']['d2_coherent']
    assert not row['upper_identity_including_failure']
    assert row['success_effect_commuting_action_count'] == 4
    assert row['fixed_mark_probability_witness'] == {
        'upper_basis_index': 0,
        'before': {'success': ['1', '0'], 'failure': ['0', '0']},
        'after': {'success': ['3/4', '0'], 'failure': ['1/4', '0']}}
    assert not packet['scope']['universal_compiler_obstruction']


def test_dephasing_keeps_the_effect_but_loses_the_surviving_reference_state(packet):
    coherent = packet['services']['d4_coherent']
    dephased = packet['services']['d4_dephased']
    assert dephased['success_effect_upper'] == coherent['success_effect_upper']
    assert not dephased['upper_identity_including_failure']
    assert coherent['entangled_reference_fixture_output'][0][0][3] == ['1/2', '0']
    assert dephased['entangled_reference_fixture_output'][0][0][3] == ['0', '0']


@pytest.mark.parametrize('attack', [
    'historical_mark', 'same_character_inverse', 'galois_intertwiner', 'wrong_scale',
    'lower_attachment', 'selected_sign', 'physical_selection', 'universal_obstruction',
    'spatial_state_encoding', 'compiled_word_claim', 'missing_failure', 'dimension',
    'reference_coherence', 'diagonal_only', 'source_commit', 'source_blob',
    'source_sha', 'source_byte_count', 'historical_hash_repair', 'receipt_extra_field',
])
def test_hostile_promotions_and_mutations_fail_closed(packet, attack):
    bad = deepcopy(packet)
    four = bad['services']['d4_coherent']
    if attack == 'historical_mark':
        bad['marking']['current_index'] = 1
        bad['marking']['current_traces'] = bad['marking']['historical_traces']
    elif attack == 'same_character_inverse':
        p = bad['marking']['native_action']
        bad['marking']['native_action'] = [p.index(i) for i in range(12)]
    elif attack == 'galois_intertwiner':
        for row in bad['mathematical_attachment']['T_upper_from_spatial_loads']:
            for pair in row:
                pair[1] = str(-Q(pair[1]))
    elif attack == 'wrong_scale':
        bad['mathematical_attachment']['TT_star_scale'] = ['1', '0']
    elif attack == 'lower_attachment':
        bad['mathematical_attachment']['hom_slow_to_lower'] = 1
    elif attack in ('selected_sign', 'physical_selection', 'universal_obstruction', 'spatial_state_encoding', 'compiled_word_claim'):
        key = {'selected_sign': 'target_character_used_for_selection',
               'physical_selection': 'physical_golden_embedding_selected',
               'universal_obstruction': 'universal_compiler_obstruction',
               'spatial_state_encoding': 'source_spatial_quantum_preparation_encoding_derived',
               'compiled_word_claim': 'compiled_marked_action_words_certified'}[attack]
        bad['scope'][key] = True
    elif attack == 'missing_failure':
        four['outcomes'].remove('failure')
    elif attack == 'dimension':
        four['output_dimension'] = 4
    elif attack == 'reference_coherence':
        four['entangled_reference_fixture_output'][0][0][3] = ['0', '0']
        four['entangled_reference_fixture_output'][0][3][0] = ['0', '0']
    elif attack == 'diagonal_only':
        four['all_input_matrix_units_checked'] = 6
    elif attack.startswith('source_'):
        key = {'source_commit': 'commit', 'source_blob': 'blob', 'source_sha': 'sha256', 'source_byte_count': 'bytes'}[attack]
        row = bad['source_pins'][producer.COMPILER]
        row[key] = row[key]+1 if key == 'bytes' else '0'*len(row[key])
    elif attack == 'historical_hash_repair':
        bad['historical_galois_drift']['hashes_updated'] = True
    else:
        bad['new_attachment'] = 'target chosen'
    reject(bad)


@pytest.mark.parametrize('mutation', ['range(d, 5)', 'range(0, 6)'])
def test_changed_source_export_is_rejected_even_before_receipt_comparison(mutation):
    text = (producer.ROOT/producer.COMPILER).read_text()
    text = text.replace('range(d, 6)', mutation)
    with pytest.raises(ValueError, match='AST contract'):
        verifier.exact_export_source(text)


def test_actual_compiler_missing_failure_branch_breaks_trace_preservation():
    original = producer.compiler.transfer_kraus
    with patch.object(producer.compiler, 'transfer_kraus', side_effect=lambda d, coherent: original(d, coherent)[:-1]):
        with pytest.raises(ValueError, match='trace preservation'):
            producer.service_report(4, True, [], 4)


def test_existing_whole_law_guard_replays_on_the_audited_prose(packet):
    # The guard must reject drift, not the audited text: every pin now matches,
    # and the independent pin controls reject any single wrong digest.
    import galois_source_response_control as historical
    drift = packet['historical_galois_drift']
    assert drift['mismatches'] == []
    assert sum(row['expected'] == row['actual'] for row in drift['pins'].values()) == 4
    result = historical.verify()
    assert result['plus'] == result['minus']


def test_duplicate_and_nonfinite_receipts_are_not_admitted(tmp_path):
    path = tmp_path/'bad.json'
    for text in ('{"verdict":0,"verdict":1}', '{"value":NaN}'):
        path.write_text(text)
        with pytest.raises(ValueError):
            verifier.m.load(path)
