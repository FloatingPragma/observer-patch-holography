"""Typed replay, source custody and complete-observation boundaries."""
from copy import deepcopy
from fractions import Fraction as F

import pytest

from . import packet as p


@pytest.fixture(scope='module')
def candidate():
    return p.produce()


def test_owner_replay_and_complete_source_boundaries(candidate):
    assert p.verify()
    assert candidate['scope']['global_density_trace_distance_bound_claimed'] is False
    assert candidate['scope']['physical_vacuum_or_particle_identification'] is False
    assert candidate['scope']['arbitrary_initial_quantum_state_and_untouched_reference'] is True
    for row in candidate['sparse_federation_cross_effect_controls']:
        assert row['degree_upper'] <= 4
    budgets = [row for graph in candidate['federation_preparation_budgets'] for row in graph['budgets']]
    assert any(row['single_program_admitted'] for row in budgets)
    assert any(not row['single_program_admitted'] for row in budgets)
    for row in budgets:
        assert F(row['normalized_score_error_upper']) <= F(1, 1024)
        assert row['complete_attempts'] == row['warmup_attempts']+row['lag_attempts']+4
        assert row['global_state_comparison_halving_blocks'] >= row['halving_blocks']


@pytest.mark.parametrize('change', ('bool', 'float', 'degree', 'error', 'limit', 'hash'))
def test_receipt_rejects_typed_scientific_and_custody_mutations(candidate, tmp_path, change):
    value = deepcopy(candidate)
    if change == 'bool':
        value['all_domain_theorem']['degree_upper'] = True
    elif change == 'float':
        value['all_domain_theorem']['degree_upper'] = 4.0
    elif change == 'degree':
        value['all_domain_theorem']['degree_upper'] = 2
    elif change == 'error':
        value['all_input_H66_prepared_response']['budget']['normalized_score_error_upper'] = '0'
    elif change == 'limit':
        value['scope']['all_budgets_fit_native_single_program_limit'] = True
    else:
        value['source_bindings']['model/native_correlation_preparation/preparation.py'] = '0'*64
    path = tmp_path/'changed.json'; path.write_text(p.canonical(value))
    with pytest.raises(ValueError):
        p.verify(path)


def test_duplicate_receipt_keys_rejected_before_replay(candidate, tmp_path):
    path = tmp_path/'duplicate.json'
    path.write_text('{"schema":"wrong",'+p.canonical(candidate)[1:])
    with pytest.raises(ValueError):
        p.verify(path)
