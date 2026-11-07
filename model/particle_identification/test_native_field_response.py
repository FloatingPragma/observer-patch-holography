"""Reject cross-source, cross-sector and normalization substitutions."""
from copy import deepcopy
from fractions import Fraction as F

import pytest

from tools import strict_json
from ..native_correlation_preparation import packet as prep
from ..native_exchange_clock import packet as clock
from .native_field_response import compose


def retained():
    return strict_json.load(prep.OUTPUT), strict_json.load(clock.OUTPUT)


def test_composition_retains_same_sector_and_native_resource_boundary():
    result = compose(*retained())
    assert len(result['responses']) == 9
    for row in result['responses']:
        expected = F(row['normalized_preparation_error_upper']) + F(row['normalized_clock_error_upper'])
        assert F(row['combined_normalized_reference_score_error_upper']) == expected
        assert row['complete_attempts'] == row['warmup_attempts'] + row['lag_attempts'] + 4
        assert row['single_program_admitted'] == (row['level'] < 2)
        assert bool(row['complete_program_digest']) == row['single_program_admitted']
        assert row['expanded_history_executed'] is False
    assert result['physical_particle_identification'] is False


@pytest.mark.parametrize('field,value', [
    ('ports', 241), ('seams', 631), ('ports', 240.0),
    ('probe_seams', [[1, 82], [0, 6]]), ('level', False),
    ('native_time', '1/1000')])
def test_mismatched_source_is_rejected(field, value):
    prepared, response = retained(); response = deepcopy(response)
    response['evaluated_federations'][0][field] = value
    with pytest.raises(ValueError):
        compose(prepared, response)


@pytest.mark.parametrize('field,value', [('global_occupation', 3),
                                       ('variance', '1/33')])
def test_mismatched_occupation_or_variance_is_rejected(field, value):
    prepared, response = retained(); response = deepcopy(response)
    response['evaluated_federations'][0]['occupations'][0][field] = value
    with pytest.raises(ValueError):
        compose(prepared, response)


def test_mismatched_preparation_normalization_is_rejected():
    prepared, response = retained(); prepared = deepcopy(prepared)
    prepared['federation_preparation_budgets'][0]['budgets'][0]['normalized_score_error_upper'] = '1/10000'
    with pytest.raises(ValueError):
        compose(prepared, response)
