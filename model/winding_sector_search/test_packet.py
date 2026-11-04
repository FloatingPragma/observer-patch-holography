"""Retain all trials, eigenvalues, source custody and physical scope."""
from copy import deepcopy

import pytest

from tools import strict_json
from . import packet


@pytest.fixture(scope='module')
def fresh():
    return packet.report()


def test_retained_receipt_matches_complete_fresh_search(fresh):
    packet.compare(strict_json.load(packet.OUTPUT), fresh)
    search = fresh['finite_search']
    assert search['completed_case_count'] == 28
    assert len(search['all_initial_cases']) == 24
    assert search['volume_followup']['chosen_after_initial_search'] is True
    for row in search['all_initial_cases']+search['volume_followup']['all_cases']:
        assert len(row['all_eigenvalues']) == row['dimension']
        assert row['dimension'] == 2*(2*row['radius']+1)**3


@pytest.mark.parametrize('mutation', ['trial', 'level', 'index', 'scope', 'type', 'source', 'nan', 'threshold'])
def test_cannot_omit_candidates_change_index_or_promote_particle_assignment(fresh, monkeypatch, mutation):
    monkeypatch.setattr(packet, 'report', lambda: fresh)
    damaged = deepcopy(fresh)
    if mutation == 'trial':
        damaged['finite_search']['all_initial_cases'].pop()
    elif mutation == 'level':
        damaged['finite_search']['all_initial_cases'][0]['all_eigenvalues'][0] += .1
    elif mutation == 'index':
        damaged['degree_family'][0]['source_continuum_chiral_index'] = -4
    elif mutation == 'scope':
        damaged['scope']['same_current_particle_assignment'] = True
    elif mutation == 'type':
        damaged['scope']['all_source_spin_and_valley_rows'] = 1
    elif mutation == 'source':
        damaged['source_bindings'].pop(next(iter(damaged['source_bindings'])))
    elif mutation == 'threshold':
        damaged['finite_search']['all_initial_cases'][0]['numerical_zero_threshold'] = 1.005e-8
    else:
        damaged['finite_search']['all_initial_cases'][0]['all_eigenvalues'][0] = float('nan')
    with pytest.raises(ValueError, match='fresh calculation'):
        packet.validate(damaged)
