"""Complete trial retention, exact interval bounds and physical scope."""
from copy import deepcopy

import pytest

from tools import strict_json
from . import packet


@pytest.fixture(scope='module')
def fresh():
    return packet.report()


def test_complete_retained_catalogue_has_a_fresh_certificate(fresh):
    assert packet.canonical(strict_json.load(packet.OUTPUT)) == packet.canonical(fresh)
    assert fresh['case_count'] == 28
    assert fresh['exact_two_zero_kernel_cases'] == 16
    assert fresh['certified_invertible_cases'] == fresh['unresolved_kernel_cases'] == 6
    for row in fresh['complete_cases']:
        assert len(row['ordered_energy_intervals']) == row['dimension']
        assert len(row['full_Fock_contact_cases']) == 3


@pytest.mark.parametrize('mutation', ['missing', 'energy', 'index', 'scope', 'source', 'type'])
def test_claim_omissions_and_promotions_fail_closed(fresh, monkeypatch, mutation):
    monkeypatch.setattr(packet, 'report', lambda: fresh)
    value = deepcopy(fresh)
    if mutation == 'missing': value['complete_cases'].pop()
    elif mutation == 'energy': value['complete_cases'][0]['ordered_energy_intervals'][0][0] = '0'
    elif mutation == 'index': value['complete_cases'][0]['exact_kernel_dimension'] = 6
    elif mutation == 'scope': value['scope']['observed_particle_identification'] = True
    elif mutation == 'source': value['source_bindings'].pop(next(iter(value['source_bindings'])))
    else: value['scope']['all_finite_source_eigenvalues_retained'] = 1
    with pytest.raises(ValueError, match='fresh certification'):
        packet.validate(value)
