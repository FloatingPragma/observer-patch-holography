"""Candidate coverage, ancestry and scope across independent calculations."""
from copy import deepcopy
import json

import pytest

from . import model as m


@pytest.fixture(scope='module')
def value():
    return m.build()


def test_all_target_and_source_ids_are_preserved(value):
    source = json.loads((m.ROOT / m.INPUTS['dictionary']).read_text())
    assert [r['target'] for r in value['target_rows']] == [r['target'] for r in source['target_rows']]
    assert value['source_excitation_ids'] == [r['id'] for r in source['source_excitation_rows']]
    assert value['protected_catalogue_ids'] == [r['id'] for r in source['protected_catalogue_rows']]


def test_generations_share_family_tests_without_acquiring_flavors(value):
    rows = {r['target']: r for r in value['target_rows']}
    for family in (('electron', 'muon', 'tau'), ('up', 'charm', 'top'),
                   ('down', 'strange', 'bottom'),
                   ('electron_neutrino', 'muon_neutrino', 'tau_neutrino')):
        assert len({rows[name]['family'] for name in family}) == 1
        assert all(rows[name]['test_receipt'] == m.INPUTS['charged'] for name in family)
    assert all(rows[name]['test_receipt'] == m.INPUTS['bosons']
               for name in ('photon', 'gluon', 'W_plus', 'W_minus', 'Z', 'Higgs'))
    assert all(row['source_flavor_index'] is None for row in rows.values())


@pytest.mark.parametrize('field', ['physical_particle_map', 'mass_predictions_GeV',
                                 'fit_parameters_used', 'physical_gate_promotions'])
def test_unproved_scope_promotions_fail(value, field):
    changed = deepcopy(value)
    changed[field] = ['electron']
    with pytest.raises(ValueError, match='physical particles'):
        m.validate(changed)


def test_forged_nested_consequence_fails():
    nested = {'selected_member_bound': {'result': {'mass': 1}, 'result_sha256': '0' * 64}}
    with pytest.raises(ValueError, match='consequence digest'):
        m.validate_result_digests(nested)


def test_changed_producer_custody_fails(monkeypatch):
    monkeypatch.setattr(m, 'sha', lambda _: '0' * 64)
    with pytest.raises(ValueError, match='stale candidate producer'):
        m.read_inputs()


def test_same_source_binding_rejects_a_different_clock():
    member = m.contract.member(0)
    result = {'member_id': member['member_id'], **member['ids']}
    m.validate_primary_binding({'selected_member_bound': {'result': result}})
    result['clock_id'] = 'different-clock'
    with pytest.raises(ValueError, match='source, state, reference or clock'):
        m.validate_primary_binding({'selected_member_bound': {'result': result}})


def test_relabeling_a_family_cannot_redirect_its_evidence(value):
    changed = deepcopy(value)
    row = next(r for r in changed['target_rows'] if r['target'] == 'electron')
    row.update(family='photon', test_receipt=m.INPUTS['bosons'],
               source_test_summary=m.SOURCE_TESTS['photon'])
    with pytest.raises(ValueError, match='wrong source owner'):
        m.validate(changed)


def test_machine_and_readable_coverage_agree(value):
    table = m.render(value)
    for row in value['target_rows']:
        assert '| ' + row['target'] + ' |' in table
    assert value['physical_particle_map'] == value['mass_predictions_GeV'] == []


def test_calculated_detailed_results_cannot_drift_from_the_receipt(value):
    changed = deepcopy(value)
    changed['calculated_source_results']['composites']['cyclic_radix'] = 5
    with pytest.raises(ValueError, match='result content'):
        m.validate(changed)


def test_resealed_summary_cannot_replace_authenticated_calculation(value):
    changed = deepcopy(value)
    changed['calculated_source_results']['composites']['cyclic_radix'] = 5
    changed['content_sha256'] = m.content_digest(changed)
    with pytest.raises(ValueError, match='authenticated producers'):
        m.validate(changed)


def test_member_order_does_not_change_the_selected_calculation():
    values, _ = m.read_inputs()
    expected = m.calculated_summary(values)
    values['charged']['members'].reverse()
    assert m.calculated_summary(values) == expected


def test_cross_family_support_cannot_be_removed(value):
    changed = deepcopy(value)
    row = next(r for r in changed['target_rows'] if r['target'] == 'electron_neutrino')
    row['supporting_receipts'] = []
    changed['content_sha256'] = m.content_digest(changed)
    with pytest.raises(ValueError, match='supporting evidence'):
        m.validate(changed)
