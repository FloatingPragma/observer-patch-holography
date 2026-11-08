"""Admission boundaries and live source classification checks."""
import copy
import subprocess
import sys

import pytest

from . import packet, pauli


@pytest.fixture(scope='module')
def produced():
    return packet.build()


def test_retained_catalogue_replays(produced):
    retained = packet.strict_json.load(packet.HERE/'catalog.json')
    assert packet.canonical(retained) == packet.canonical(produced)
    assert (packet.HERE/'CATALOG.md').read_text() == packet.render(produced)
    for law in ('C', 'E'):
        assert produced['bare_complete_branch_classification'][law]['full78']['commutant_complex_dimension'] == 1
        assert produced['bare_complete_branch_classification'][law]['binary66']['commutant_complex_dimension'] == 1
        assert produced['ordered_cover_complete_branch_classification'][law]['commutant_complex_dimension'] == 2


def test_source_matrix_mutation_rejected(monkeypatch):
    original = pauli.source.matrix
    def changed():
        h, scale = original()
        h = h.copy(); h[0, 0] += 1
        return h, scale
    monkeypatch.setattr(pauli.source, 'matrix', changed)
    with pytest.raises(ArithmeticError, match='actual source Hamiltonian'):
        pauli.source_polynomial()


@pytest.mark.parametrize('forgery', ['physical_mass','global_completeness','species','Pauli_count'])
def test_claim_promotion_rejected(produced, monkeypatch, forgery):
    # Independent immutable reference avoids another expensive full source
    # replay for each pure admission mutation.
    expected = copy.deepcopy(produced)
    monkeypatch.setattr(packet, 'build', lambda: expected)
    changed = copy.deepcopy(produced)
    if forgery == 'physical_mass':
        changed['entries'][2]['physical_mass'] = 1
    elif forgery == 'global_completeness':
        changed['catalogue_is_exhaustive_over_all_OPH_sources_or_defects'] = True
    elif forgery == 'species':
        changed['entries'][2]['identified_Standard_Model_species'] = ['electron']
    else:
        changed['selected_spin_complete_Pauli_classification']['selected_field_charge_words_modulo_phase'].append([0, 1])
    with pytest.raises(ValueError, match='catalogue classification'):
        packet.validate(changed)


def test_no_double_counting_or_physical_assignments(produced):
    entries = produced['entries']
    assert len({row['id'] for row in entries}) == len(entries)
    assert all(not row['identified_Standard_Model_species'] and row['physical_mass'] is None for row in entries)
    assert not produced['physical_particle_map']
    assert len(next(row for row in entries if row['id']=='ordered_exchange_charge')['same_object_realizations']) == 2


def test_joint_family_has_literal_observables_and_distinct_state_clock(produced):
    entry=next(row for row in produced['entries'] if row['id']=='joint_positive_field_canonical_sector')
    attached=produced['joint_positive_field_source_attachment']
    assert attached['catalogue_entry']==entry['id']
    assert attached['canonical_pairs']==[['Q3','P2'],['Q2','P3']]
    assert attached['source_observables']['P2']=='-p sum G(px) S_x^2'
    assert attached['source_observables']['P3']=='+p sum G(px) S_x^3'
    assert not attached['ordered_zero_field_history_identity_proved']
    assert not attached['repair_identity_generates_the_H_B_clock']
    assert attached['owning_scientific_producer_replayed']
    assert attached['complete_finite_history_example']['all_native_failure_and_pointer_records_retained']
    assert not attached['physical_species'] and not attached['physical_masses']


@pytest.mark.parametrize('forgery',['missing_second_pair','old_limit_identity','physical_species','unread_selection'])
def test_joint_owner_premise_or_scope_mutation_rejected(monkeypatch,forgery):
    original=packet.strict_json.load(packet.joint_field.OUTPUT)
    changed=copy.deepcopy(original)
    key={'missing_second_pair':'second_canonical_pair_on_same_joint_history_subsequence',
         'old_limit_identity':'equality_to_ordered_zero_field_state_clock_limit_proved',
         'physical_species':'physical_particle_assignment',
         'unread_selection':'sequence_forced_by_unread_repair'}[forgery]
    changed['scope'][key]=False if forgery=='missing_second_pair' else ('electron' if forgery=='physical_species' else True)
    # The owner is a separate scientific admission; neither changed output nor
    # matching retained bytes can transfer its unproved physical premises.
    monkeypatch.setattr(packet.joint_field,'build',lambda:changed)
    monkeypatch.setattr(packet.strict_json,'load',lambda path:changed)
    with pytest.raises(ValueError,match='joint-field'):
        packet.joint_field_attachment()


def test_fresh_full_producer_binds_runtime_imports():
    script = '''
from pathlib import Path
import sys
from model.protected_defect_catalog import packet
result = packet.build()
missing = []
for module in tuple(sys.modules.values()):
    filename = getattr(module, '__file__', None)
    if filename:
        path = Path(filename).resolve()
        if path.is_relative_to(packet.ROOT) and path.suffix == '.py':
            name = path.relative_to(packet.ROOT).as_posix()
            if name not in result['source_bindings']:
                missing.append(name)
assert not missing, sorted(set(missing))
'''
    subprocess.run([sys.executable, '-c', script], cwd=packet.ROOT,
                   capture_output=True, text=True, check=True)
