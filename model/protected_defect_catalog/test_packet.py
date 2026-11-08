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
