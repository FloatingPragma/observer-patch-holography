"""Exact receipt scope, edge counting and independent-process replay controls."""
from copy import deepcopy
import subprocess
import sys

import pytest

from . import packet


@pytest.fixture(scope='module')
def fresh():
    return packet.report()


def test_complete_source_receipt_and_scientific_boundary(fresh):
    packet.validate(fresh)
    assert len(fresh['exact_degree_and_edges']['exact_radial_shells']) == 8
    assert fresh['canonical_physical_premises_discharged'] == []
    assert not any(p.startswith('model/completion/') for p in fresh['source_bindings'])


@pytest.mark.parametrize('damage', ['degree', 'edge', 'energy', 'repair', 'lifetime', 'particle', 'type', 'source'])
def test_exact_scope_and_claim_mutations_fail_closed(fresh, monkeypatch, damage):
    monkeypatch.setattr(packet, 'report', lambda: fresh)
    value = deepcopy(fresh)
    if damage == 'degree':
        value['exact_degree_and_edges']['exact_radial_shells'][0]['degree'] = 0
    elif damage == 'edge':
        value['exact_degree_and_edges']['exact_radial_shells'][0]['source_edge_count'] -= 1
    elif damage == 'energy':
        value['same_source_energy']['charged_shell_energy_lower'] = '0'
    elif damage == 'repair':
        value['scope']['interleaved_Hamiltonian_evolution_is_repair_identity'] = True
    elif damage == 'lifetime':
        value['scope']['absolute_crossing_floor_is_lifetime_or_activation_barrier'] = True
    elif damage == 'particle':
        value['scope']['all_quantum_defects_or_particles_excluded'] = True
    elif damage == 'type':
        value['scope']['native_repair_selects_winding_or_spin_direction'] = 0
    else:
        value['source_bindings'].pop(next(iter(value['source_bindings'])))
    with pytest.raises(ValueError, match='fresh calculation'):
        packet.validate(value)


def test_independent_cold_replay(tmp_path, fresh):
    path = tmp_path/'lattice_texture.json'
    path.write_text(packet.canonical(fresh))
    result = subprocess.run([sys.executable, '-W', 'error', '-m', 'model.lattice_texture_selection',
                             '--check', str(path)], cwd=packet.ROOT, capture_output=True, text=True, timeout=180)
    assert result.returncode == 0, result.stdout+result.stderr
    assert 'scope replay passed' in result.stdout
