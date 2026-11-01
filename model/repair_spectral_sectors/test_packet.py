"""Cold-source custody and hostile complete-response/accounting mutations."""
from copy import deepcopy

import pytest

from tools import strict_json
from . import packet


@pytest.fixture(scope='module')
def replay():
    return packet.report()


def test_retained_packet_cold_replays_every_result(replay):
    retained = strict_json.load(packet.OUTPUT)
    packet.compare(retained, replay)
    packet.validate_accounting(retained)
    assert replay['exact_active_cyclicity']['specialized_rank'] == 90
    assert replay['exact_dark_eta_spectrum']['multiplicity'] == 15


@pytest.mark.parametrize('key,value', [('integrated_gate_accepted', True),
                                     ('physical_comparison_count_added', 1),
                                     ('primitive_or_repair_law_changed', True),
                                     ('native_pulse_source_executed', True),
                                     ('all_Fock_common_algebra_classified', True),
                                     ('controller_work_or_heat_identified', True)])
def test_execution_physical_or_scope_promotions_reject(replay, key, value):
    changed = deepcopy(replay); changed[key] = value
    with pytest.raises(ValueError):
        packet.compare(changed, replay)


def test_exact_witness_and_all_retained_modes_are_required(replay):
    changed = deepcopy(replay)
    changed['exact_active_cyclicity']['pivot_witness'][0][2] += 1
    with pytest.raises(ValueError):
        packet.compare(changed, replay)
    changed = deepcopy(replay)
    changed['exact_dark_eta_spectrum']['other_energy_multiplicities'].pop()
    with pytest.raises(ValueError):
        packet.compare(changed, replay)
    changed = deepcopy(replay)
    changed['complete_dark_spectral_response']['complete_stopping_history_outputs'].pop()
    with pytest.raises(ValueError):
        packet.validate_accounting(changed)


def test_tiny_leaf_errors_cannot_destroy_collective_probability_mass(replay):
    changed = deepcopy(replay)
    for row in changed['complete_active_energy_response']['complete_energy_transition_matrix']:
        for j in range(len(row)):
            row[j] += 9e-11
    packet.compare(changed, replay)  # Below the separate leaf allowance.
    with pytest.raises(ValueError):
        packet.validate_accounting(changed)
    changed = deepcopy(replay)
    changed['complete_active_energy_response']['mean_K_transfer_per_energy_input'][0] += 1e-5
    with pytest.raises(ValueError):
        packet.validate_accounting(changed)


@pytest.mark.parametrize('bad', [-1e-20, float('nan'), float('inf'), True, 0, '0'])
def test_bad_numerical_evidence_rejects(bad):
    with pytest.raises(ValueError):
        packet.compare({'probability': bad}, {'probability': 0.})


def test_negative_energies_and_transfers_allowed_but_admissions_fixed():
    packet.compare({'mean_K_transfer_per_energy_input': [-.1]}, {'mean_K_transfer_per_energy_input': [-.1]})
    with pytest.raises(ValueError):
        packet.compare({'numerical_admission': 3.1e-10}, {'numerical_admission': 3e-10})


def test_changed_source_snapshot_rejects_before_worker(monkeypatch):
    changed = packet.source_bindings()
    changed['model/repair_spectral_sectors/algebra.py'] = '0'*64
    monkeypatch.setattr(packet, 'source_bindings', lambda: changed)
    with pytest.raises(ValueError, match='fresh interpreter'):
        packet.report()


def test_transitive_source_manifest_includes_actual_laws_and_exact_proof():
    bindings = packet.source_bindings()
    for name in ('model/repair_spectral_sectors/algebra.py', 'model/repair_spectral_sectors/response.py',
                 'model/repair_spectral_sectors/test_algebra.py', 'model/repair_spectral_sectors/test_packet.py',
                 'model/repair_spectral_sectors/contract.json', 'model/repair_spectral_sectors/README.md',
                 'model/common_spin_protection/results.json', 'model/charged_defect/results.json',
                 'model/charged_defect/dictionary.py', 'model/common_spin_protection/spin.py',
                 'model/field_symmetry/symmetry.py', 'model/particle_transport/algebra.py',
                 'paper/tex_fragments/REPAIR_SPECTRAL_SECTORS.tex', 'model/requirements.txt'):
        assert name in bindings
    assert 'model/repair_spectral_sectors/results.json' not in bindings


def test_duplicate_keys_and_extra_claims_reject(replay):
    with pytest.raises(ValueError):
        strict_json.loads('{"accepted": false, "accepted": true}')
    with pytest.raises(ValueError):
        packet.compare({**replay, 'extra': True}, replay)
