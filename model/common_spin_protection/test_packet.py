"""Cold-source custody and hostile retained-spectrum/history mutations."""
from copy import deepcopy

import pytest

from tools import strict_json
from . import packet


@pytest.fixture(scope='module')
def replay():
    return packet.report()


def test_retained_packet_cold_replays_every_current_result(replay):
    retained = strict_json.load(packet.OUTPUT)
    packet.compare(retained, replay)
    packet.validate_joint_probabilities(retained)
    assert replay['same_full_field_and_repair_checks']['complete_repair_outcomes_checked'] == 156
    assert replay['complete_finite_spectrum']['absolute_physical_masses_predicted'] == 0


@pytest.mark.parametrize('key,value', [('integrated_gate_accepted', True),
                                     ('physical_comparison_count_added', 1),
                                     ('primitive_or_repair_law_changed', True),
                                     ('native_pulse_source_executed', True)])
def test_physical_or_execution_promotions_reject(replay, key, value):
    changed = deepcopy(replay); changed[key] = value
    with pytest.raises(ValueError):
        packet.compare(changed, replay)


def test_spectrum_counts_complete_outputs_and_alignment_are_not_optional(replay):
    changed = deepcopy(replay)
    changed['complete_finite_spectrum']['full_charge_spin_spectrum'][0]['energies'].append(0.)
    with pytest.raises(ValueError):
        packet.compare(changed, replay)
    changed = deepcopy(replay)
    changed['complete_joint_history']['complete_joint_outcomes'].pop()
    with pytest.raises(ValueError):
        packet.validate_joint_probabilities(changed)
    changed = deepcopy(replay)
    changed['general_block_alignment']['exact_parity_cases'][1]['aligned'] = True
    with pytest.raises(ValueError):
        packet.compare(changed, replay)


def test_small_per_leaf_errors_cannot_corrupt_complete_probability_accounting(replay):
    changed = deepcopy(replay)
    for row in changed['complete_joint_history']['complete_joint_outcomes']:
        row['probability'] += 9e-11
    # Every separate delta is smaller than the replay tolerance, but the
    # joint probability mass is not.  Semantic completeness must still fail.
    with pytest.raises(ValueError):
        packet.validate_joint_probabilities(changed)


@pytest.mark.parametrize('bad', [-1e-20, float('nan'), float('inf'), True, 0, '0'])
def test_nonfinite_negative_or_wrong_typed_numerical_evidence_rejects(bad):
    with pytest.raises(ValueError):
        packet.compare({'probability': bad}, {'probability': 0.})


def test_negative_energies_are_admitted_but_negative_gaps_and_weights_are_not():
    packet.compare({'energies': [-3.], 'ground_energy': -3.}, {'energies': [-3.], 'ground_energy': -3.})
    for name in ('hamiltonian_model_gaps', 'all_highest_weight_residues'):
        with pytest.raises(ValueError):
            packet.compare({name: [-1e-20]}, {name: [0.]})
    with pytest.raises(ValueError):
        packet.compare({'numerical_admission': 2.1e-10}, {'numerical_admission': 2e-10})


def test_changed_source_snapshot_rejects_before_worker(monkeypatch):
    changed = packet.source_bindings()
    changed['model/common_spin_protection/spin.py'] = '0'*64
    monkeypatch.setattr(packet, 'source_bindings', lambda: changed)
    with pytest.raises(ValueError, match='fresh interpreter'):
        packet.report()


def test_source_closure_binds_consumed_laws_exact_proof_tests_and_parent_receipts():
    bindings = packet.source_bindings()
    for name in ('model/common_spin_protection/spin.py', 'model/common_spin_protection/history.py',
                 'model/common_spin_protection/spectrum.py', 'model/common_spin_protection/test_spin.py',
                 'model/common_spin_protection/test_spectrum.py', 'model/common_spin_protection/test_history.py',
                 'model/common_spin_protection/test_packet.py', 'model/common_spin_protection/contract.json',
                 'model/charged_defect/results.json', 'model/charged_defect/dictionary.py',
                 'model/field_symmetry/results.json', 'model/field_symmetry/symmetry.py',
                 'model/field_stress/energy.py', 'model/particle_transport/algebra.py',
                 'paper/tex_fragments/COMMON_FIELD_SPIN_PROTECTION.tex', 'model/requirements.txt'):
        assert name in bindings
    assert 'model/common_spin_protection/results.json' not in bindings


def test_duplicate_keys_and_extra_claims_reject(replay):
    with pytest.raises(ValueError):
        strict_json.loads('{"accepted": false, "accepted": true}')
    with pytest.raises(ValueError):
        packet.compare({**replay, 'extra': True}, replay)
