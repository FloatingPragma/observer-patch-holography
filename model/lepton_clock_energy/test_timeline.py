from fractions import Fraction as F
from itertools import product

import numpy as np

from . import timeline as t


def test_literal_source_and_reference_isometry():
    blocks = t.source_blocks()
    for kind, original in [('weak', t.lepton.weak_instrument(float(t.KICK))),
                           ('mass', t.lepton.mass_instrument(float(t.THETA), float(t.KICK)))]:
        for label in t.lepton.LABELS:
            assert np.linalg.norm(blocks['lift'] @ blocks['bank'][kind][label]
                                  - original[label].reshape(3600, 3)) < 3e-12
    # Entangled family/reference input, with every quantum archive kept.
    state = np.array([[1, 2j], [3-1j, -2], [1j, 4]], complex)
    state /= np.linalg.norm(state)
    reference_output = np.zeros((2, 2), complex)
    for returns in product(t.lepton.LABELS, repeat=2):
        archive = t.archive_columns(returns)
        vector = archive[:, :, None]*state[None, :, :]
        reference_output += np.einsum('akr,aks->rs', vector, vector.conj())
    assert np.linalg.norm(reference_output-state.T@state.conj()) < 3e-12


def test_exact_native_full_state_and_complete_stopping_law():
    for law in ('C', 'E'):
        rows, error = t.stopped_clock(law, 3, F(1, 2))
        assert len(rows) == 52
        assert sum(F(r['probability_exact']) for r in rows) == 1
        assert error < 1e-12
        assert any('idle' in r['word'] and F(r['probability_exact']) > 0 for r in rows)
    result = t.joint_run()
    assert result['complete_record_count'] == 1272
    assert F(result['timeout_probability_exact']) == F(15806, 16875)
    assert any(r['probability'] == 0 for r in result['complete_joint_records'])


def test_waiting_and_correlated_reset_are_observable():
    event = [t.joint_run(initial_idles=i) for i in (0, 1)]
    assert event[0]['detector_probabilities'] == event[1]['detector_probabilities']
    assert F(event[1]['mean_controller_slots_exact'])-F(event[0]['mean_controller_slots_exact']) == 1
    slot = [t.joint_run(initial_idles=i, work_policy='slot') for i in (0, 1)]
    assert np.linalg.norm(np.array(slot[0]['detector_probabilities'])-slot[1]['detector_probabilities']) > 1e-5
    assert t.no_reset_control()['return_total_variation'] > .0023


def test_generator_obstruction_does_not_preclude_timeline():
    result = t.generator_and_scale()
    for row in result['all_six_affine_dictionary_obstructions']:
        lo, hi = map(F, row['second_difference_interval'])
        assert lo > 0 or hi < 0
    assert result['distinct_generators_do_not_prevent_a_shared_timeline']
    assert not result['absolute_GeV_or_SI_scale_derived']
    assert max(r['phase_identity_error'] for r in result['dimensionful_rescaling_controls']) < 5e-16
