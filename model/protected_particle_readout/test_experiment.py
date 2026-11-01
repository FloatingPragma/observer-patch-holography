from fractions import Fraction as F

import numpy as np

from ..cube_eta_protection import repair
from . import experiment


def test_actual_full_field_detector_is_attached_to_both_certified_quadratures(readout):
    value, arrays = readout
    result = value['executed_native_readout']
    assert len(arrays) == 9
    for setting in result['pointer_settings']:
        low, high = map(F, setting['rigorous_actual_source_pointer_zero_interval'])
        assert float(low) <= setting['pointer_probabilities'][0] <= float(high)
        assert setting['complete_output_error_diagnostic'] < 3e-6
        sample = setting['simulated_counts']
        assert sample['pointer_zero']+sample['all_other_complete_outcomes'] == 1 << 24
        low, high = map(F, sample['source_plus_sampling_interval'])
        assert low <= F(sample['pointer_zero'], 1 << 24) <= high
    # Even after the declared sampling and host errors, the imaginary fringe
    # differs from the exact matched contact-free prediction.
    low, _ = map(F, result['pointer_settings'][1]['simulated_counts']['source_plus_sampling_interval'])
    assert low > F(1, 2)


def test_every_native_repair_output_is_reconstructible_with_all_zero_slots(readout):
    value, _ = readout
    rows = value['executed_native_readout']['all130_repair_records']
    index = 0
    for law in ('C', 'E'):
        for seam_index, _ in enumerate(repair.SEAMS):
            labels = [('native', coin, raw) for coin in (0, 1) for raw in range(4 if law == 'C' else 8)]+[('outside',)]
            for label in labels:
                vector = experiment.repair_output(law, seam_index, label)
                assert vector.shape == (131072, 1)
                probability = float(np.vdot(vector, vector).real)
                np.testing.assert_allclose(probability, rows[index]['probability'], atol=1e-42, rtol=2e-10)
                index += 1
    assert index == 130


def test_native_zero_time_and_exact_contact_free_controls(readout):
    result = readout[0]['executed_native_readout']
    zero = result['executed_zero_field_duration_control']
    assert abs(zero[0]['all_six_pointer_probabilities'][0]-1) < 1e-11
    assert abs(zero[1]['all_six_pointer_probabilities'][0]-.5) < 1e-11
    free = result['exact_contact_free_control']
    assert free['exact_Gamma_H0_Gamma_equals_minus_H0']
    assert free['free_pi_over_two_pointer_zero_probability'] == '1/2'
    assert free['additional_native_U0_word_executed'] is False
    assert abs(free['free_characteristic_at_time_one_diagnostic'][1]) < 1e-14
