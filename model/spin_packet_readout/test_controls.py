from .controls import exact_control


def test_complete_interference_reference_and_clock_countercontrols():
    result = exact_control()
    assert result['all_four_Kraus_effects_sum_to_identity']
    assert result['reference_marginal_exactly_preserved']
    assert result['maximally_entangled_input_output_probabilities'] == ['1/4']*4
    assert result['one_sided_evolution_countercontrol_rejected']
    assert result['opposite_filter_phase_countercontrol_rejected']
    assert result['subsequent_two_bin_probe_all_eight_effects_complete']


def test_regional_input_control_is_not_unknown_reference_control():
    assert exact_control()['same_marginal_joint_reference_trace_distance'] == '2'
