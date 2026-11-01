from fractions import Fraction as F

from . import experiment as e


def test_sixteen_certified_settings_resolve_all_eight_weights_from_counts():
    result = e.execution()
    assert len(result['settings']) == 8
    assert result['total_simulated_reset_trials'] == 77309411328
    assert result['weight_reconstruction']['all_eight_weights_strictly_positive']
    assert result['weight_reconstruction']['source_weights_used_in_reconstruction'] is False
    assert all(row['independent_exact_weight_contained'] for row in result['independent_residue_comparison'])
    assert F(result['weight_reconstruction']['weight_error_upper']) < F(3, 10000)
    for time in result['settings']:
        for row in time['pointer_settings']:
            assert F(row['total_systematic_probability_bound']) < F(1, 65536)
            assert row['pointer_zero_count']+row['all_other_complete_outcomes_count'] == e.TRIALS


def test_zero_time_word_and_scope_are_retained():
    result = e.execution()
    rows = result['settings'][0]['pointer_settings']
    assert rows[0]['sampling_probability_binary64'] == 1.
    assert rows[1]['sampling_probability_binary64'] == .5
    assert result['actual_raw_native_long_word_executed'] is False
    assert result['unknown_frequency_support_inferred'] is False
    assert result['physical_particle_species_assigned'] == 0
