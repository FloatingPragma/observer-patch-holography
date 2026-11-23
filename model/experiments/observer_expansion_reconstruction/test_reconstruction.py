import json

import pytest

from . import measure


@pytest.fixture(scope='module')
def receipt():
    return measure.run()


def test_registered_events_and_replay(receipt):
    assert receipt == json.loads((measure.HERE / 'results.json').read_text())
    for comparison in receipt['layered_comparison']:
        aperture = comparison['aperture']
        assert aperture['raw_writes'] == [aperture['aperture_size']] * 4


def test_acceleration_and_independent_clock_controls(receipt):
    truth = {'dust': 0.5, 'coasting': 0.0, 'de_sitter': -1.0}
    for comparison in receipt['layered_comparison']:
        for panel in comparison['panels']:
            records = panel['histories']
            for name, target in truth.items():
                reader = records[name]['reader']
                assert all(abs(x - target) <= 0.12 for x in reader['deceleration_parameters'])
                assert reader['clock_intervals_accepted'] == [True] * 3
            alias = records['dust_with_varying_density']
            assert alias['marker_totals'] == records['de_sitter']['marker_totals']
            assert alias['proper_clock_intervals'] == records['dust']['proper_clock_intervals']
            assert not all(alias['reader']['clock_intervals_accepted'])
            assert records['exact_static_control']['deceleration_parameters'] == [None, None]


def test_publication_counts_reuse_the_clock_trigger(receipt):
    assert [entry['cells'] for entry in receipt['interval_echo_clocks']] == [20, 80, 320]
    for entry in receipt['interval_echo_clocks']:
        for owner in ('home', 'traveller'):
            read = entry[owner]
            assert read['count_mismatches'] == read['trigger_mismatches'] == 0
            assert 0 < read['positives'] == read['reunion_ticks'] < entry['services']


def test_first_interval_forecasts_held_expansion(receipt):
    from fractions import Fraction

    for comparison in receipt['layered_comparison']:
        for trial in comparison['amplitude_trials']:
            forecast = trial['forecast']
            if trial['preparation'] == 'dust':
                assert not forecast['passes_constant_rate_forecast']
                assert forecast['de_sitter_lambda_in_clock_units'] is None
            else:
                truth = float(Fraction(trial['prepared_initial_rate']))
                assert abs(forecast['hubble_rate'] / truth - 1) < 0.02
                assert forecast['passes_constant_rate_forecast']
                assert abs(forecast['de_sitter_lambda_in_clock_units'] / (3 * truth**2) - 1) < 0.041
