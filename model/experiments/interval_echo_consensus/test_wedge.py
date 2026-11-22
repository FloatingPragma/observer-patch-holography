import json
from fractions import Fraction

from . import __main__ as wedge
from . import phase
from . import collar_reference


def test_two_unit_channels_and_balanced_input_control():
    result = wedge.run()
    saved = json.loads((wedge.HERE / 'results.json').read_text())
    assert json.loads(json.dumps(result)) == saved
    assert [row['middle_double_order_trace_distance'] for row in result['priority_rows']] == [1, 0, 0]
    assert [r['initial'] for r in result['legacy_ready_on_balanced_blank']] == [[1, 1, 0], [0, 1, 1]]
    assert all(r['opposite_phases_with_transcript_distance'] == 0 for r in result['priority_rows'])


def test_addressed_window_preserves_within_winner_phase():
    result = phase.run()
    assert result == json.loads((wedge.HERE / 'phase.json').read_text())
    for row in result['priority_rows']:
        assert row['priority_probe_rail_probabilities'] == [0, 1]
        assert row['both_probes_rail_probabilities'] == [.5, .5]


def test_local_reference_forecasts_held_magnitude():
    observed = collar_reference.measure()
    assert observed["public_records"]["magnitude_two_ratio"] == Fraction(132651, 20000)
    assert observed["equal_parity_splitting_matches"]
    assert observed["all_walk_lengths_ratio_floor"] == Fraction(3456, 625)
    assert json.loads(json.dumps(observed, default=str)) == json.loads(
        (wedge.HERE / "collar_reference.json").read_text()
    )
