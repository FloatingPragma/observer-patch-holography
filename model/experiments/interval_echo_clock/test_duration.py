import json


from .duration import HERE, history, prediction


def test_persistent_rate_and_transient_offset_depend_on_the_trajectory():
    _, _, once = prediction(*history('repeated_motion', 1))
    for cycles in (4, 16, 60):
        _, _, ongoing = prediction(*history('repeated_motion', cycles))
        _, _, parked = prediction(*history('one_trip_then_rest', cycles))
        assert ongoing == cycles*once
        assert parked == once
    raw, accepted, gap = prediction(*history('at_rest', 60))
    assert abs(raw-accepted) < 1e-14
    assert gap == 0


def test_unresolved_terminal_read_can_coexist_with_a_measured_count_difference():
    rows = json.loads((HERE/'duration.json').read_text())['histories']
    panel = next(r for r in rows if r['history']=='slow_repeated_motion' and r['opportunities']==434)['panels'][1]
    assert panel['estimated_ratio_gap'] > 4*panel['paired_standard_error']
    assert not panel['terminal_split_resolved']



def test_parked_patch_retains_each_local_mark_and_reunites(capsys):
    from model.__main__ import main

    assert main(['run', '--law', 'interval_echo', '--clock-trip', '--rest-after', '12', '2', '6']) == 0
    trip = json.loads(capsys.readouterr().out)
    parked = trip['local_events'][8:]
    assert trip['rest_after'] == len(parked) == 12
    assert all(not event['handoff'] and event['primary'] == event['tick'] for event in parked)
    assert trip['reunion']['traveller'] == sum(event['tick'] for event in trip['local_events'])
