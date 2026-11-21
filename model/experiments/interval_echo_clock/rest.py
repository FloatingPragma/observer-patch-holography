"""Verify full material custody and the parked clock's retained offset on 480-opportunity histories."""
from argparse import ArgumentParser
from hashlib import sha256
import json
from pathlib import Path

from model.candidate_laws.interval_echo.clock import material_replay
from model.executable_model import federation

HERE = Path(__file__).resolve().parent


def count_omissions(events):
    return sum(e['primary'] - e['tick'] for e in events)


def run():
    cases = []
    for seed in (31, 47):
        first_leg_pair = None
        for cycles, hold in ((60, 0), (1, 472)):
            trip = material_replay(federation(0), (2, 6)*cycles, age_bits=5, seed=seed, rest_after=hold)
            events = trip['local_events']
            if len(events) != 480 or sum(e['handoff'] for e in events) != 2*cycles:
                raise AssertionError('wrong physical itinerary')
            if first_leg_pair is None:
                first_leg_pair = events[:8]
            elif first_leg_pair != events[:8]:
                raise AssertionError('the matched first trip diverged')
            if hold and any(e['handoff'] or e['tick'] != e['primary'] for e in events[8:]):
                raise AssertionError('parked patch lost or suppressed an ordinary repair')
            posted = sum(e['tick'] for e in events)
            if posted != trip['home_retained_traveller'] or posted != trip['reunion']['traveller']:
                raise AssertionError('the two readers disagree at reunion')
            cases.append(dict(seed=seed, cycles=cycles, rest_opportunities=hold,
                              gap_after_first_trip=count_omissions(events[:8]), gap_at_final_read=count_omissions(events),
                              reunion=trip['reunion'], intercell_pi_contacts=trip['gates']['intercell_pi_contacts']))
    return dict(protocol_sha256=sha256((HERE/'REST_PROTOCOL.md').read_bytes()).hexdigest(), material=cases)


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--record', action='store_true')
    args = parser.parse_args()
    receipt = run()
    path = HERE/'rest.json'
    if args.record:
        with path.open('x', encoding='utf-8') as f:
            f.write(json.dumps(receipt, indent=2)+'\n')
    elif receipt != json.loads(path.read_text()):
        raise SystemExit('material histories disagree with the record')
    print(json.dumps(receipt, indent=2))
