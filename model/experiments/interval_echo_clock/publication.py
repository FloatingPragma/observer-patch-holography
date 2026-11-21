"""Check the retained local tick and its two native clock readers on the frozen material journeys."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import clock
from model.executable_model import federation
from model.executable_model.carrier import seam_digest

HERE = Path(__file__).resolve().parent
PROTOCOL = HERE / 'PUBLICATION_PROTOCOL.md'
SCHEDULE = (8, 1, 2, 4, 16, 2)


def independent_probe():
    """Compute the endpoint from 99-state exposure columns, independently of the transported state."""
    law = clock.law
    h0, h1 = law.local_term_parts()
    H = h0 + clock.KAPPA * h1
    birth = np.eye(law.DIM)[:, law.SRC.index((2, 0))]
    pulse = expm(-1j * np.pi / 16 * H)
    shifted = expm(-1j * np.pi / 2 * H)
    row = law.SRC.index((2, 0))
    tables = []
    state = birth
    for _ in range(sum(SCHEDULE) + 1):
        tables.append([abs(state[row])**2, abs((shifted @ state)[row])**2])
        state = pulse @ state
    return np.array(tables)


def run():
    probability = independent_probe()
    specimens = []
    transitions = Counter()
    residuals = []
    for level in range(3):
        car = federation(level)
        for width in (3, 4):
            for seed in (31, 2026112301):
                records = [clock.material_replay(car, SCHEDULE, age_bits=width, seed=seed, analyzer=a)
                           for a in (0, 1)]
                first, second = records
                events = first['local_events']
                if events != second['local_events']:
                    raise AssertionError("terminal analyzer changed earlier records")
                moving = [e['tick'] for e in events]
                home = [e['home_tick'] for e in events]
                for i, (event, m, h) in enumerate(zip(events, moving, home, strict=True)):
                    if m != event['oscillator_tick'] or h != event['home_oscillator_tick'] or h != event['home_primary']:
                        raise ArithmeticError('a mechanical and coherent listener disagree on publication')
                    if i and moving[i-1] and not m:
                        reason = 'rejected_positive_after_one' if event['primary'] else 'quiet_after_one'
                        transitions[reason] += 1
                    if i in (7, sum(SCHEDULE[:-2]) + 15) and event['handoff'] and event['primary']:
                        transitions[f'positive_boundary_width_{width}_at_{i+1}'] += 1
                probe_readings = {}
                for name, bits in (('traveller', moving), ('home', home)):
                    if np.cumsum(bits).tolist() != [e[name] for e in events]:
                        raise ArithmeticError('five-position hand left its public events')
                    entries = [row['oscillators'][name] for row in records]
                    for entry in entries:
                        residuals.extend(abs(probability[sum(bits)] - entry['first_rail_probabilities']))
                        if entry['owners'] != [entry['bit']]*2:
                            raise ArithmeticError('a terminal shared bit was not agreed')
                    probe_readings[name] = {'public_count': sum(bits), 'rails': entries[0]['rails'],
                                            'probabilities': entries[0]['first_rail_probabilities'],
                                            'published_bits_at_settings': [e['bit'] for e in entries]}
                specimens.append({'cells': car.cells, 'seams_sha256': seam_digest(car), 'width': width, 'seed': seed,
                                  'physical_event_digits': first['elapsed_event_registers'],
                                  'primary': [e['primary'] for e in events],
                                  'moving_public': moving, 'home_public': home,
                                  'handoffs_at_services': [i + 1 for i, e in enumerate(events) if e['handoff']],
                                  'endpoint': probe_readings, 'gates': first['gates']})
    if max(residuals) > 1e-10 or min(transitions[k] for k in ('quiet_after_one', 'rejected_positive_after_one')) < 1:
        raise ArithmeticError('endpoint residual or frozen clearing control failed')
    return {'protocol_sha256': hashlib.sha256(PROTOCOL.read_bytes()).hexdigest(), 'schedule': list(SCHEDULE),
            'exact_pointer_residual_below': 1e-10, 'observed_transition_counts': dict(transitions),
            'specimens': specimens}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, help='write a new receipt')
    args = p.parse_args()
    result = run()
    if args.output:
        with args.output.open('x') as f:
            f.write(json.dumps(result, separators=(',', ':'), allow_nan=False) + '\n')
    elif result != json.loads((HERE / 'publication.json').read_text()):
        raise SystemExit('public clock records differ from receipt')
    print(result['observed_transition_counts'])
