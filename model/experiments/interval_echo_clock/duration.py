"""Long-window ensemble elapsed counts and resolution of the native two-quadrature terminal reader."""
import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from model.candidate_laws.interval_echo import clock, cohorts
from model.executable_model import federation
from .experiment import PULSE, joint_marks, native_probe

HERE = Path(__file__).resolve().parent
BUDGETS = (128, 8192)
SEED = 2026112357


def history(kind, cycles):
    if kind == 'repeated_motion':
        return 8*cycles, Counter({2: cycles, 6: cycles})
    if kind == 'one_trip_then_rest':
        return 8*cycles, Counter({2: 1, 6: 1})
    if kind == 'slow_repeated_motion':
        return 62*cycles, Counter({31: 2*cycles})
    if kind == 'at_rest':
        return 8*cycles, Counter()
    raise ValueError(kind)


def prediction(total, closings):
    step = np.exp(1j*PULSE)
    ordinary = (1 + step)/2
    accepted = ordinary**(total - sum(closings.values()))
    omitted = 0.
    for q, count in closings.items():
        suppression = cohorts.suppression(q, 5)
        tick = (1 - suppression)/2
        accepted *= (1 + tick*(step-1))**count
        omitted += count * suppression/2
    return ordinary**total, accepted, omitted


def pointer_probabilities(fringe):
    return np.clip(np.array([(1+fringe.real)/2, (1-fringe.imag)/2]), 0, 1)


def check_one_trip():
    support, weights = joint_marks((False, True, False, False, False, False, False, True))
    raw, accepted, omitted = prediction(*history('repeated_motion', 1))
    for col, fringe in enumerate((raw, accepted)):
        exact_pointer = native_probe()[:, support[:, col]] @ weights
        np.testing.assert_allclose(pointer_probabilities(fringe), exact_pointer, rtol=0, atol=1e-12)
    np.testing.assert_allclose(omitted, (support[:,0] - support[:,1]) @ weights, rtol=0, atol=1e-12)


def panel(total, closings, raw, accepted, n, key):
    rng = np.random.default_rng(np.random.SeedSequence([SEED, key, n]))
    # Given the itinerary, only a positive-and-suppressed closing repair contributes to the paired gap.
    missing = np.zeros(n, dtype=int)
    for q, number in closings.items():
        missing += rng.binomial(number, cohorts.suppression(q, 5)/2, size=n)
    home = rng.binomial(total, .5, size=n)
    ratio = missing.sum() / home.sum()
    se = np.std(missing-ratio*home, ddof=1) / (np.sqrt(n)*np.mean(home))
    clicks = rng.binomial(n, [pointer_probabilities(raw), pointer_probabilities(accepted)])
    freq = clicks/n
    pooled = freq.mean(axis=0)
    margin = 4*np.sqrt(2*pooled*(1-pooled)/n)
    return dict(shots_per_analyzer=n, paired_omitted_marks=int(missing.sum()), home_marks=int(home.sum()),
                estimated_ratio_gap=round(float(ratio),10), paired_standard_error=round(float(se),10),
                rail_clicks=clicks.tolist(), four_se_quadrature_margin=margin.round(10).tolist(),
                terminal_split_resolved=bool(np.any(np.abs(freq[0]-freq[1]) > margin)))


def run():
    check_one_trip()
    rows = []
    for kind in ('repeated_motion','one_trip_then_rest','slow_repeated_motion','at_rest'):
        cycles = (1, 4, 7) if kind == 'slow_repeated_motion' else (1, 4, 16, 32, 60)
        for repeats in cycles:
            total, closings = history(kind, repeats)
            raw, accepted, omitted = prediction(total, closings)
            rows.append(dict(history=kind, cycles_or_windows=repeats, opportunities=total,
                             expected_count_gap=round(omitted,12), expected_ratio_gap=round(2*omitted/total,12),
                             fringe_order=['raw','accepted'],
                             expected_fringe=[[round(z.real,12), round(z.imag,12)] for z in (raw,accepted)],
                             visibility=[round(abs(z),12) for z in (raw,accepted)],
                             panels=[panel(total, closings, raw, accepted, n, len(rows)) for n in BUDGETS]))
    material = []
    for cycles in (1,4):
        trip = clock.material_replay(federation(0), (2,6)*cycles, age_bits=5, seed=SEED)
        posted = sum(e['tick'] for e in trip['local_events'])
        if posted != trip['reunion']['traveller'] or posted != trip['home_retained_traveller']:
            raise ArithmeticError('material reunion differs from the accepted history')
        material.append(dict(cycles=cycles, observed_primary_marks=sum(e['primary'] for e in trip['local_events']),
                             reunion=trip['reunion']))
    return dict(protocol_sha256=sha256((HERE/'DURATION_PROTOCOL.md').read_bytes()).hexdigest(),
                seed=SEED, material=material, histories=rows)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record', action='store_true')
    args = parser.parse_args()
    result = run()
    path = HERE/'duration.json'
    if args.record:
        with path.open('x', encoding='utf-8') as f:
            f.write(json.dumps(result, indent=2, allow_nan=False)+'\n')
    elif result != json.loads(path.read_text()):
        raise SystemExit('duration replay differs from the receipt')
    for r in result['histories']:
        print(r['history'], r['opportunities'], r['expected_count_gap'], r['expected_ratio_gap'],
              [p['terminal_split_resolved'] for p in r['panels']], r['visibility'])
