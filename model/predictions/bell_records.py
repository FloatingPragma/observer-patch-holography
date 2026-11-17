"""A classical recorded-branch boundary and a published Bell-data replay.

The Bell bound is reused mathematics, not a new OPH prediction. Its application
requires factorized local responses, setting independence and the declared
event-ready eligibility contract. The complete quantum observer is outside
that null. No native clock is identified with experimental spacelike separation.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime
from fractions import Fraction as F
import hashlib
from itertools import product
import json
from math import comb
from pathlib import Path

from model.observer import native
from model.observer.scalar_response.experiment import compare_receipt
from tools import strict_json
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / 'model/predictions/data/hensen2015.txt'
DATA_SHA256 = '5a6b1853554ad508d35c795415d64922ed9e5bc0469f6052ad222486ad3566ff'
DATA_DOI = 'https://doi.org/10.4121/uuid:6e19e9b2-4a2d-40b5-8dd3-a660bf3c0a31'
SOURCE_PATHS = ('model/observer/native.py', 'model/predictions/bell_records.py')
# Published calibration allowance; it is not inferred from these Bell records.
TAU = F(27, 2500000)
SETTINGS = tuple(product((0, 1), repeat=2))
OUTCOMES = tuple(product((1, -1), repeat=2))


def local_strategies():
    """All extremal local binary response tables; shared randomness mixes them."""
    result = []
    for a0, a1, b0, b1 in product((-1, 1), repeat=4):
        responses = ((a0, a1), (b0, b1))
        wins = tuple(int((-1)**(a*b)*responses[0][a]*responses[1][b] == 1)
                     for a, b in SETTINGS)
        result.append((responses, wins))
    return tuple(result)


def native_local_example():
    """Execute local native repairs on two disjoint carrier triangles.

    Settings select seams (0,1)/(0,5) and (3,8)/(3,9). Initial local
    occupations are disclosed and signed currents update the two counters.
    No inter-region repair occurs between the setting and readout. This
    example checks the native-to-classical mapping, not the laboratory apparatus.
    """
    regions = (((0, 1), (0, 5)), ((3, 8), (3, 9)))
    counts = Counter()
    for initial in native.states(binary=True):
        for a, b in SETTINGS:
            for ca, cb in SETTINGS:
                left, ra = native.native(initial, regions[0][a], ca)
                final, rb = native.native(left, regions[1][b], cb)
                reverse, _ = native.native(initial, regions[1][b], cb)
                reverse, _ = native.native(reverse, regions[0][a], ca)
                x, y = initial[0]-ra[1], initial[3]-rb[1]
                if reverse != final or (x, y) != (final[0], final[3]):
                    raise ArithmeticError('local counters or disjoint source chronology changed')
                counts[a, b, 2*x-1, 2*y-1] += 1
    correlations = tuple(sum(x*y*counts[a, b, x, y] for x, y in OUTCOMES)/F(264)
                         for a, b in SETTINGS)
    score = sum((-1)**(a*b)*c for (a, b), c in zip(SETTINGS, correlations, strict=True))
    return {'source_states': 66, 'complete_source_histories': sum(counts.values()),
            'settings': [list(row) for row in SETTINGS],
            'outcomes': [list(row) for row in OUTCOMES],
            'counts': [[counts[a, b, x, y] for x, y in OUTCOMES] for a, b in SETTINGS],
            'CHSH': str(score), 'winning_probability': str(F(1, 2)+score/8),
            'all_current_readbacks_and_disjoint_orders_agree': True,
            'C_and_E_share_diagonal_law': True,
            'laboratory_apparatus_identified': False}


def parse_data(raw):
    """Read every deposited row; event selection is a separate operation."""
    if type(raw) is not bytes:
        raise ValueError('original dataset bytes required')
    rows = []
    for line in raw.decode('ascii').splitlines():
        parts = line.split(',')
        if len(parts) != 17:
            raise ValueError('each dataset row needs a timestamp and sixteen integers')
        datetime.strptime(parts[0], '%Y-%m-%d %H:%M:%S.%f')
        if any(not p.isascii() or not p.isdigit() for p in parts[1:]):
            raise ValueError('nonnegative decimal integer source columns required')
        row = tuple(map(int, parts[1:]))
        if any(row[i] not in (0, 1) for i in (3, 5, 6, 7)):
            raise ValueError('original binary photon-channel and setting values required')
        rows.append(row)
    if not rows:
        raise ValueError('empty experimental dataset')
    return tuple(rows)


def eligible(row):
    """The authors' fixed eligibility rule; no setting or spin outcome enters."""
    starts = (5426350, 5425700)
    click1, channel1, click2, channel2 = row[2:6]
    return (starts[channel1] <= click1 < starts[channel1]+52450
            and starts[channel2]+250000 <= click2 < starts[channel2]+302450
            and channel1 != channel2
            and all(row[i] == 0 or row[i] > 250 for i in (14, 15))
            and row[12] == row[13] == 0)


def score_data(rows):
    counts = Counter()
    indices, wins = [], 0
    for index, row in enumerate(rows):
        if not eligible(row):
            continue
        a, b = row[6:8]
        x, y = (1 if 10620 < row[i] <= 14320 else -1 for i in (10, 11))
        counts[a, b, x, y] += 1
        wins += int((-1)**(a*b)*x*y == 1)
        indices.append(index)
    totals = tuple(sum(counts[a, b, x, y] for x, y in OUTCOMES) for a, b in SETTINGS)
    if not all(totals):
        raise ValueError('every setting needs at least one eligible record for the CHSH estimate')
    correlations = tuple(F(sum(x*y*counts[a, b, x, y] for x, y in OUTCOMES), n)
                         for (a, b), n in zip(SETTINGS, totals, strict=True))
    chsh = sum((-1)**(a*b)*v for (a, b), v in zip(SETTINGS, correlations, strict=True))
    return {'deposited_candidate_rows': len(rows), 'eligible_trials': len(indices),
            'ineligible_candidate_rows': len(rows)-len(indices), 'wins': wins,
            'eligible_zero_based_rows': indices,
            'counts_by_setting_and_outcome': [[counts[a, b, x, y] for x, y in OUTCOMES]
                                             for a, b in SETTINGS],
            'setting_totals': list(totals), 'correlations': list(map(str, correlations)),
            'CHSH_estimate': str(chsh), 'game_winning_fraction': str(F(wins, len(indices)))}


def tail(n, wins, probability):
    """Exact fixed-n binomial domination for conditional win bounds, with memory."""
    if (type(n) is not int or not 1 <= n <= 1000 or type(wins) is not int
            or not 0 <= wins <= n or type(probability) is not F or not 0 <= probability <= 1):
        raise ValueError('bounded trial counts and an exact probability required')
    a, d = probability.numerator, probability.denominator
    return F(sum(comb(n, k)*a**k*(d-a)**(n-k) for k in range(wins, n+1)), d**n)


def outward(value, bits=80):
    grid = 1 << bits
    numerator = value.numerator*grid
    return [str(F(numerator//value.denominator, grid)),
            str(F(-(-numerator//value.denominator), grid))]


def run():
    raw = DATA.read_bytes()
    if hashlib.sha256(raw).hexdigest() != DATA_SHA256:
        raise ValueError('the deposited data bytes differ from the pinned source')
    rows = parse_data(raw)
    result = score_data(rows)
    strategies = local_strategies()
    q = F(3, 4)+3*(TAU+TAU**2)
    p = tail(result['eligible_trials'], result['wins'], q)
    if max(sum(wins) for _, wins in strategies) != 3:
        raise ArithmeticError('the complete deterministic local response bound changed')
    return {'schema': 'oph.predictions.bell_records.v1',
        'comparison': 'retrospective published Bell null replay, not a prospective OPH prediction',
        'data_source': {'doi': DATA_DOI, 'paper': 'https://doi.org/10.1038/nature15759',
            'archive_url': 'https://ndownloader.figshare.com/files/24056582',
            'archive_md5': '342f29f8288c46575818acd2acebd535',
            'original_file': 'bell_open_data.txt', 'retained_sha256': DATA_SHA256,
            'reuse_terms': 'https://doi.org/10.4121/resource:terms_of_use'},
        'classical_record_boundary': {'deterministic_local_strategies': len(strategies),
            'win_counts': dict(sorted(Counter(str(sum(w)) for _, w in strategies).items())),
            'uniform_setting_CHSH_upper': 2, 'uniform_setting_win_upper': '3/4',
            'native_example': native_local_example()},
        'published_selection_analysis': result,
        'conditional_null_test': {'published_setting_allowance_tau': str(TAU),
            'conditional_win_upper': str(q), 'fixed_n_tail_interval': outward(p),
            'uniform_setting_tail_interval': outward(tail(result['eligible_trials'], result['wins'], F(3, 4))),
            'tail_below_one_twentieth': p < F(1, 20),
            'device_memory_allowed': True, 'IID_device_trials_required': False,
            'optional_stopping_certificate_supplied': False,
            'setting_independence_locality_and_eligibility_authenticated_by_this_code': False},
        'scope': {'diagonal_local_factorized_response_extension_excluded_under_stated_null': p < F(1, 20),
            'whole_OPH_excluded': False, 'quantum_C_or_E_selected': False,
            'new_native_empirical_confirmation': False, 'lab_clock_calibrated': False,
            'native_repairs_identified_with_laboratory_spacelike_separation': False},
        'source_bindings': bound_pins({name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                            for name in SOURCE_PATHS})}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--check', type=Path)
    modes.add_argument('--write', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare_receipt(strict_json.load(args.check), result)
        print('Published Bell records and conditional classical-branch boundary replayed.')
    elif args.write:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write('\n')
    else:
        print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


if __name__ == '__main__':
    main()
