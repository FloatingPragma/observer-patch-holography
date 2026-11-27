"""Decode every reached branch and both physical publications from counts."""

import json
from math import log, sqrt
from pathlib import Path


def material(bits):
    value = 0
    for k, rails in enumerate(bits):
        pair = tuple(rails)
        if pair not in (('10', '01'), ('01', '10')):
            raise ValueError('damaged balanced publication')
        value += (pair == ('10', '01')) << k
    return value


def check(data):
    n, k, delta = data['cohort_size'], data['total_look_cap'], data['simultaneous_failure_budget']
    radius = sqrt(2 * log(2*k / delta) / n) + 1e-9
    consumed = offered = 0
    for entry in data['panels']:
        panel = entry['panel']
        next_edge = (0, 1)
        previous = []
        signs, estimates = [], []
        for i, row in enumerate(panel['records']):
            raw = row['private_prefix_then_terminal_counts']
            counts = [sum(raw[v::3]) for v in range(3)]
            if counts != row['terminal_counts'] or sum(counts) != n or row['cohort'] != i:
                raise ValueError('cohort population differs')
            if i < 3:
                previous.append(next_edge)
                offered += n*(4-i)
            if (previous != list(map(tuple, row['word'])) or len(raw) != 3**(len(previous)+1)
                    or row['station_owner'] != (1, 3, 1, 3)[i]):
                raise ValueError('word or observing station differs')
            mu = (counts[2] - counts[0]) / n
            label = 'high' if mu > radius else 'low' if mu < -radius else None
            signs.append(label)
            estimates.append(mu)
            if i < 2:
                if label is None and i+1 != len(panel['records']):
                    raise ValueError('unresolved decision continued')
                next_edge = (((0, 5) if label == 'low' else (0, 7)) if i == 0 else
                             ((0, 1) if label == 'low' else (0, 5)))
        if len(signs) == 3 or not signs:
            raise ValueError('incomplete terminal comparison')
        final = None
        if len(signs) == 4:
            lo = max(-1., max(estimates[-2:]) - radius)
            hi = min(1., min(estimates[-2:]) + radius)
            if signs[-1] is not None and signs[-1] == signs[-2] and lo <= hi:
                final = [lo, hi]
        label = signs[-1] if final else None
        if panel['public']['sign'] != label or panel['final_interval'] != final:
            raise ValueError('publication disagrees with the two recorded intervals')
        for flags, words in zip(panel['public']['flags_at_both_owners'], panel['public']['score_words_at_both_owners'],
                                strict=True):
            if [material([pair]) for pair in flags] != [int(label == 'low'), int(label == 'high')]:
                raise ValueError('damaged material sign')
            expected = [int(n + n*z) for z in estimates[-2:]] if label else [0, 0]
            if [material(word) for word in words] != expected:
                raise ValueError('damaged terminal scores')
        consumed += n*len(signs)
        if panel['survivor_copies'] != n*(4-len(signs)):
            raise ValueError('incorrect surviving cohort budget')
    if consumed != data['measured_replicas'] or offered != data['repair_executions']:
        raise ValueError('aggregate budgets differ')
    return dict(reached_looks=consumed//n, charged_repairs=offered, family_halfwidth=radius)


if __name__ == '__main__':
    print(json.dumps(check(json.loads(Path(__file__).with_name('results.json').read_text())), indent=2))
