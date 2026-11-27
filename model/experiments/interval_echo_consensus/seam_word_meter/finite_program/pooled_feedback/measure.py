"""Run frozen continuum feedback banks and compare their declared-record estimates."""

import argparse
import json
from math import pi, sqrt
from pathlib import Path

import numpy as np

from ... import apparatus as app, chart
from .service import Node, flag_layout, program, run
from .publication import layout as score_layout

ACTIVE = (0, 1, 5, 7)
REST = (-.87, .33, -.19, .71, -.41, .29, -.61, .11)
FIRST = (-sqrt(2)/2, .95, -.7, .55)
TIE = (-.8, .8, -.55, .55)
FAMILIES = dict(irrational=FIRST, reflected=tuple(-v for v in FIRST), deep_low=(-.8, -.4, -.6, .8),
                second_zero=(-.75, .5, .125, .75), terminal_zero=(-.4, .8, -.4, .6),
                matched_a=TIE, matched_b=tuple(-v for v in TIE), near_low=(-.005, -.005, .8, -.8),
                near_high=(.005, .005, .8, -.8))


def values(active):
    rest = iter(REST)
    indexed = dict(zip(ACTIVE, active))
    return tuple(indexed[j] if j in indexed else next(rest) for j in range(app.N))


def independently_read_expected(panel, initial, p):
    rows = []
    for record in panel['records']:
        x = np.array(initial)
        for a, b in record['word']:
            u, v = x[[a, b]]
            x[a], x[b] = (1-p)*u+p*v, p*u+(1-p)*v
        mu = float(x[record['port']])
        measured = record['decision_read']
        lo, hi = measured['interval']
        rows.append(dict(cohort=record['cohort'], target=mu, covered=lo <= mu <= hi,
                         absolute_error=abs(measured['estimate'] - mu)))
    return rows


def report():
    seeds = iter(np.random.SeedSequence(180027).spawn(2*len(FAMILIES)))
    rows = []
    for name, angles in FAMILIES.items():
        initial = values(angles)
        for law, theta, p in [('native', pi/2, .5), ('third_turn', pi/3, .25)]:
            panel = run(program(), initial, next(seeds), angle=theta)
            rows.append(dict(family=name, law=law, panel=panel,
                             withheld_evaluation=independently_read_expected(panel, initial, p)))
    spare = np.random.SeedSequence(180028)
    invalid = list(values(FIRST))
    invalid[0] = 1.1
    refusals = [run(program(), invalid, spare), run(program(), values(FIRST), spare, replicas=4*32768-1),
                run(Node((2, 5), 0, 0), values(FIRST), spare)]
    return dict(seed=180027, simultaneous_failure_budget=.01, total_look_cap=72, cohort_size=32768,
                chart_layout=chart.stock_layout(), observer_layout=app.layout(), publication_layout=flag_layout(), score_layout=score_layout(),
                panels=rows, preflight_controls=refusals, reserved_replicas=sum(r['panel']['reserved_copies'] for r in rows),
                measured_replicas=sum(r['panel']['measured_copies'] for r in rows),
                held_replicas=sum(r['panel']['survivor_copies'] for r in rows),
                repair_executions=sum(r['panel']['scalar_repair_executions'] for r in rows))


if __name__ == '__main__':
    cli = argparse.ArgumentParser()
    cli.add_argument('--output', type=Path, required=True)
    opts = cli.parse_args()
    opts.output.write_text(json.dumps(report(), indent=1, allow_nan=False) + '\n')
