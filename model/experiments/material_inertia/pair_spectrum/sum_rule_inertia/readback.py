"""Solver-free inference from terminal counts and numerical certificates."""

import json
from math import asin, log, sqrt, sin
from pathlib import Path

HERE = Path(__file__).resolve().parent


def check(ok, message):
    if not ok:
        raise ArithmeticError(message)


def population(record):
    counts = record['counts'][:-1]
    n = sum(counts)
    error = sqrt(log(2 * 200 * 20 / .01) / (2 * n))
    result = []
    for c, eps in zip(counts, record['amplitude_systematic_bounds'], strict=True):
        p = c / n
        result.append([max(0, sqrt(max(0, p - error)) - eps)**2,
                       min(1, (sqrt(min(1, p + error)) + eps)**2)])
    return result


def response(record, pop, edges, profile):
    counts = record['counts']
    n = sum(counts[:4])
    tau = record['tau']
    seen = (counts[0] - counts[1]) / n / tau
    err = sqrt(2 * log(2 * 200 * 20 / .01) / n) / tau
    err += (record['operator_approximation_bound']
            + record['response_operator_width'] * record['coefficient_error']
            + record['rest_energy_width'])
    argument = record['band_energy_width'] * tau
    check(0 < argument < 1.5708, 'unbounded sine inversion')
    sinc = sin(argument) / argument
    moment = [max(0, seen - err), (seen + err) / sinc]
    d = [0., 0.]
    for a, b in edges:
        square = (profile[a] - profile[b])**2
        for i in (0, 1):
            d[i] += sqrt(pop[a][i] * pop[b][i]) * square
    check(d[0] > 0, 'unresolved population geometry')
    return [moment[0] / d[1], moment[1] / d[0]]


def infer(raw, key):
    rows = []
    sources = raw[key]
    calibration = None
    for source in sources:
        pop = population(source['population'])
        for law, panels in source['response'].items():
            for i, panel in enumerate(panels):
                intervals = [response(r, pop, raw['owner_edges'], raw['profiles'][i]['values'])
                             for r in panel['records']]
                interval = [max(p[0] for p in intervals), min(p[1] for p in intervals)]
                check(interval[0] <= interval[1], 'clock intervals disjoint')
                if calibration is None:
                    calibration = interval
                    outcome = 'calibration'
                else:
                    outcome = ('reject' if interval[1] < .95 * calibration[0]
                               or interval[0] > 1.05 * calibration[1] else 'compatible')
                rows.append(dict(source=source['label'], law=law, shell=panel['shell'], center=panel['center'],
                                 kinetic_coefficient_interval=interval, outcome=outcome))
    return dict(calibration=calibration, operational_mass_interval=[1 / (2 * calibration[1]), 1 / (2 * calibration[0])],
                relative_law_tolerance=.05, held=rows)


def close(a, b):
    if isinstance(b, dict):
        check(a.keys() == b.keys(), 'readback keys')
        for k in b:
            close(a[k], b[k])
    elif isinstance(b, list):
        check(len(a) == len(b), 'readback length')
        for c, d in zip(a, b, strict=True):
            close(c, d)
    elif isinstance(b, float):
        check(abs(a - b) <= 2e-12 * max(1e-9, abs(b)), f'count readback differs: {a} vs {b}')
    else:
        check(type(a) is type(b) and a == b, 'count readback verdict differs')


def clock(raw):
    rows = []
    for r in raw['clock_standard']['records']:
        n = sum(r['counts'])
        p = r['counts'][1] / n
        err = sqrt(log(2 * 200 * 20 / .01) / (2 * n))
        gap = [2 * asin(sqrt(max(0., p - err))) / r['commanded_appointment'],
               2 * asin(sqrt(min(1., p + err))) / r['commanded_appointment']]
        close(gap, r['native_gap_interval'])
        check(gap[0] <= 1 <= gap[1], 'reference-clock gap no longer calibrated')
        rows.append(gap)
    check(all(max(rows[0][0], r[0]) <= min(rows[0][1], r[1]) for r in rows[1:]), 'held reference-clock failure')
    return rows


def main():
    raw = json.loads((HERE / 'results.json').read_text())
    clock(raw)
    close(infer(raw, 'native_sources'), raw['native_inference'])
    close(infer(raw, 'scalar_controls'), raw['scalar_inference'])
    print('Recomputed calibration, held scalar-inertia failures and controls from counts.')


if __name__ == '__main__':
    main()
