"""Read acceleration from calibrated marker and clock records on a finite layered aperture."""

import argparse
from collections import Counter
import gzip
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np

from model.candidate_laws.interval_echo import clock
from model.executable_model import federation


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCE_LOGS = ROOT / 'evidence/source_net_causal_poset/carrier_source_net_logs'
MEAN_COUNT = 262144
CLOCK_UNITS = 10**9
SEED = 2026112313
HISTORIES = ('dust', 'coasting', 'de_sitter')


def layered_aperture(q):
    path = SOURCE_LOGS / f'q{q}_event_log.json.gz'
    with gzip.open(path, 'rt') as stream:
        log = json.load(stream)
    writes = {(row['write'][0], row['write'][1]): row for row in log['events']}
    if len(writes) != len(log['events']):
        raise ValueError('duplicate register versions')
    initial = (row for row in log['events'] if row['event'][0] == 1)
    receiver = max(initial, key=lambda row: (len(row['reads']), -row['event'][1]))['event'][1]
    aperture = {address for address, _ in writes[receiver, 2]['reads']}
    for layer in range(1, log['rounds'] + 1):
        record = writes[receiver, layer + 1]
        if set(map(tuple, record['reads'])) != {(address, layer) for address in aperture}:
            raise ValueError('changing aperture or parent version')
        for parent in record['reads']:
            if writes[tuple(parent)]['event'] != [layer - 1, parent[0]]:
                raise ValueError('incoming version does not reference its preceding local write')
    local_writes = Counter(layer for layer, address in (row['event'] for row in writes.values())
                           if address in aperture)
    return dict(q=q, receiver=receiver, aperture_size=len(aperture),
                raw_writes=[local_writes[i] for i in range(log['rounds'] + 1)],
                source_sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def history(name, eta, rate=1):
    scaled = rate * eta
    if name == 'dust':
        return (1 + scaled / 2)**2, (2 / (3 * rate)) * ((1 + scaled / 2)**3 - 1)
    if name == 'de_sitter':
        return 1 / (1 - scaled), -math.log1p(-scaled) / rate
    if name == 'coasting':
        return math.exp(scaled), math.expm1(scaled) / rate
    if name == 'static':
        return 1.0, eta
    raise ValueError(name)


def marker_counts(name, step_denominator, aperture_size, q, *, changing_density=False, initial_rate=1, stream=None):
    marker_history = 'de_sitter' if changing_density else name
    seed_parts = [SEED, q, step_denominator, HISTORIES.index(marker_history)]
    if stream is not None:
        seed_parts.append(stream)
    independent_draws = np.random.default_rng(np.random.SeedSequence(seed_parts))
    counts = []
    for j in range(4):
        epoch = Fraction(j, step_denominator)
        epoch *= initial_rate
        if name == 'coasting':
            scale = math.exp(float(epoch))
        elif name == 'dust':
            scale = (1 + epoch / 2)**2
        else:
            scale = 1 / (1 - epoch)
        density = MEAN_COUNT
        if changing_density:
            if name != 'dust':
                raise ValueError('the density intervention requires the dust geometry')
            de_sitter = 1 / (1 - epoch)
            density *= (de_sitter / scale)**4
        flat_cell_volume = Fraction(1, q**3)
        conformal_tick = Fraction(1, step_denominator)
        preparations = q**3 * step_denominator
        rate = float(preparations * density * conformal_tick * flat_cell_volume * scale**4)
        sites = independent_draws.poisson(rate, size=aperture_size)
        counts.append(int(sites.sum()))
    return counts


def clock_intervals(name, step_denominator, initial_rate=1):
    ticks = [round(CLOCK_UNITS * history(name, j / step_denominator, initial_rate)[1]) for j in range(4)]
    return [right - left for left, right in zip(ticks, ticks[1:])]


def decode(counts, clock_ticks, aperture_size, conformal_step):
    if len(counts) != 4 or len(clock_ticks) != 3 or min(counts) <= 0 or min(clock_ticks) <= 0:
        raise ValueError('four positive totals and three positive clock intervals are required')
    scales = [(total / (MEAN_COUNT * aperture_size))**0.25 for total in counts]
    scale_errors = [scale / (4 * math.sqrt(total)) for scale, total in zip(scales, counts)]
    intervals = [ticks / CLOCK_UNITS for ticks in clock_ticks]
    slopes = [(right - left) / dt for left, right, dt in zip(scales, scales[1:], intervals)]
    acceleration = []
    for i in (1, 2):
        left, right = intervals[i - 1:i + 1]
        speed = (right * slopes[i - 1] + left * slopes[i]) / (left + right)
        second = 2 * (slopes[i] - slopes[i - 1]) / (left + right)
        acceleration.append(round(-scales[i] * second / speed**2, 7) if speed else None)
    residuals = []
    tolerances = []
    for i, observed in enumerate(intervals):
        predicted = conformal_step * (scales[i] + scales[i + 1]) / 2
        margin = 9 * conformal_step**3 / 12
        margin += 6 * conformal_step * max(scale_errors[i:i + 2]) + 1 / CLOCK_UNITS
        residuals.append(abs(observed - predicted))
        tolerances.append(margin)
    return dict(scale_estimates=[round(x, 9) for x in scales],
                deceleration_parameters=acceleration,
                clock_absolute_residuals=[round(x, 10) for x in residuals],
                clock_margins=[round(x, 10) for x in tolerances],
                clock_intervals_accepted=[r <= m for r, m in zip(residuals, tolerances)])


def observe(q, size, denominator, profile, *, alias=False):
    totals = marker_counts(profile, denominator, size, q, changing_density=alias)
    times = clock_intervals(profile, denominator)
    return dict(marker_totals=totals, proper_clock_intervals=times,
                reader=decode(totals, times, size, 1 / denominator))


def measure_constant_expansion(totals, clock_ticks):
    # Ratios remove the absolute marker normalization.
    relative_scales = [(n / totals[1])**0.25 for n in totals]
    intervals = [dt / CLOCK_UNITS for dt in clock_ticks]
    estimate = -math.log(relative_scales[0]) / intervals[0]
    errors = [math.exp(estimate * sum(intervals[1:j])) / relative_scales[j] - 1
              for j in (2, 3)]
    compatible = max(abs(e) for e in errors) < 0.002
    return dict(hubble_rate=round(estimate, 8),
                fractional_errors_on_held_scales=[round(e, 8) for e in errors],
                passes_constant_rate_forecast=compatible,
                de_sitter_lambda_in_clock_units=round(3 * estimate**2, 8) if compatible else None)


def amplitude_trials(q, size):
    trials = []
    menu = [('de_sitter', Fraction(1, 2)), ('de_sitter', Fraction(3, 4)),
            ('de_sitter', Fraction(5, 4)), ('dust', Fraction(1))]
    for i, (profile, initial_rate) in enumerate(menu):
        counts = marker_counts(profile, 16, size, q, initial_rate=initial_rate, stream=101 + i)
        intervals = clock_intervals(profile, 16, float(initial_rate))
        trials.append(dict(preparation=profile, prepared_initial_rate=str(initial_rate),
                           marker_totals=counts, proper_clock_intervals=intervals,
                           forecast=measure_constant_expansion(counts, intervals)))
    return trials


def native_publications():
    panels = []
    for level in range(3):
        car = federation(level)
        trip = clock.material_replay(car, (8, 1, 2, 4, 16, 2), age_bits=3, seed=SEED)
        result = dict(cells=car.cells, services=len(trip['local_events']))
        for owner, digit, oscillator in (('home', 'home_tick', 'home_oscillator_tick'),
                                          ('traveller', 'tick', 'oscillator_tick')):
            events = trip['local_events']
            bits = [e[digit] for e in events]
            elapsed = [0] + [e[owner] for e in events]
            increments = [b - a for a, b in zip(elapsed, elapsed[1:])]
            result[owner] = dict(positives=sum(bits), quiet=len(bits) - sum(bits),
                                 reunion_ticks=trip['reunion'][owner],
                                 count_mismatches=sum(x != y for x, y in zip(bits, increments)),
                                 trigger_mismatches=sum(e[digit] != e[oscillator] for e in events))
        panels.append(result)
    return panels


def run():
    comparisons = []
    for q in (5, 8):
        source = layered_aperture(q)
        size = source['aperture_size']
        panels = []
        for denominator in (8, 16):
            readings = {name: observe(q, size, denominator, name) for name in HISTORIES}
            readings['dust_with_varying_density'] = observe(q, size, denominator, 'dust', alias=True)
            static = [MEAN_COUNT * size] * 4
            readings['exact_static_control'] = decode(static, clock_intervals('static', denominator),
                                                      size, 1 / denominator)
            panels.append(dict(conformal_step=f'1/{denominator}', histories=readings))
        comparisons.append(dict(aperture=source, panels=panels, amplitude_trials=amplitude_trials(q, size)))
    return dict(protocol_sha256=hashlib.sha256((HERE / 'PROTOCOL.md').read_bytes()).hexdigest(),
                flat_count_mean_per_site=MEAN_COUNT, proper_clock_ticks_per_time_unit=CLOCK_UNITS,
                seed=SEED, layered_comparison=comparisons, interval_echo_clocks=native_publications())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='write an independent receipt to a new path')
    args = parser.parse_args()
    result = run()
    if args.output:
        with args.output.open('x') as stream:
            stream.write(json.dumps(result, indent=2, allow_nan=False) + '\n')
    elif result != json.loads((HERE / 'results.json').read_text()):
        raise SystemExit('expansion reconstruction differs from the published receipt')
    print('Reconstructed four comparison panels and checked three material-clock boards.')
