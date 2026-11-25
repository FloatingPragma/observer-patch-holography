"""Read every proper template half-turn on the finite energy-filtered trimer."""

import argparse
import gzip
import json
from math import log, sqrt
from pathlib import Path

import numpy as np

from model.experiments.interval_echo_encounter.gateway import counts_law, sample
from ..measure import compile_turn
from . import measure as exchange

STARTS = 2**32
SEED = 807332
RECEIPT = exchange.HERE / 'half_turns.json.gz'


def record(refused, admission, conditional, rng):
    probabilities = np.r_[refused, admission * conditional]
    exchange.check(np.all(probabilities >= 0) and abs(probabilities.sum() - 1) < 5e-13,
                   'Incomplete trimer detector')
    weights = counts_law(probabilities)
    error = float(np.abs(weights / weights.sum() - probabilities).sum())
    exchange.check(error < 3e-11, 'Sampling law moved outside its stated allowance')
    return sample(weights, STARTS, rng)


def measure():
    sec, _, _, _, _, _, finite, source = exchange.source.source()
    words = np.asarray(sec.states)
    masses = np.asarray(source['filter_acceptances'])
    refused = masses[:-1] - masses[1:]
    sharp = np.zeros(len(words), complex)
    sharp[sec.index[tuple(source['source_word'])]] = 1
    identity = tuple(range(12))
    turns = [p for p in sorted(exchange.rotations()) if tuple(p[i] for i in p) == identity]
    exchange.check(len(turns) == 16 and turns[0] == identity, 'Incomplete proper half-turn census')

    rng = np.random.default_rng(SEED)
    panels = []
    for perm in turns:
        contacts = compile_turn(perm)
        moved = words.copy()
        unequal = np.zeros(len(words), int)
        for a, b in contacts:
            unequal += np.any(moved == a, axis=1) ^ np.any(moved == b, axis=1)
            moved = np.where(moved == a, b, np.where(moved == b, a, moved))
        exchange.check(np.array_equal(moved, np.asarray(perm)[words]), 'Incorrect physical rotation')
        destinations = [sec.index[tuple(sorted(w))] for w in moved]
        settings = {}
        for name, initial, failures, acceptance in (
            ('selected', finite, refused, masses[-1]),
            ('polar', sharp, np.zeros(4), 1.),
        ):
            turned = np.empty_like(initial)
            turned[destinations] = initial
            address = np.abs(turned)**2
            phase_counts = {}
            for contrast, gate_phase in (('native_8', 1), ('changed_7_5', (1j)**(unequal % 4))):
                turned[destinations] = initial * gate_phase
                phase_ports = np.array([np.linalg.norm((initial + sign * turned) / 2)**2
                                        for sign in (1, -1)])
                phase_counts[contrast] = record(failures, acceptance, phase_ports, rng)
            settings[name] = dict(address_counts=record(failures, acceptance, address, rng),
                                  phase_counts=phase_counts)
        panels.append(dict(ports=list(perm), contacts=[list(edge) for edge in contacts], **settings))
    return dict(starts=STARTS, seed=SEED, words=words.tolist(), panels=panels,
                record_order='four source refusals then all address triples, or plus/minus',
                conditional_sampling_l1_upper=3e-11 / masses[-1] * 2,
                source_admission=float(masses[-1]))


def address_separation(data):
    """Lower bounds on total variation from individual, simultaneously covered bins."""
    panels = data['panels']
    comparisons = len(panels) * len(data['words']) * 2  # Selected and sharp sources.
    coverage_log = log(2 * comparisons / 0.01)
    bounds = {}
    for name in ('selected', 'polar'):
        intervals = []
        for panel in panels:
            counts = np.asarray(panel[name]['address_counts'][4:])
            total = counts.sum()
            exchange.check(total > 0, 'No prepared source for this comparison')
            radius = sqrt(coverage_log / (2 * total)) + data['conditional_sampling_l1_upper']
            intervals.append((counts / total, radius))
        idle, idle_error = intervals[0]
        bounds[name] = [float(np.max(np.abs(active - idle)) - error - idle_error)
                        for active, error in intervals[1:]]
    return bounds


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    data = measure()
    serialized = json.dumps(data, separators=(',', ':'), allow_nan=False).encode()
    if args.output:
        with args.output.open('xb') as stream:
            stream.write(gzip.compress(serialized, mtime=0))
    else:
        old = json.loads(gzip.decompress(RECEIPT.read_bytes()))
        exchange.compare(data, old, 'half_turns')
        lower = min(address_separation(old)['selected'])
        exchange.check(lower > 0, 'A proper half-turn is unresolved')
        print(f'All 15 proper half-turns distinguished by address counts, lower bound {lower:.8f}.')


if __name__ == '__main__':
    main()
