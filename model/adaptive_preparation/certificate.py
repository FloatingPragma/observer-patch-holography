"""Exact finite-attempt preparation from the actual native feedback kernel."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import numpy as np

from model.observer import native
from model.completion import observer
from model.completion.family import FamilyParameters
from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
BLOCK_LENGTH = 64
MINORIZATION = F(2, 5)


def native_counts():
    """Enumerate actual fair-coin branches; every row has denominator 120."""
    states = tuple((x, m) for x in native.states(binary=True) for m in (0, 1))
    index = {s: i for i, s in enumerate(states)}
    counts = [[0] * len(states) for _ in states]
    if len(observer.SEAMS) != 30:
        raise ValueError('thirty native seams required')
    for i, (x, memory) in enumerate(states):
        preferred = native.F if memory else native.E
        for seam in observer.SEAMS:
            weight = 1 + 30 * (seam == preferred)
            for coin in (0, 1):
                y, (_, current) = native.native(x, seam, coin)
                counts[i][index[y, int(current > 0)]] += weight
    validate_counts(counts, 120)
    return states, counts


def validate_counts(counts, denominator):
    if (type(denominator) is not int or denominator < 1 or not counts
            or any(len(row) != len(counts) for row in counts)
            or any(type(v) is not int or v < 0 for row in counts for v in row)
            or any(sum(row) != denominator for row in counts)):
        raise ValueError('nonnegative integer stochastic counts required')


def block_minorization(counts, denominator, steps, lower_bound):
    validate_counts(counts, denominator)
    if type(steps) is not int or steps < 1:
        raise ValueError('positive integer block length required')
    if not isinstance(lower_bound, F) or not 0 < lower_bound <= 1:
        raise ValueError('positive exact minorization at most one required')
    # Object dtype retains Python integers, with no overflow or floating tolerance.
    power = np.linalg.matrix_power(np.array(counts, dtype=object), steps)
    block_denominator = denominator ** steps
    if any(sum(row) != block_denominator for row in power):
        raise ArithmeticError('block kernel lost normalization')
    minima = [min(int(v) for v in power[:, j]) for j in range(len(counts))]
    mass = F(sum(minima), block_denominator)
    if mass < lower_bound:
        raise ValueError('claimed minorization exceeds exact common-column mass')
    return power, minima, mass


def reference(states, counts):
    rows = strict_json.load(ROOT / 'model/completion/stationary_reference_results.json')['stationary_joint']
    if [(tuple(r['state']), r['memory']) for r in rows] != list(states):
        raise ValueError('stationary reference domain or order changed')
    probabilities = [F(row['probability']) for row in rows]
    if min(probabilities) < 0 or sum(probabilities) != 1:
        raise ValueError('stationary reference is not normalized')
    if any(sum(p * counts[i][j] for i, p in enumerate(probabilities)) != 120 * probabilities[j]
           for j in range(len(states))):
        raise ValueError('stationary reference fails the native kernel residual')
    return probabilities


def full_record_quotient(states, counts):
    """Check every source/bit basis input of both complete endpoint instruments."""
    index = {s: i for i, s in enumerate(states)}
    checked = 0
    for c_weight in (F(0), F(1)):
        params = FamilyParameters(c_weight=c_weight, repair_probability=F(1),
                                  exploration_probability=F(1, 2), memory_capacity=1,
                                  overflow='export_oldest')
        for memory in (0, 1):
            record = observer.Record('repair', 'feedback', native.E, 'C', 0, memory, memory)
            branches = observer.branches((record,), params)
            projected = {i: {} for i, (_, m) in enumerate(states) if m == memory}
            for branch in branches:
                if not branch.weight:
                    continue
                if branch.memory != (branch.record,) or branch.exported != record or branch.halted:
                    raise ValueError('last-record instrument or retained export changed')
                for source, target in branch.pairs:
                    x, y = observer.STATES[source], observer.STATES[target]
                    if (x, memory) not in index:
                        continue
                    i = index[x, memory]
                    j = index[y, int(branch.record.current > 0)]
                    projected[i][j] = projected[i].get(j, F(0)) + branch.weight
            for i, row in projected.items():
                if [row.get(j, F(0)) for j in range(len(states))] != [F(v, 120) for v in counts[i]]:
                    raise ValueError('full record instrument differs from native bit quotient')
                checked += 1
    return checked


def budget(accuracy, *, steps=BLOCK_LENGTH, minorization=MINORIZATION):
    if not isinstance(accuracy, F) or not 0 < accuracy < 1:
        raise ValueError('exact accuracy strictly between zero and one required')
    if type(steps) is not int or steps < 1 or not isinstance(minorization, F) or not 0 < minorization <= 1:
        raise ValueError('invalid contraction certificate')
    error, blocks = F(1), 0
    while error > accuracy:
        error *= 1 - minorization
        blocks += 1
    return {'target_total_variation': str(accuracy), 'blocks': blocks,
            'source_bit_attempts': steps * blocks, 'full_last_record_attempts': steps * blocks + 1,
            'total_variation_bound': str(error)}


def report():
    states, counts = native_counts()
    probabilities = reference(states, counts)
    checked = full_record_quotient(states, counts)
    power, minima, mass = block_minorization(counts, 120, BLOCK_LENGTH, MINORIZATION)
    # The common row measure is minima/sum(minima). Scaling it by 2/5
    # leaves a nonnegative residual of row mass 3/5 in every block row.
    if any(MINORIZATION * F(minima[j], sum(minima)) > F(int(power[i, j]), 120 ** BLOCK_LENGTH)
           for i in range(len(states)) for j in range(len(states))):
        raise ArithmeticError('minorization residual is negative')
    sources = ['model/observer/native.py', 'model/completion/observer.py',
               'model/completion/observer_reference.py', 'model/completion/family.py',
               'model/observer/quantum_mixing.py', 'model/repair_scheduler/coins.py',
               'model/completion/stationary_reference_results.json']
    sources += [p.relative_to(ROOT).as_posix() for p in sorted(HERE.iterdir()) if p.suffix in ('.py', '.md')]
    return {'schema': 'oph.adaptive_preparation.v1',
            'source_bindings': {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in sources},
            'kernel': {'states': len(states), 'stationary_support': sum(p > 0 for p in probabilities),
                       'transition_denominator': 120, 'block_attempts': BLOCK_LENGTH,
                       'common_column_numerators': minima, 'block_denominator': 120 ** BLOCK_LENGTH,
                       'common_column_mass': str(mass), 'certified_minorization': str(MINORIZATION),
                       'contraction_per_block': str(1 - MINORIZATION)},
            'complete_record_basis_inputs_checked': checked,
            'budgets': [budget(F(1, 10 ** n)) for n in (3, 6, 9)],
            'scope': {'initial_states': 'Every classical distribution on binary66 times one current-sign bit.',
                      'observer': 'r=1, exploration=1/2, one last record, export_oldest; C, E and fixed tagged mixtures.',
                      'distance': 'Total variation on the finite source/bit or source/last-record marginal.',
                      'archive': 'All older records are exported; the bound traces their marginal and does not erase them.',
                      'quantum': 'No contraction of arbitrary off-diagonal source states or quantum references is claimed.',
                      'physical': 'Native attempt counts; no physical duration, modular apparatus or particle dictionary.',
                      'formalization': 'Exact executable certificate and elementary contraction proof; no additional Lean theorem claimed.'}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    value = report()
    path = HERE / 'results.json'
    if args.check and strict_json.load(path) != value:
        raise ValueError('adaptive preparation receipt differs from exact replay')
    if args.write:
        path.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'minorization': value['kernel']['certified_minorization'], 'budgets': value['budgets']}, indent=2))
