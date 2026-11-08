"""Exact finite continuity certificate attached to complete native readout."""
import argparse
from fractions import Fraction as F
import hashlib
from itertools import product
from pathlib import Path

from model.native_fair_local_readout import packet as inherited
from tools import strict_json
from . import history as h

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results.json'


def source_bindings():
    paths = set(inherited.source_bindings())
    paths.update(p.relative_to(ROOT).as_posix() for p in HERE.iterdir()
                 if p.suffix in ('.py', '.md'))
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def produce():
    before = source_bindings()
    p = h.r.plan((1, 2), ((0, 1), (1, 2)))
    edges = ((0, 1), (1, 2), (0, 1))
    weights = (F(2), F(-3, 2))
    blocks = {}
    for initial in range(8):
        for coins in product((0, 1), repeat=3):
            state = initial
            events = []
            for tick, (edge, coin) in enumerate(zip(edges, coins, strict=True), 1):
                a, b = edge
                rec = h.r.record(edge, coin, ((state >> a) & 1)-((state >> b) & 1))
                state += rec.current*((1 << b)-(1 << a))
                events.append((tick, rec))
            blocks.setdefault(tuple(events), []).append((initial, state))
    column_weights = [F()]*8
    decoded = erasures = 0
    for events, pairs in blocks.items():
        observed = h.observe(p, events, 3)
        result = h.ward(p, observed, weights)
        decoded += observed['absolute_history_reconstructed']
        erasures += not observed['absolute_history_reconstructed']
        for initial, last in pairs:
            actual_delta = sum((weights[i]*(((last >> v) & 1)-((initial >> v) & 1))
                                for i, v in enumerate(p['ports'])), F())
            if result['net_change'] != actual_delta:
                raise ArithmeticError('exact fine-branch operator continuity failed')
            column_weights[initial] += F(1, 8)
            if observed['absolute_history_reconstructed']:
                if (initial >> 1) & 3 != observed['decoded']['initial_mask']:
                    raise ArithmeticError('initial projection identity failed')
    if column_weights != [F(1)]*8:
        raise ArithmeticError('complete branch sum changed')
    value = {
        'schema': 'oph.native_local_charge_history.v1',
        'source_bindings': before,
        'finite_control': {'input_columns': 8, 'live_fine_words': len(blocks),
                           'decoded_words': decoded, 'erasure_words': erasures,
                           'complete_column_weights': [str(x) for x in column_weights],
                           'branch_operator_continuity_checked_exactly': True},
        'history_error': 'same complete local-readout error; no snapshot-count factor',
        'reference_correlations_allowed': True,
        'physical_current_identification': False,
        'physical_particle_identification': False,
        'stress_energy_tensor_derived': False,
        'expanded_large_source_execution_asserted': False,
    }
    if source_bindings() != before:
        raise RuntimeError('source changed during production')
    return value


def check():
    if inherited.canonical(strict_json.load(OUTPUT)) != inherited.canonical(produce()):
        raise ValueError('charge-history receipt differs from fresh exact replay')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.write:
        OUTPUT.write_text(inherited.canonical(produce()))
    else:
        check()
    print('Native local charge-history certificate passed.')
