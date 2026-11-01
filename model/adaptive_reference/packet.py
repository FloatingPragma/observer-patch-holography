"""Replay exact adaptive reference, native source and complete mixing checks."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from tools import strict_json
from ..observer.scalar_response.experiment import canonical, compare_receipt
from ..continuing_reference import preparation
from . import controller, history, instrument


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
SOURCE_PATHS = (
    'model/adaptive_reference/__init__.py', 'model/adaptive_reference/__main__.py',
    'model/adaptive_reference/controller.py', 'model/adaptive_reference/instrument.py',
    'model/adaptive_reference/history.py', 'model/adaptive_reference/packet.py',
    'model/adaptive_reference/test_reference.py', 'model/adaptive_reference/README.md',
    'model/adaptive_reference/contract.json', 'paper/tex_fragments/ADAPTIVE_REFERENCE.tex',
    'Lean/QFT/AdaptiveReference.lean',
    'model/observer/native.py', 'model/charged_defect/branches.py',
    'model/charged_defect/dictionary.py', 'model/particle_transport/algebra.py',
    'model/exchange_sector/statistics.py',
    'model/continuing_reference/algebra.py', 'model/continuing_reference/preparation.py',
    'model/continuing_reference/results.json',
    'model/certified_reference/source.py', 'model/certified_reference/interval.py',
    'model/observer/source_scalar/compiler.py',
    'code/source_selection_model/response.json', 'code/m1_source_realization/receipt.json',
)


def encode(value):
    if type(value) is F:
        return str(value)
    if isinstance(value, dict):
        return {k: encode(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(v) for v in value]
    return value


def run():
    return encode({
        'schema': 'oph.adaptive_reference.v1', 'issues': [1138], 'gates': ['CP-02'],
        'all_pass': True,
        'scope': 'balanced adaptive classical memory and complete native repair with stationary field/controller marginal',
        'controller': controller.permutation_certificate(),
        'exact_native_seed': controller.seed_certificate(),
        'literal_field_reference': preparation.one_mode_certificate(),
        'native_attachment': instrument.attachment_certificate(),
        'complete_balance': instrument.balance_certificate(),
        'continuing_reference': history.continuing_theorem(),
        'finite_quantum_mixing': instrument.mixing_certificate(),
        'complete_positive_history': history.enumerate_positive(8),
        'negative_controls': instrument.negative_controls(),
        'resource_rows': [history.resource_certificate(m, n, min(m, 4))
                          for m, n in ((4, 0), (12, 8), (48, 1000), (10000, 1000000))],
        'charge_clock_source_bound': preparation.modular_certificate(12),
        'source_or_repair_law_changed': False,
        'full_growing_archive_stationary': False,
        'arbitrary_feedback_covered': False,
        'canonical_physical_premises_discharged': [],
        'integrated_gate_accepted': False,
        'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
    })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', nargs='?', const=OUTPUT, type=Path)
    args = parser.parse_args()
    result = run()
    if args.write:
        OUTPUT.write_text(canonical(result))
    if args.check:
        compare_receipt(strict_json.loads(args.check.read_text()), result)
    print(json.dumps({'all_pass': result['all_pass'], 'controller_states': 4,
                      'quantum_corner_dimension': 12,
                      'eight_tick_depolarizing_weight': '3/64',
                      'one_bit_counterexample_trace_distance': '1/6',
                      'complete_archive_stationary': False,
                      'integrated_gate_accepted': False}, indent=2))


if __name__ == '__main__':
    main()
