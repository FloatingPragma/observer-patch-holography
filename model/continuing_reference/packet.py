"""Replay the scalable source-reference theorem, bounds and negative controls."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from tools import strict_json
from ..observer.scalar_response.experiment import canonical, compare_receipt
from . import algebra, preparation


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
SOURCE_PATHS = (
    'model/continuing_reference/__init__.py', 'model/continuing_reference/__main__.py',
    'model/continuing_reference/algebra.py', 'model/continuing_reference/preparation.py',
    'model/continuing_reference/packet.py', 'model/continuing_reference/test_reference.py',
    'model/continuing_reference/README.md', 'model/continuing_reference/contract.json',
    'model/certified_reference/source.py', 'model/certified_reference/interval.py',
    'model/observer/source_scalar/compiler.py', 'model/observer/native.py',
    'model/charged_defect/branches.py', 'model/charged_defect/dictionary.py',
    'model/particle_transport/algebra.py', 'model/exchange_sector/statistics.py',
    'model/common_spin_protection/spin.py', 'model/common_spin_protection/spectrum.py',
    'model/field_stress/energy.py', 'model/field_symmetry/symmetry.py',
    'model/full_fock_field/field.py', 'model/repair_spectral_sectors/algebra.py',
    'model/repair_spectral_sectors/results.json',
    'code/source_selection_model/response.json', 'code/m1_source_realization/receipt.json',
    'Lean/QFT/ContinuingReference.lean',
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
        'schema': 'oph.continuing_reference.v1',
        'issues': [1138], 'gates': ['CP-02'], 'all_pass': True,
        'scope': 'finite full-Fock source reference and complete record-independent continuing process',
        'source_preparation': preparation.one_mode_certificate(),
        'resource_rows': [preparation.resource_certificate(m, min(m, 4)) for m in (4, 12, 48, 10000)],
        'source_execution_diagnostic': preparation.numerical_diagnostic(),
        'modular_source_rows': [preparation.modular_certificate(m) for m in (4, 12, 48)],
        'exact_number_and_unitality_premises': algebra.exact_number_premises(),
        'continuing_reference': algebra.history_theorem(),
        'fixed_sector_reference_classification': algebra.fixed_reference_classification(),
        'record_feedback_counterexample': algebra.feedback_counterexample(),
        'separate_aligned_spin_target': algebra.aligned_spin_reference(),
        'full_space_stationarity_diagnostics': algebra.numerical_stationarity(),
        'source_law_or_native_repair_changed': False,
        'successful_field_substate_silently_normalized': False,
        'physical_clock_or_energy_temperature_identified': False,
        'canonical_physical_premises_discharged': [],
        'integrated_gate_accepted': False,
        'source_sha256': {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in SOURCE_PATHS},
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
    print(json.dumps({
        'all_pass': result['all_pass'],
        'exact_stationary_field_weight': result['continuing_reference']['reference_success_block'],
        'one_mode_complete_error_bound': float(F(result['source_preparation']['complete_one_mode_isometry_error'])),
        'active_energy_variance': result['fixed_sector_reference_classification']['active_K_variance'],
        'feedback_counterexample_distance': result['record_feedback_counterexample']['reference_trace_distance'],
        'integrated_gate_accepted': result['integrated_gate_accepted'],
    }, indent=2))


if __name__ == '__main__':
    main()
