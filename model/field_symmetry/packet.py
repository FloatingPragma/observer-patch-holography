"""Current complete-field symmetry evidence, with owned fresh replay."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

from tools import strict_json
from tools.build_completion_program import integration_source_paths
from . import readout, symmetry


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results.json'


def source_bindings():
    # Reuse the completion builder's static import-closure walker, including
    # eager initializers and function-local imports; no new custody parser.
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': 'model/field_symmetry/README.md',
        'source_receipts': ['model/interacting_field/results.json'],
        'independent_checks': [],
    })
    paths.update(('model/field_symmetry/contract.json',
                  'code/source_selection_model/response.json',
                  'code/m1_source_realization/receipt.json'))
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


_IMPORTED_BINDINGS = source_bindings()


def report():
    bindings = source_bindings()
    if bindings != _IMPORTED_BINDINGS:
        raise ValueError('source changed since import; replay in a fresh interpreter')
    # A caller can have warmed the parent producers before importing this
    # packet. Bind fresh data to fresh numerical products, not stale caches.
    readout.compiler._source.cache_clear()
    readout.compiler._native.cache_clear()
    readout.field.gaussian._local.cache_clear()
    readout.field.source._contact.cache_clear()
    readout.compiler.source_summary()
    certificate = symmetry.certificate()
    rows = [readout.execute(plane)['report'] for plane in range(3)]
    rows += [readout.execute(centered=True)['report'], readout.execute(interaction=0.)['report']]
    bounds = [F(row['numerical_dilation_budget']) for row in rows]
    # ||actual-ideal|| <= b implies probability difference <= b(2+b),
    # also when the numerical successful vector is not exactly normalized.
    allowance = max(b*(2+b) for b in bounds)
    lower = F(529, 2304)-allowance
    if lower <= F(22, 100):
        raise ArithmeticError('complete source signal no longer resolved')
    result = {
        'schema': 'oph.field_symmetry.v1',
        'campaign_id': 'OPH-FIELD-SYMMETRY-2026-10-31',
        'source_bindings': bindings,
        'symmetry': certificate,
        'complete_source_experiments': rows,
        'signal': {'ideal_probability_interval': ['529/2304', '1/4'],
                   'numerically_admitted_probability_lower': str(lower),
                   'state_specific_ideal_split_error': '0',
                   'bound_includes_fixed_numerical_admission_ceiling': True,
                   'outward_floating_or_hardware_certificate': False},
        'integrated_gate_accepted': False,
        'physical_comparison_count_added': 0,
        'physical_species_identified': False,
        'native_topological_defect_identified': False,
    }
    if bindings != source_bindings():
        raise ValueError('source changed during replay')
    return result


def compare(actual, expected, path='receipt'):
    """Strict typed tree comparison; every numerical leaf must be finite.

    Raw residuals are observations, not authenticated historical certificates.
    Fixed ceilings and independently replayed current execution own admission.
    """
    if type(actual) is not type(expected):
        raise ValueError(f'{path}: exact JSON type mismatch')
    if type(expected) is dict:
        if actual.keys() != expected.keys():
            raise ValueError(f'{path}: complete key inventory mismatch')
        for key in expected:
            compare(actual[key], expected[key], f'{path}.{key}')
    elif type(expected) is list:
        if len(actual) != len(expected):
            raise ValueError(f'{path}: complete list inventory mismatch')
        for i, (a, e) in enumerate(zip(actual, expected, strict=True)):
            compare(a, e, f'{path}[{i}]')
    elif type(expected) is float:
        if not math.isfinite(actual) or not math.isfinite(expected):
            raise ValueError(f'{path}: nonfinite diagnostic')
        # All exported floats are nonnegative probabilities, costs or norms.
        if actual < 0 or expected < 0:
            raise ValueError(f'{path}: negative probability/cost/residual')
        if not math.isclose(actual, expected, rel_tol=2e-10, abs_tol=2e-11):
            raise ValueError(f'{path}: numerical replay mismatch')
        if 'probabilit' in path and (actual > 1 or expected > 1):
            raise ValueError(f'{path}: probability exceeds one')
        if path.endswith('complete_probability_residual') and actual > readout.OUTPUT_CEILING:
            raise ValueError(f'{path}: probability accounting failed admission')
        if path.endswith('raw_primitive_residual') and actual > float(readout.PULSE_CEILING):
            raise ValueError(f'{path}: primitive failed admission')
    elif actual != expected:
        raise ValueError(f'{path}: replay value mismatch')


def verify(path=OUTPUT):
    retained = strict_json.load(Path(path))
    fresh = report()
    compare(retained, fresh)
    # Apply each retained execution's own bound, before tolerant comparison
    # can excuse a scientifically unacceptable stored value.
    for row in retained['complete_source_experiments']:
        if row['complete_output_error'] > float(F(row['numerical_dilation_budget'])):
            raise ValueError('retained complete output exceeds its budget')
        mass = math.fsum(row['terminal_probabilities'])+row['first_export_failure_probability']
        if abs(mass-1) > readout.OUTPUT_CEILING:
            raise ValueError('retained complete probability accounting failed')
    return fresh


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', nargs='?', const=str(OUTPUT))
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    if args.check is not None and args.write:
        parser.error('choose write or check')
    result = verify(args.check) if args.check is not None else report()
    if args.write:
        OUTPUT.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps({'schema': result['schema'], 'verified': args.check is not None,
                      'source_experiments': len(result['complete_source_experiments']),
                      'physical_species_identified': False}, indent=2))
