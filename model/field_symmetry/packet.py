"""Current complete-field symmetry evidence, with owned fresh replay."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile

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


def _produce():
    """Worker calculation; public evidence is produced in a fresh interpreter."""
    bindings = source_bindings()
    if bindings != _IMPORTED_BINDINGS:
        raise ValueError('source changed since import; replay in a fresh interpreter')
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


def report():
    """Own cold replay, including code imported before this packet was loaded.

    A distinct empty bytecode-cache location matters: -B alone still reads
    timestamp-valid stale pyc files. Parent and worker compare the complete
    manifest before and after execution; no in-process fallback is allowed.
    """
    bindings = source_bindings()
    if bindings != _IMPORTED_BINDINGS:
        raise ValueError('source changed since import; replay in a fresh interpreter')
    bootstrap = ('import runpy,sys; sys.path.insert(0,sys.argv.pop(1)); '
                 'runpy.run_module("model.field_symmetry",run_name="__main__")')
    with tempfile.TemporaryDirectory(prefix='oph-field-symmetry-') as cache:
        child = subprocess.run(
            [sys.executable, '-I', '-B', '-W', 'error', '-X', 'utf8',
             '-X', f'pycache_prefix={cache}', '-c', bootstrap, str(ROOT), '--worker'],
            cwd=ROOT, capture_output=True, text=True, encoding='utf-8', timeout=180)
    if child.returncode:
        raise ValueError('fresh source replay failed: '+child.stderr[-4000:])
    result = strict_json.loads(child.stdout)
    keys = {'schema', 'campaign_id', 'source_bindings', 'symmetry',
            'complete_source_experiments', 'signal', 'integrated_gate_accepted',
            'physical_comparison_count_added', 'physical_species_identified',
            'native_topological_defect_identified'}
    if (type(result) is not dict or result.keys() != keys
            or result['schema'] != 'oph.field_symmetry.v1'):
        raise ValueError('fresh source replay returned an incomplete payload')
    if result['source_bindings'] != bindings or source_bindings() != bindings:
        raise ValueError('source changed across fresh replay')
    compare(result, result)  # Also reject nonfinite diagnostics from child JSON.
    _worker_inventory(result)
    return result


def _worker_inventory(result):
    """Validate the fixed campaign's output protocol, not a second producer."""
    flags = ('integrated_gate_accepted', 'physical_species_identified',
             'native_topological_defect_identified')
    if (result['campaign_id'] != 'OPH-FIELD-SYMMETRY-2026-10-31'
            or any(result[key] is not False for key in flags)
            or type(result['physical_comparison_count_added']) is not int
            or result['physical_comparison_count_added'] != 0):
        raise ValueError('fresh replay changed the campaign or physical scope')
    certificate = result['symmetry']
    if (type(certificate) is not dict
            or certificate.get('schema') != 'oph.full_field_spin_eta_symmetry.v1'
            or type(certificate.get('frame')) is not dict
            or not certificate['frame']
            or type(certificate.get('full_carrier_checks')) is not list
            or len(certificate['full_carrier_checks']) != 2
            or type(certificate.get('full_carrier_irreps')) is not dict
            or certificate['full_carrier_irreps'].keys() != {'1', '2', '4', '8'}
            or any(type(rows) is not list or not rows
                   for rows in certificate['full_carrier_irreps'].values())
            or any(type(row) is not dict or row.get('complete_Fock_dimension') != dimension
                   for row, dimension in zip(certificate['full_carrier_checks'], (256, 65536), strict=True))
            or certificate.get('physical_particle_or_topological_defect_identified') is not False
            or certificate.get('interacting_bulk_gap_or_localized_spinor_proved') is not False):
        raise ValueError('fresh replay omitted the complete symmetry inventory')
    signal = result['signal']
    if (type(signal) is not dict or signal.keys() != {
            'ideal_probability_interval', 'numerically_admitted_probability_lower',
            'state_specific_ideal_split_error', 'bound_includes_fixed_numerical_admission_ceiling',
            'outward_floating_or_hardware_certificate'}
            or signal['bound_includes_fixed_numerical_admission_ceiling'] is not True
            or signal['outward_floating_or_hardware_certificate'] is not False):
        raise ValueError('fresh replay omitted the signal inventory')
    rows = result['complete_source_experiments']
    cases = [(0, 'H', 1.), (1, 'H', 1.), (2, 'H', 1.), (0, 'K=H-U*N/2', 1.), (0, 'H', 0.)]
    if type(rows) is not list or len(rows) != len(cases):
        raise ValueError('fresh replay omitted an experiment')
    for row, (plane, generator, interaction) in zip(rows, cases, strict=True):
        if (type(row) is not dict or type(row.get('plane')) is not int
                or type(row.get('interaction')) is not float
                or (row['plane'], row.get('generator'), row.get('interaction')) != (plane, generator, interaction)
                or type(row.get('sides')) is not list
                or any(type(side) is not int for side in row['sides'])
                or row.get('sides') != [1 if axis == plane else 2 for axis in range(3)]
                or type(row.get('steps')) is not int or row['steps'] != 1
                or row.get('duration') != '1'
                or type(row.get('terminal_probabilities')) is not list
                or len(row['terminal_probabilities']) != 256
                or type(row.get('read_one_probabilities')) is not list
                or len(row['read_one_probabilities']) != 4):
            raise ValueError('fresh replay changed the complete experiment inventory')


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
        # These gates use the fresh child's ceilings, never potentially stale
        # readout globals in the calling interpreter.
        for residual, ceiling in (('raw_primitive_residual', 'numerical_primitive_ceiling'),
                                  ('complete_probability_residual', 'numerical_probability_mass_ceiling'),
                                  ('full_field_word_success_error', 'full_field_word_success_budget')):
            if residual in expected and actual[residual] > float(F(expected[ceiling])):
                raise ValueError(f'{path}.{residual}: failed fresh numerical admission')
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
    elif actual != expected:
        raise ValueError(f'{path}: replay value mismatch')


def _ideal_terminal_probabilities(current):
    """Owned product law, independent of retained probabilities and targets."""
    volume = math.prod(current['sides'])
    phase = (0. if current['generator'] == 'K=H-U*N/2' else
             -current['interaction']*float(F(current['duration'])))
    empty, paired = math.cos(phase/2)**2, math.sin(phase/2)**2
    ideal = []
    for address in range(1 << (2*volume)):
        local = [(address >> (2*site)) & 3 for site in range(volume)]
        ideal.append(0. if any(value not in (0, 3) for value in local) else
                     math.prod(paired if value == 3 else empty for value in local))
    return ideal


def _admit_retained_probabilities(retained, current):
    """Necessary complete-output amplitude gate, with numerical guards.

    For any phases, ||psi-ideal||^2 is at least the squared Hellinger
    amplitude distance of their terminal probabilities.  Orthogonal abort
    outputs add their full probability.  This uses the fresh fixed budget;
    a retained raw residual is not an authority for a different budget.

    The small additive amplitude guard covers the four product factors,
    square roots and accumulation in this 256-address numerical check.
    It is not an outward arithmetic certificate and never floors a cell.
    Marginal sums use the standard gamma_n bound for nonnegative summation,
    allowing the producer's sum and the separate fsum's final rounding.
    """
    probabilities = retained['terminal_probabilities']
    failure = retained['first_export_failure_probability']
    ideal = _ideal_terminal_probabilities(current)
    amplitude_distance = math.sqrt(math.fsum(
        (math.sqrt(value)-math.sqrt(target))**2
        for value, target in zip(probabilities, ideal, strict=True))+failure)
    epsilon = 2.**-52
    amplitude_guard = 64*epsilon*math.sqrt(len(ideal))
    budget = float(F(current['numerical_dilation_budget']))
    if amplitude_distance > budget+amplitude_guard:
        raise ValueError('retained complete distribution exceeds its amplitude budget')
    mass = math.fsum(probabilities)+failure
    if abs(mass-1) > float(F(current['numerical_probability_mass_ceiling'])):
        raise ValueError('retained complete probability accounting failed')
    terms = len(probabilities)//2
    gamma = (terms+2)*epsilon/(1-(terms+2)*epsilon)
    for site, retained_marginal in enumerate(retained['read_one_probabilities']):
        marginal = math.fsum(value for address, value in enumerate(probabilities)
                             if (address >> (2*site)) & 1)
        # Scale by this marginal itself: unrelated probability mass must not
        # hide a fabricated tiny read. The ulp also handles subnormal sums.
        sum_guard = gamma*marginal+math.ulp(marginal)
        if abs(retained_marginal-marginal) > sum_guard:
            raise ValueError('retained read marginal differs from its complete terminal law')


def verify(path=OUTPUT):
    retained = strict_json.load(Path(path))
    fresh = report()
    compare(retained, fresh)
    # Apply each retained execution's own bound, before tolerant comparison
    # can excuse a scientifically unacceptable stored value.
    for row, current in zip(retained['complete_source_experiments'],
                            fresh['complete_source_experiments'], strict=True):
        if row['complete_output_error'] > float(F(current['numerical_dilation_budget'])):
            raise ValueError('retained complete output exceeds its budget')
        _admit_retained_probabilities(row, current)
    return fresh


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', nargs='?', const=str(OUTPUT))
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        if args.check is not None or args.write:
            parser.error('worker only produces a fresh payload')
        print(json.dumps(_produce(), allow_nan=False))
        return
    if args.check is not None and args.write:
        parser.error('choose write or check')
    result = verify(args.check) if args.check is not None else report()
    if args.write:
        OUTPUT.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps({'schema': result['schema'], 'verified': args.check is not None,
                      'source_experiments': len(result['complete_source_experiments']),
                      'physical_species_identified': False}, indent=2))
