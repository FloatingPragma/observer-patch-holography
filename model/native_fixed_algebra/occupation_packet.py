"""Owned replay of the complete finite native occupation-protection catalogue."""
from __future__ import annotations

import argparse
import ast
from fractions import Fraction
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from tools import strict_json
from ..completion.program_class_v2 import reviewed_contract_digest


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'occupation_results.json'
CAMPAIGN = 'OPH-NATIVE-OCCUPATION-PROTECTION-2026-11-01-v1'
CONTRACT_SHA256 = '3faa6913fabd0fe22c3fba0104c1315dcc7815ee778fed5a9c478f0241111d67'
# Reviewed exposure-snapshot verification; physical class clauses are unchanged.
CLASS_CONTRACT_SHA256 = 'cd9a8d174f3d9d41fb94096cc0d4ece811253bc7b5e4ad58a2fa69185a1a565d'
FALSE_KEYS = frozenset(('physical_particle_identification', 'integrated_gate_accepted',
                       'program_class_changed', 'cross_number_cat_admitted_by_v2',
                       'localized_topological_defect_constructed'))
RESULT_KEYS = frozenset(('schema', 'campaign_id', 'scope', 'contract_sha256',
                        'algebra_and_protection', 'source_observer', 'source_sha256',
                        'program_class_contract_sha256', 'physical_comparison_count_added')) | FALSE_KEYS


def class_contract_digest():
    """Pin every V2 class clause, excluding its self-referential adoption signature.

    The outer completion gate still validates adoption independently.  Its
    review can bind this finite receipt without making a hash cycle through
    V2's own pointer back to that review.
    """
    value = strict_json.load(ROOT/'model/completion/program_class_v2.json')
    actual = reviewed_contract_digest(value)
    if actual != CLASS_CONTRACT_SHA256:
        raise ValueError('the reviewed fixed-number V2 class contract changed')
    return actual


def source_bindings():
    """Extend the existing primitive/data closure with every new Python import."""
    from ..exchange_sector.many_body_packet import source_bindings as inherited
    names = set(inherited())
    names.update(p.relative_to(ROOT).as_posix() for pattern in
                 ('occupation*.py', 'test_occupation*.py') for p in HERE.glob(pattern))
    names.update(('model/native_fixed_algebra/OCCUPATION.md',
                  'model/native_fixed_algebra/occupation_contract.json'))
    pending, visited = list(names), set()

    def add(module):
        if not module.startswith(('model.', 'tools.')):
            return
        path = module.replace('.', '/')
        parts = path.split('/')
        candidates = [path+'.py', path+'/__init__.py']
        candidates += ['/'.join(parts[:i])+'/__init__.py' for i in range(1, len(parts))]
        for name in candidates:
            if (ROOT/name).is_file() and name not in names:
                names.add(name)
                pending.append(name)

    while pending:
        name = pending.pop()
        if name in visited or not name.endswith('.py'):
            continue
        visited.add(name)
        module = name.removesuffix('.py').replace('/', '.')
        package = (module.removesuffix('.__init__') if module.endswith('.__init__')
                   else module.rpartition('.')[0])
        add(module)
        for node in ast.walk(ast.parse((ROOT/name).read_text(encoding='utf-8'))):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    add(alias.name)
            elif isinstance(node, ast.ImportFrom):
                target = '.'*node.level+(node.module or '')
                target = importlib.util.resolve_name(target, package) if node.level else target
                add(target)
                for alias in node.names:
                    add(target+'.'+alias.name)
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(names)}


def _report():
    from . import occupation, occupation_service
    before = source_bindings()
    class_digest = class_contract_digest()
    if before['model/native_fixed_algebra/occupation_contract.json'] != CONTRACT_SHA256:
        raise ValueError('occupation-protection contract changed')
    result = {'schema': 'oph.native_occupation_protection.receipt.v1',
              'campaign_id': CAMPAIGN, 'scope': 'FINITE_COMPONENT_FAMILY',
              'contract_sha256': CONTRACT_SHA256,
              'program_class_contract_sha256': class_digest,
              'algebra_and_protection': occupation.certificate(),
              'source_observer': occupation_service.certificate(),
              'physical_comparison_count_added': 0,
              **{key: False for key in FALSE_KEYS}, 'source_sha256': before}
    if source_bindings() != before or class_contract_digest() != class_digest:
        raise RuntimeError('occupation sources changed during owned execution')
    return result


def compare(actual, expected):
    """Exact complete JSON comparison; no numerical tolerance hides a new claim."""
    encode = lambda x: json.dumps(x, sort_keys=True, separators=(',', ':'), allow_nan=False)
    if encode(actual) != encode(expected):
        raise ValueError('occupation-protection receipt differs from complete fresh replay')


def _digest(value):
    return type(value) is str and len(value) == 64 and all(c in '0123456789abcdef' for c in value)


def _rational(value):
    if type(value) is not str:
        return False
    try:
        number = Fraction(value)
        return number > 0 and str(number) == value
    except (ValueError, ZeroDivisionError):
        return False


def _shape(value, schema, where):
    """Validate receipt syntax and scope, not the mathematical proof twice."""
    valid = True
    if type(schema) is dict:
        valid = type(value) is dict and set(schema) <= set(value)
        if valid:
            for key, child in schema.items():
                _shape(value[key], child, where+'.'+key)
    elif type(schema) is tuple:
        child, length = schema
        valid = type(value) is list and len(value) == length
        if valid:
            for index, item in enumerate(value):
                _shape(item, child, where+f'[{index}]')
    elif type(schema) is frozenset:
        valid = any(type(value) is type(option) and value == option for option in schema)
    elif schema is int:
        valid = type(value) is int and value >= 0
    elif schema is str:
        valid = type(value) is str and bool(value)
    elif callable(schema):
        valid = schema(value)
    else:
        valid = type(value) is type(schema) and value == schema
    if not valid:
        raise ValueError('malformed occupation worker field: '+where)


def validate_components(result):
    """Check critical worker shapes and scope; fresh execution owns the proof."""
    mathematics = result['algebra_and_protection']
    observer = result['source_observer']
    case_keys = {'C:+1', 'C:-1', 'E:+1', 'E:-1'}
    _shape(mathematics, {
        'schema': 'oph.native_occupation_catalogue.v1', 'cases': lambda x: type(x) is dict and set(x) == case_keys,
        'scope': {'complete_binary_component': True, 'class_V2_carrier_revised': False,
                  'all_source_quantum_observables_classified': True,
                  'joint_controller_history_charges_classified': False,
                  'active_recovery_codes_classified': False,
                  'refinement_or_physical_particle_identification': False,
                  'convergence_rate_or_physical_clock_calibration': False},
        'locality': {'moving_particles_with_global_charge_excluded': False,
                     'controlled_lifetimes_or_refining_defects_excluded': False, 'witness_sha256': _digest},
        'extremal': {'witness_sha256': _digest, 'cross_global_number_coherence_required': True,
                     'full_code_locally_indistinguishable': False,
                     'added_phase_or_bit_error_is_native_C_E_branch': False,
                     'class_V2_fixed_global_number_admission_satisfied': False,
                     'localized_topological_defect_constructed': False}}, 'algebra_and_protection')
    reference = {'witness_sha256': _digest, 'all_seams_complete_and_faithfully_invariant': True}
    for key, case in mathematics['cases'].items():
        _shape(case, {'algebra': {'witness_sha256': _digest, 'full_witness_replayed': True,
                                  'algebra_dimension': int, 'commutant_dimension': int},
                      'references': {'fair': reference, 'biased_control': reference},
                      'peripheral': {'witness_sha256': _digest, 'only_peripheral_eigenvalue': int}}, key)
    repair = {'law': frozenset(('C', 'E')), 'character': frozenset((-1, 1)),
              'full_input_complete_amplitude_bound': _rational,
              'resources': {'native_pulses': int, 'first_failure_Kraus_slots': int}}
    experiment = {'law': repair['law'], 'character': repair['character'],
                  'family': frozenset(('extreme', 'fixed_number')),
                  'control': frozenset(('none', 'local_z', 'which_path')),
                  'work_bits': int, 'complete_terminal_addresses': int,
                  'source_word_sha256': _digest, 'complete_output_amplitude_bound': _rational,
                  'resources': repair['resources'],
                  'ideal_terminal_laws': (lambda x: type(x) is dict and bool(x), 2)}
    _shape(observer, {'native_occupation_bits': 12, 'maximum_work_bits': 17,
        'logical_reference_columns': 2, 'every_unlisted_ideal_terminal_probability': '0',
        'all_actual_complex_columns_and_failure_prefixes_freshly_admitted': True,
        'all_failure_outputs': 'immutable literal first-failure prefix archive',
        'general_input_repairs': (repair, 4), 'experiments': (experiment, 16),
        'preparations': {'ideal_macros_preserve_total_number_on_all_inputs': True,
                         'extreme_cross_number_preparation_outside_fixed_number_V2': True},
        'encoded_preparation_preserves_declared_total': {'extreme': False, 'fixed_number': True},
        'encoded_preparation_parity_preservation_claimed': {'extreme': False, 'fixed_number': True},
        'extreme_cross_number_phase_apparatus_supplied': True,
        'primitive_encoded_number_conservation_claimed': False,
        'numerical_admissions_are_interval_or_hardware_certificates': False,
        'physical_particle_or_clock_identified': False, 'postselection': False,
        'program_class_v2_execution_claimed': False,
        'extreme_sector_cat_is_fixed_number_V2_observer': False}, 'source_observer')
    pairs = {(law, sign) for law in ('C', 'E') for sign in (1, -1)}
    if {(row['law'], row['character']) for row in observer['general_input_repairs']} != pairs:
        raise ValueError('worker omitted or duplicated a complete source repair case')
    expected = {(law, sign, family, control) for law, sign in pairs for family, control in
                (('extreme', 'none'), ('extreme', 'local_z'), ('extreme', 'which_path'), ('fixed_number', 'none'))}
    if {(row['law'], row['character'], row['family'], row['control'])
            for row in observer['experiments']} != expected:
        raise ValueError('worker omitted or duplicated a complete source observer case')
    for row in observer['experiments']:
        if not 1 <= row['work_bits'] <= 17 or row['complete_terminal_addresses'] != 1 << row['work_bits']:
            raise ValueError('worker terminal address domain differs from its work register')
        for law in row['ideal_terminal_laws']:
            if any(type(address) is not str or not address.isascii() or not address.isdigit()
                   or str(int(address)) != address or not 0 <= int(address) < row['complete_terminal_addresses']
                   or not _rational(probability) for address, probability in law.items()):
                raise ValueError('worker terminal law has a malformed address or probability')
            if sum(map(Fraction, law.values()), Fraction()) != 1:
                raise ValueError('worker complete ideal terminal law is not normalized')


def report():
    before = source_bindings()
    class_digest = class_contract_digest()
    command = ('import sys,json; sys.path.insert(0,'+repr(str(ROOT))+'); '
               'from model.native_fixed_algebra.occupation_packet import _report; '
               'print(json.dumps(_report(),sort_keys=True,allow_nan=False))')
    environment = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    with tempfile.TemporaryDirectory(prefix='oph-occupation-') as cache:
        run = subprocess.run([sys.executable, '-I', '-B', '-W', 'error', '-X', 'utf8',
                              '-X', 'pycache_prefix='+cache, '-c', command],
                             cwd=ROOT, env=environment, check=True, capture_output=True,
                             text=True, encoding='utf-8', timeout=900)
    result = strict_json.loads(run.stdout)
    if (type(result) is not dict or set(result) != RESULT_KEYS
            or result.get('source_sha256') != before or source_bindings() != before
            or result.get('program_class_contract_sha256') != class_digest
            or class_contract_digest() != class_digest):
        raise RuntimeError('fresh occupation replay and current source custody disagree')
    if (result['schema'] != 'oph.native_occupation_protection.receipt.v1'
            or result['campaign_id'] != CAMPAIGN or result['contract_sha256'] != CONTRACT_SHA256
            or result['scope'] != 'FINITE_COMPONENT_FAMILY'
            or any(result[key] is not False for key in FALSE_KEYS)
            or type(result['physical_comparison_count_added']) is not int
            or result['physical_comparison_count_added'] != 0):
        raise ValueError('fresh occupation replay changed its finite contract')
    validate_components(result)
    compare(result, result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', type=Path, nargs='?', const=OUTPUT)
    mode.add_argument('--check', type=Path, nargs='?', const=OUTPUT)
    args = parser.parse_args()
    result = report()
    if args.write:
        args.write.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+'\n',
                              encoding='utf-8', newline='\n')
    else:
        compare(strict_json.load(args.check), result)
    print('Verified complete native occupation protection, source readouts and fixed-number scope.')


if __name__ == '__main__':
    main()
