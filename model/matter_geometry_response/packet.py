"""Cold replay of the actual-source electric-energy geometry response."""
import argparse
import ast
from copy import deepcopy
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from tools import strict_json
from . import geometry, source

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')


def canonical(value):
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n'


def source_bindings():
    """Repository import closure plus explicitly loaded pinned carrier files."""
    paths = {str(path.relative_to(ROOT)) for path in Path(__file__).parent.iterdir()
             if path.suffix in ('.py', '.md')}
    paths.update(('model/requirements.txt', 'evidence/native_geometric_source_20260925/spec.json',
        'code/source_selection_model/response.json', 'code/m1_source_realization/receipt.json',
        'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
        'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py'))
    pending = list(paths)

    def add_module(name):
        if not name.startswith(('model.', 'tools.')):
            return
        stem = name.replace('.', '/')
        parts = stem.split('/')
        candidates = [stem+'.py', stem+'/__init__.py']
        candidates.extend('/'.join(parts[:i])+'/__init__.py' for i in range(1, len(parts)))
        for candidate in candidates:
            if (ROOT/candidate).is_file() and candidate not in paths:
                paths.add(candidate)
                pending.append(candidate)

    while pending:
        path = pending.pop()
        if not path.endswith('.py'):
            continue
        module = path.removesuffix('.py').replace('/', '.')
        if not module.startswith(('model.', 'tools.')):
            continue
        package = module.removesuffix('.__init__') if module.endswith('.__init__') else module.rpartition('.')[0]
        add_module(module)
        for node in ast.walk(ast.parse((ROOT/path).read_text())):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    add_module(alias.name)
            elif isinstance(node, ast.ImportFrom):
                target = '.'*node.level+(node.module or '')
                target = importlib.util.resolve_name(target, package) if node.level else target
                add_module(target)
                for alias in node.names:
                    add_module(target+'.'+alias.name)
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(paths)}


def scope():
    return {
        'status': 'SOURCE_ELECTRIC_ENERGY_RESPONSIVE_RANK_THREE_GEOMETRY',
        'same_fixed_instrument_for_all_matter_inputs': True,
        'source_interaction_pointer_and_fair_coin_words_executed': True,
        'unchanged_native_integer_repair_used': True,
        'entire_M78_native_instrument_compiled_into_M6_processor': False,
        'scheduler_table_and_quantum_operation_schedule_are_specified': True,
        'fresh_probe_preparation_and_independence_are_inputs': True,
        'all_terminal_quantum_work_and_classical_records_retained': True,
        'uniform_isolated_rank_three_band_derived': True,
        'finite_energy_dependent_Gram_contrast_certified': True,
        'mean_iteration_for_arbitrary_correlated_probe_reuse': False,
        'whole_common_source_Hamiltonian_evolved_by_interaction_subword': False,
        'dimensionless_response_Gram_identified_with_Lorentzian_metric': False,
        'source_state_map_to_spatial_fermion_core_or_defect_family': False,
        'joint_stress_entropy_area_or_Einstein_premises_instantiated': False,
        'spatial_or_continuum_refinement_proved': False,
        'specified_24000000_attempt_observation_executed': False,
        'physical_gate_accepted': False,
    }


def observation_certificate():
    partial = sum((F(6)**k/math.factorial(k) for k in range(21)), F(0))
    if partial <= 400:
        raise ArithmeticError('finite observation confidence bound failed')
    return {'fresh_attempts_per_energy_endpoint': 12000000,
            'probability_error': '1/2000', 'Gram_error_per_endpoint': '1/100',
            'contrast_lower_on_success': '11/200',
            'failure_probability_strict_upper': '1/100',
            'exp6_rational_lower': str(partial),
            'fresh_fair_coin_words_total': 216000000,
            'large_observation_executed_by_receipt': False}


def _produce():
    pins = source_bindings()
    value = {'schema': 'oph.matter_geometry_response.v1', 'all_pass': True,
             'issues': [1144], 'gates': ['CP-08'],
             'quantum_source': source.quantum_execution(),
             'native_source': source.native_execution(),
             'geometry': geometry.certificate(),
             'finite_observation': observation_certificate(),
             'scope': scope(), 'source_bindings': pins}
    if source_bindings() != pins:
        raise ValueError('source changed during the geometry response replay')
    return value


def report():
    pins = source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-matter-geometry-') as directory:
        temporary = Path(directory)
        output = temporary/'receipt.json'
        subprocess.run([sys.executable, '-W', 'error', '-X',
                        'pycache_prefix='+str(temporary/'pycache'),
                        '-m', 'model.matter_geometry_response.packet', '--worker', str(output)],
                       cwd=ROOT, env=dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1'),
                       check=True, capture_output=True, text=True, timeout=300)
        value = strict_json.load(output)
    if source_bindings() != pins or value['source_bindings'] != pins:
        raise ValueError('source changed across the cold geometry response replay')
    return value


def _comparable(value):
    """Keep exact custody while admitting four bounded numerical diagnostics."""
    comparable = deepcopy(value)
    for resource in ('charged_endpoint_preparation_resources',
                     'interaction_and_scalar_preparation_resources',
                     'one_fair_coin_word_resources', 'pointer_resources'):
        try:
            row = comparable['quantum_source'][resource]
            leakage = row['maximum_numerical_valid_code_leakage']
        except (KeyError, TypeError) as error:
            raise ValueError('missing declared numerical leakage diagnostic: '+resource) from error
        if type(leakage) is not float or not math.isfinite(leakage) or not 0 <= leakage <= 1e-9:
            raise ValueError('numerical leakage exceeds its declared admission: '+resource)
        # Source exponentials and singular-value diagnostics vary with BLAS.
        # Both receipts must pass the existing numerical admission; no other
        # path, exact result, resource count or source binding is normalized.
        row['maximum_numerical_valid_code_leakage'] = 0.0
    return comparable


def validate(value, expected):
    if value.get('scope') != scope():
        raise ValueError('actual-source geometry scope or freshness requirement changed')
    if canonical(_comparable(value)) != canonical(_comparable(expected)):
        raise ValueError('geometry response receipt differs from fresh replay')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', nargs='?', const=OUTPUT, type=Path)
    group.add_argument('--check', nargs='?', const=OUTPUT, type=Path)
    group.add_argument('--worker', type=Path)
    args = parser.parse_args()
    if args.worker:
        args.worker.write_text(canonical(_produce()))
        return
    value = report()
    if args.write:
        args.write.write_text(canonical(value))
    else:
        validate(strict_json.load(args.check), value)
    print('source electric-energy geometry response: PASS; finite rank-three Gram contrast >3/40')


if __name__ == '__main__':
    main()
