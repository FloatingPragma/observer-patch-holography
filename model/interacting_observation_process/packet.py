"""Source-bound cold replay of the interacting multi-time observation comparison."""
import argparse
import ast
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from tools import strict_json
from ..compact_interacting_response import packet as compact_packet
from . import bounds, experiment


ROOT = Path(__file__).resolve().parents[2]
PACKAGE = Path(__file__).resolve().parent
OUTPUT = PACKAGE/'results.json'
PARENT = 'model/compact_interacting_response/results.json'


def source_bindings():
    paths = set(compact_packet.source_bindings()) | {PARENT}
    paths.update(p.relative_to(ROOT).as_posix() for p in PACKAGE.iterdir()
                 if p.suffix in ('.py', '.md'))
    pending, seen = list(paths), set()
    def add(module):
        if not module.startswith(('model.', 'tools.')):
            return
        parts = module.split('.')
        candidates = [module.replace('.', '/')+'.py', module.replace('.', '/')+'/__init__.py']
        candidates += ['/'.join(parts[:i])+'/__init__.py' for i in range(1, len(parts))]
        for path in candidates:
            if (ROOT/path).is_file() and path not in paths:
                paths.add(path)
                pending.append(path)
    while pending:
        path = pending.pop()
        if path in seen or not path.endswith('.py'):
            continue
        seen.add(path)
        module = path[:-3].replace('/', '.')
        package = module.removesuffix('.__init__') if module.endswith('.__init__') else module.rpartition('.')[0]
        add(module)
        for node in ast.walk(ast.parse((ROOT/path).read_text(encoding='utf-8'))):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    add(alias.name)
            elif isinstance(node, ast.ImportFrom):
                target = '.'*node.level+(node.module or '')
                if node.level:
                    target = importlib.util.resolve_name(target, package)
                add(target)
                for alias in node.names:
                    add(target+'.'+alias.name)
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def produce():
    pins = source_bindings()
    parent = strict_json.load(ROOT/PARENT)
    if parent['source_bindings'] != compact_packet.source_bindings():
        raise ValueError('the inherited compact source receipt is stale')
    inherited = parent['uniform_core_limit']
    if (parent['all_pass'] is not True or inherited['localized_annihilators'] != 180
            or inherited['core_Fock_norm_upper'] != str(bounds.COMPACT_NORM)
            or inherited['contraction_upper'] != '1/12'
            or inherited['uniform_in_exterior_volume'] is not True
            or parent['scope']['complete_record_limit_uniform_in_exterior_volume'] is not True
            or parent['rational_ground_and_pole_certificate']['actual_memory_gain_below_one_and_contraction_below_one_over12'] is not True):
        raise ValueError('the source does not supply the claimed common finite inputs')
    result = {'schema': 'oph.interacting_observation_process.v1', 'issues': [1138, 1141, 1142],
              'all_pass': True, 'exact_resource_certificate': bounds.certificate(),
              'complete_record_execution': experiment.certificate(),
              'archive_rate_countercontrol': experiment.archive_countercontrol(),
              'inherited_source': {'receipt': PARENT, 'norm_and_contraction_inputs_checked': True,
                                   'all_predecessor_science_replayed': False},
              'scope': {'same_compact_interacting_memory_generator': True,
                        'arbitrary_reference_adaptive_instrument_comparison': True,
                        'actual_memory_continues_across_interventions': True,
                        'all_native_and_detector_records_retained': True,
                        'explicit_finite_resource_schedule': True,
                        'uniform_in_exterior_volume': True,
                        'uniform_rare_conditional_state_error': False,
                        'all_microscopic_adaptive_policies_converge': False,
                        'native_compilation_of_added_detectors': False,
                        'physical_species_or_masses_derived': False,
                        'interacting_chiral_spatial_continuum': False,
                        'Einstein_register_realized': False,
                        'integrated_gate_accepted': False},
              'canonical_physical_premises_discharged': [], 'physical_comparison_count_added': 0,
              'source_bindings': pins}
    if pins != source_bindings():
        raise ValueError('source changed during observation-process replay')
    return result


def report():
    pins = source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-interacting-observation-') as folder:
        path = Path(folder)/'worker.json'
        env = dict(os.environ, PYTHONUTF8='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
        result = subprocess.run([sys.executable, '-W', 'error', '-X',
                                 'pycache_prefix='+str(Path(folder)/'pycache'),
                                 '-m', 'model.interacting_observation_process.packet', '--worker', str(path)],
                                cwd=ROOT, env=env, text=True, capture_output=True, timeout=300)
        if result.returncode:
            raise RuntimeError('observation-process worker failed:\n'+result.stdout+result.stderr)
        value = strict_json.load(path)
    if pins != source_bindings() or value['source_bindings'] != pins:
        raise ValueError('source changed across cold observation-process replay')
    return value


def validate(value, expected=None):
    expected = report() if expected is None else expected
    if canonical(value) != canonical(expected):
        raise ValueError('interacting observation receipt differs from complete reconstruction')


def verify(path=OUTPUT):
    validate(strict_json.load(Path(path)))
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', nargs='?', const=OUTPUT, type=Path)
    group.add_argument('--check', nargs='?', const=OUTPUT, type=Path)
    group.add_argument('--worker', type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    value = produce() if args.worker else report()
    if args.check:
        validate(strict_json.load(args.check), value)
    else:
        (args.worker or args.write).write_text(canonical(value), encoding='utf-8')
    if not args.worker:
        print('Interacting adaptive observation certificate OK')


if __name__ == '__main__':
    main()
