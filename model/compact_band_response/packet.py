"""Source-bound replay of compact-band work and thermal response."""
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
from . import response

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = Path(__file__).resolve().parent
OUTPUT = PACKAGE/'results.json'


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def source_bindings():
    from model.compact_defect_pair.packet import source_bindings as compact_bindings
    from model.interacting_defect_response.packet import source_bindings as radial_bindings
    names = set(compact_bindings()) | set(radial_bindings())
    names.update(str(p.relative_to(ROOT)) for p in PACKAGE.iterdir()
                 if p.is_file() and p.suffix in ('.py', '.md'))
    names.add('model/compact_interacting_response/PROOF.md')
    pending, visited = list(names), set()

    def add_module(name):
        if not name.startswith(('model.', 'tools.')):
            return
        parts = name.split('.')
        candidates = [name.replace('.', '/')+'.py', name.replace('.', '/')+'/__init__.py']
        candidates += ['/'.join(parts[:i])+'/__init__.py' for i in range(1, len(parts))]
        for path in candidates:
            if (ROOT/path).is_file() and path not in names:
                names.add(path); pending.append(path)

    while pending:
        path = pending.pop()
        if path in visited or not path.endswith('.py'):
            continue
        visited.add(path)
        module = path[:-3].replace('/', '.')
        package = module.removesuffix('.__init__') if module.endswith('.__init__') else module.rpartition('.')[0]
        add_module(module)
        for node in ast.walk(ast.parse((ROOT/path).read_text())):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    add_module(alias.name)
            elif isinstance(node, ast.ImportFrom):
                target = '.'*node.level+(node.module or '')
                if node.level:
                    target = importlib.util.resolve_name(target, package)
                add_module(target)
                for alias in node.names:
                    add_module(target+'.'+alias.name)
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(names)}


def scope():
    return {'same_compact_source_hopping_and_interaction_force': True,
            'moving_localization_projector_included': True,
            'all_complete_memory_record_Hamiltonians_have_same_drive_force': True,
            'finite_Kato_pairing_and_native_current_compatibility': True,
            'source_spectral_complement_retained': True,
            'rational_volume_uniform_derivative_bounds': True,
            'thermal_readout_uses_actual_compact_band_splitting': True,
            'Kato_transporter_is_autonomous_native_evolution': False,
            'hopping_drive_identified_with_metric_variation': False,
            'full_core_Gibbs_state_executed': False,
            'conserved_relativistic_stress_constructed': False,
            'regional_area_entropy_or_Einstein_stationarity_derived': False,
            'integrated_gate_accepted': False}


def produce():
    pins = source_bindings()
    result = {'schema': 'oph.compact_band_response.v1', 'all_pass': True,
              'issues': [1140, 1145],
              'moving_compact_source': response.source_control(),
              'exact_missing_force_control': response.force_control(),
              'same_source_thermal_response': response.thermal_control(),
              'scope': scope(), 'canonical_physical_premises_discharged': [],
              'source_bindings': pins}
    if pins != source_bindings():
        raise ValueError('source changed during response calculation')
    return result


def report():
    pins = source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-compact-work-') as temporary:
        path = Path(temporary)/'worker.json'
        child = subprocess.run([sys.executable, '-W', 'error', '-m',
                                'model.compact_band_response.packet', '--worker', str(path)],
                               cwd=ROOT, env=dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1'),
                               capture_output=True, text=True, timeout=300)
        if child.returncode:
            raise RuntimeError('compact-work replay failed:\n'+child.stdout+child.stderr)
        result = strict_json.load(path)
    if result['source_bindings'] != pins or source_bindings() != pins:
        raise ValueError('source changed across response replay')
    return result


def validate(value, replay=None):
    if canonical(value) != canonical(report() if replay is None else replay):
        raise ValueError('compact-work receipt differs from complete replay')
    if value.get('source_bindings') != source_bindings():
        raise ValueError('compact-work source manifest is stale')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', nargs='?', const=str(OUTPUT))
    group.add_argument('--check', nargs='?', const=str(OUTPUT))
    group.add_argument('--worker')
    args = parser.parse_args()
    if args.worker:
        Path(args.worker).write_text(canonical(produce()))
    elif args.write:
        Path(args.write).write_text(canonical(report()))
        print('wrote compact-band work and thermal response')
    else:
        validate(strict_json.load(Path(args.check)))
        print('compact-band work and thermal response verified')


if __name__ == '__main__':
    main()
