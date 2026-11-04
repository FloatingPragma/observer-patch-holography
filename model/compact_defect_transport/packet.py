"""Replay the complete guarded law and every finite transport control."""
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
from . import bounds, source, transport

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = Path(__file__).resolve().parent
OUTPUT = PACKAGE/'results.json'


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def source_bindings():
    from model.common_source.chiral_candidate.action import SOURCE_FILES
    names = set(SOURCE_FILES) | {'model/requirements.txt', 'requirements.txt',
                                'model/compact_defect_pair/PROOF.md'}
    names.update(str(p.relative_to(ROOT)) for p in PACKAGE.iterdir()
                 if p.is_file() and p.suffix in ('.py', '.md'))
    pending, visited = list(names), set()

    def add_module(name):
        if not name.startswith(('model.', 'tools.')):
            return
        stem = name.replace('.', '/')
        parts = stem.split('/')
        candidates = [stem+'.py', stem+'/__init__.py']
        candidates += ['/'.join(parts[:end])+'/__init__.py' for end in range(1, len(parts))]
        for candidate in candidates:
            if (ROOT/candidate).is_file() and candidate not in names:
                names.add(candidate)
                pending.append(candidate)

    while pending:
        name = pending.pop()
        if name in visited or not name.endswith('.py'):
            continue
        visited.add(name)
        module = name.removesuffix('.py').replace('/', '.')
        package = module.removesuffix('.__init__') if module.endswith('.__init__') else module.rpartition('.')[0]
        add_module(module)
        for node in ast.walk(ast.parse((ROOT/name).read_text())):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    add_module(alias.name)
            elif isinstance(node, ast.ImportFrom):
                target = '.'*node.level+(node.module or '')
                target = importlib.util.resolve_name(target, package) if node.level else target
                add_module(target)
                for alias in node.names:
                    add_module(target+'.'+alias.name)
    if any(name.startswith('docs/') for name in names):
        raise ValueError('research notes are not source evidence')
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(names)}


def scope():
    return {'chosen_guarded_local_texture_repair_law': True,
            'all_continuing_proposal_and_rejection_records_retained': True,
            'exact_shell_degree_protection_on_declared_law': True,
            'uniform_four_mode_gap_and_quench_leakage_bound': True,
            'polar_band_transport_preserves_pair_tensor_and_uniform_native_currents': True,
            'same_band_number_interaction_transported': True,
            'complete_finite_local_quadratic_gate_factorization': True,
            'band_protection_applies_at_record_endpoints': True,
            'band_invariance_during_each_control_pulse': False,
            'pulse_timing_and_drift_cancellation_are_supplied': True,
            'unchanged_native_CE_pulse_compilation': False,
            'autonomous_texture_formation_or_moving_core_dynamics': False,
            'local_Gauss_dressing_or_physical_species_identification': False,
            'relativistic_dispersion_or_interacting_continuum': False,
            'conserved_gravitational_stress': False,
            'integrated_gate_accepted': False}


def produce():
    pins = source_bindings()
    value = {'schema': 'oph.compact_defect_transport.v1', 'all_pass': True,
             'issues': [1139, 1141], 'exact_complete_repair': source.exact_control(),
             'rational_transport_bounds': [bounds.certificate(rounds=n) for n in (1, 4, 16)],
             'complete_source_transport_controls': transport.control(),
             'scope': scope(), 'canonical_physical_premises_discharged': [],
             'physical_comparison_count_added': 0, 'source_bindings': pins}
    if pins != source_bindings():
        raise ValueError('source changed during calculation')
    return value


def report():
    pins = source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-compact-transport-') as directory:
        output = Path(directory)/'worker.json'
        child = subprocess.run([sys.executable, '-W', 'error', '-m',
                                'model.compact_defect_transport.packet', '--worker', str(output)],
                               cwd=ROOT, env=dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1'),
                               capture_output=True, text=True, timeout=300, check=False)
        if child.returncode:
            raise RuntimeError('compact transport worker failed:\n'+child.stdout+child.stderr)
        value = strict_json.load(output)
    if value['source_bindings'] != pins or source_bindings() != pins:
        raise ValueError('source changed across cold replay')
    return value


def validate(value, replay=None):
    if canonical(value) != canonical(report() if replay is None else replay):
        raise ValueError('compact transport receipt differs from complete cold replay')
    if value['source_bindings'] != source_bindings():
        raise ValueError('source manifest is stale or incomplete')


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
        print('wrote complete guarded compact-defect transport')
    else:
        validate(strict_json.load(Path(args.check)))
        print('compact topology, complete repairs, band transport and local gates verified')


if __name__ == '__main__':
    main()
