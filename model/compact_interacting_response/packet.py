"""Replay the compact source, exact interaction spectrum and record bounds."""
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

from model.compact_defect_pair import packet as compact_packet
from model.interacting_defect_response import packet as radial_packet
from tools import strict_json
from . import algebra, finite

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = Path(__file__).resolve().parent
OUTPUT = PACKAGE/'results.json'


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def source_bindings():
    paths = set(compact_packet.source_bindings()) | radial_packet.source_paths()
    paths.update(str(path.relative_to(ROOT)) for path in PACKAGE.iterdir()
                 if path.is_file() and path.suffix in ('.py', '.md'))
    pending, visited = list(paths), set()

    def add_module(name):
        if not name.startswith(('model.', 'tools.')):
            return
        stem = name.replace('.', '/')
        parts = stem.split('/')
        candidates = [stem+'.py', stem+'/__init__.py']
        candidates += ['/'.join(parts[:end])+'/__init__.py' for end in range(1, len(parts))]
        for candidate in candidates:
            if (ROOT/candidate).is_file() and candidate not in paths:
                paths.add(candidate)
                pending.append(candidate)

    while pending:
        path = pending.pop()
        if path in visited or not path.endswith('.py'):
            continue
        visited.add(path)
        module = path.removesuffix('.py').replace('/', '.')
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
    if any(path.startswith('docs/') for path in paths):
        raise ValueError('research notes are not executable source evidence')
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(paths)}


def scope():
    return {'same_compact_finite_energy_texture_and_nonlinear_core': True,
            'complete_native_memory_mass_operator_retained': True,
            'all_twelve_uniform_native_currents_act_on_same_core': True,
            'all_twelve_currents_are_conserved_by_fixed_Higgs_Hamiltonian': False,
            'nonzero_compact_kinetic_splitting_retained': True,
            'source_lepton_slice_ground_and_poles_interval_certified': True,
            'complete_record_limit_uniform_in_exterior_volume': True,
            'same_interaction_is_strictly_site_local': False,
            'four_mode_slice_is_full_core_vacuum': False,
            'four_mode_slice_reduces_each_individual_record': False,
            'autonomous_texture_creation_or_transport_proved': False,
            'physical_rest_masses_or_observed_species_derived': False,
            'interacting_chiral_relativistic_continuum_derived': False,
            'interacting_electromagnetic_pixel_feedback_closed': False,
            'conserved_relativistic_stress_or_Einstein_response_derived': False,
            'integrated_gate_accepted': False}


def produce():
    pins = source_bindings()
    value = {'schema': 'oph.compact_interacting_response.v1', 'all_pass': True,
             'issues': [1140, 1141, 1142],
             'exact_algebra': algebra.exact_identities(),
             'actual_source_controls': finite.source_controls(),
             'rational_ground_and_pole_certificate': algebra.ground_certificate(),
             'paired_ground_control': finite.block_readout(1., .5, .1, .2),
             'complete_record_probe': finite.complete_memory_control(),
             'uniform_core_limit': algebra.compact_memory_bound(1, algebra.F(1, 128), algebra.F(1, 12), 1, 2**30),
             'finite_energy_difference': {
                 'uniform_reference': 'same pairing/memory, uniform north texture, no localized quartic',
                 'fermion_Hamiltonian_difference_norm_upper': '90*M*((2*a+3)^3-1)+8100*lambda',
                 'bound_at_M64_a1_lambda1_over128': str(algebra.F(90*64*(125-1))+algebra.F(8100, 128)),
                 'uniform_in_exterior_volume': True,
                 'constituent_rest_mass_interpretation': False},
             'scope': scope(), 'canonical_physical_premises_discharged': [],
             'physical_comparison_count_added': 0, 'source_bindings': pins}
    if pins != source_bindings():
        raise ValueError('source changed during compact interacting calculation')
    return value


def report():
    pins = source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-compact-interaction-') as directory:
        path = Path(directory)/'worker.json'
        result = subprocess.run([sys.executable, '-W', 'error', '-m',
                                 'model.compact_interacting_response.packet', '--worker', str(path)],
                                cwd=ROOT, env=dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1'),
                                text=True, capture_output=True, timeout=600)
        if result.returncode:
            raise RuntimeError('compact interaction worker failed:\n'+result.stdout+result.stderr)
        value = strict_json.load(path)
    if value['source_bindings'] != pins or source_bindings() != pins:
        raise ValueError('source changed across compact interacting replay')
    return value


def validate(value, replay=None):
    if canonical(value) != canonical(report() if replay is None else replay):
        raise ValueError('compact interacting receipt differs from complete replay')
    if value.get('source_bindings') != source_bindings():
        raise ValueError('compact interacting source manifest is stale')


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
        print('wrote compact interacting source and spectral readouts')
    else:
        validate(strict_json.load(Path(args.check)))
        print('compact source, nonlinear spectrum and complete records verified')


if __name__ == '__main__':
    main()
