"""Source-bound replay of joint positive spectral bins and native outputs."""
import argparse
import ast
import copy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from fractions import Fraction as F
import numpy as np
from tools import strict_json
from . import algebra, spectral, experiment, bounds

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).parent
OUTPUT = HERE / 'results.json'


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + '\n'


def compare(expected, actual):
    if type(expected) is not type(actual):
        raise ValueError('diagnostic type changed')
    if isinstance(expected, dict):
        if expected.keys() != actual.keys():
            raise ValueError('diagnostic fields changed')
        for key in expected:
            compare(expected[key], actual[key])
    elif isinstance(expected, list):
        if len(expected) != len(actual):
            raise ValueError('diagnostic dimension changed')
        for a, b in zip(expected, actual, strict=True):
            compare(a, b)
    elif type(expected) is float:
        if not math.isfinite(actual) or abs(expected - actual) > max(1e-40, min(2e-11, abs(expected) * 1e-7)):
            raise ValueError('native numerical diagnostic changed')
    elif expected != actual:
        raise ValueError('exact diagnostic value changed')


def bindings():
    paths = set(HERE.glob('*.py'))
    queue = list(paths)
    legacy = ROOT / 'evidence/observer_dynamics_20260925/oph-physics-sim'
    def include(path):
        if path.is_file() and path not in paths:
            paths.add(path)
            queue.append(path)
    include(ROOT / 'model/__init__.py')
    include(legacy / 'oph_exact/carrier.py')
    include(legacy / 'oph_exact/__init__.py')
    while queue:
        path = queue.pop()
        base = legacy if path.is_relative_to(legacy) else ROOT
        package = '.'.join(path.relative_to(base).with_suffix('').parts[:-1])
        for node in ast.walk(ast.parse(path.read_text())):
            names = []
            if isinstance(node, ast.Import):
                names = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom):
                name = '.' * node.level + (node.module or '')
                module = importlib.util.resolve_name(name, package) if node.level else name
                names = [module] + [module + '.' + a.name for a in node.names]
            for name in names:
                root = legacy if name.split('.')[0] in ('oph_exact', 'oph_fpe') else ROOT
                pieces = name.split('.')
                include(root.joinpath(*pieces).with_suffix('.py'))
                for j in range(1, len(pieces) + 1):
                    include(root.joinpath(*pieces[:j]) / '__init__.py')
    paths.update(ROOT / p for p in (
        'code/source_selection_model/response.json',
        'evidence/native_geometric_source_20260925/spec.json',
        'code/m1_source_realization/receipt.json',
        'model/native_spectral_detector/contract.cpp',
        'model/native_ordered_ground/witness.npz',
        'model/native_ordered_ground/README.md',
        'model/native_spin_observer/README.md',
        'model/native_joint_spin_spectroscopy/README.md',
        'model/requirements.txt',
        'paper/tex_fragments/NATIVE_JOINT_SPIN_SPECTROSCOPY.tex'))
    return {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}


def theorem():
    return {'finite_attachment': 'FA102', 'positive_complete_joint_instrument': True,
            'cube_labels': 'combined graph-reflection characters',
            'periodic_torus_translation_theorem_separate': True,
            'native_dimer_executed': True, 'native_full_cube_executed': False,
            'packet_success_weight_retained': True,
            'ground_specific_dimer_split_error_zero': True,
            'arbitrary_input_split_error_zero': False,
            'sharp_scaling_shell_proved': False, 'physical_mass_identified': False,
            'infinite_ground_marginal_modulus_supplied': False,
            'old_broad_receipts_replayed': False,
            'actual_executing_sources_and_control_data_bound': True}


def produce():
    before = bindings()
    cube = copy.deepcopy(spectral.cube_readout())
    execution, arrays = experiment.run()
    if before != bindings():
        raise ValueError('source changed during joint spectroscopy replay')
    return {'schema': 'oph.native_joint_spin_spectroscopy.v1', 'new_checks_pass': True,
            'theorem': theorem(), 'source_bindings': before, 'cube_joint_bins': cube,
            'resource_bounds': bounds.certificate(),
            'dimer_common_kernel': algebra.dimer_common_kernel(), 'execution': execution}, arrays


def validate(value, arrays):
    required = {'schema', 'new_checks_pass', 'theorem', 'source_bindings', 'cube_joint_bins', 'dimer_common_kernel', 'execution', 'resource_bounds'}
    if set(value) != required or value['schema'] != 'oph.native_joint_spin_spectroscopy.v1' or value['new_checks_pass'] is not True:
        raise ValueError('invalid complete joint spectroscopy packet')
    for key, expected in (('theorem', theorem()), ('source_bindings', bindings()), ('dimer_common_kernel', algebra.dimer_common_kernel()), ('resource_bounds', bounds.certificate())):
        if canonical(value[key]) != canonical(expected):
            raise ValueError('exact source/proof/type mismatch: ' + key)
    if canonical(value['cube_joint_bins']) != canonical(spectral.cube_readout()):
        raise ValueError('exact full cube probability certificate changed')
    if set(arrays) != {'success', 'failure_weights'}:
        raise ValueError('all native output arrays required')
    word, initial, plan = experiment.plan()
    success, failure = arrays['success'], arrays['failure_weights']
    if success.dtype != np.complex128 or success.shape != (128, 1) or failure.dtype != np.float64 or failure.shape != (len(word), 4):
        raise ValueError('full native output carrier changed')
    if not np.isfinite(success).all() or not np.isfinite(failure).all() or np.any(failure < 0):
        raise ValueError('invalid native output')
    pairs = np.array([event['kind'] == 'J' for event in word])
    if np.any(failure[pairs, 2:] != 0):
        raise ValueError('nonphysical padded native outcome')
    target, probabilities = experiment.ideal_target()
    distance = float(np.sqrt(np.linalg.norm(success - target) ** 2 + failure.sum()))
    normalization = abs(float(np.vdot(success, success).real + failure.sum()) - 1)
    if distance > float(F(plan['complete_isometry_error'])) + 2e-8 or normalization > 2e-8:
        raise ValueError('native output exceeds its certificate')
    observed = np.zeros((2, 4))
    for k in range(2):
        for j in range(4):
            start = (k << 4) + (j << 5)
            observed[k, j] = float(np.linalg.norm(success[start:start + 16]) ** 2)
    count = value['execution'].get('rounded_negative_count')
    mass = value['execution'].get('rounded_negative_mass')
    if type(count) is not int or not 0 <= count <= 4 * len(word) or type(mass) is not float or not 0 <= mass < 1e-35:
        raise ValueError('invalid negative-roundoff diagnostic')
    expected = dict(plan, ideal_joint_probabilities=probabilities.tolist(),
                    native_success_joint_probabilities=observed.tolist(),
                    complete_output_distance_diagnostic=distance,
                    full_first_failure_probability=float(failure.sum()),
                    normalization_residual=normalization,
                    rounded_negative_count=count, rounded_negative_mass=mass,
                    native_failure_slots=4 * len(word) - 2 * int(pairs.sum()))
    compare(expected, value['execution'])


def write(path, value, arrays):
    archive = path.with_suffix('.npz')
    np.savez_compressed(archive, **arrays)
    entry = {'filename': archive.name, 'sha256': hashlib.sha256(archive.read_bytes()).hexdigest()}
    path.write_text(canonical(dict(value, complete_record_archive=entry)))


def read(path):
    value = strict_json.load(path)
    entry = value.pop('complete_record_archive')
    archive = path.with_suffix('.npz')
    if set(entry) != {'filename', 'sha256'} or entry['filename'] != archive.name or entry['sha256'] != hashlib.sha256(archive.read_bytes()).hexdigest():
        raise ValueError('paired full native archive changed')
    with np.load(archive, allow_pickle=False) as loaded:
        arrays = {key: loaded[key] for key in loaded.files}
    validate(value, arrays)
    return value, arrays


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', type=Path, nargs='?', const=OUTPUT)
    group.add_argument('--check', type=Path, nargs='?', const=OUTPUT)
    group.add_argument('--worker', type=Path)
    args = parser.parse_args()
    if args.worker:
        value, arrays = produce()
        validate(value, arrays)
        write(args.worker, value, arrays)
        return
    before = bindings()
    with tempfile.TemporaryDirectory(prefix='oph-joint-spectroscopy-') as temp:
        candidate = Path(temp) / 'worker.json'
        env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
        subprocess.run([sys.executable, '-W', 'error', '-m', 'model.native_joint_spin_spectroscopy', '--worker', str(candidate)], cwd=ROOT, env=env, check=True)
        value, arrays = read(candidate)
    if before != bindings():
        raise ValueError('source changed during cold spectroscopy replay')
    if args.check:
        old, retained = read(args.check)
        for key in ('schema', 'new_checks_pass', 'theorem', 'source_bindings', 'cube_joint_bins', 'dimer_common_kernel', 'resource_bounds'):
            if canonical(old[key]) != canonical(value[key]):
                raise ValueError('cold exact replay changed: ' + key)
        compare(old['execution'], value['execution'])
        for key in arrays:
            atol, rtol = (1e-40, 1e-7) if key == 'failure_weights' else (2e-11, 0)
            if not np.allclose(arrays[key], retained[key], atol=atol, rtol=rtol):
                raise ValueError('cold complete native array changed')
    else:
        write(args.write, value, arrays)
    print('Joint spectroscopy source, positive bins and full native replay PASS')
