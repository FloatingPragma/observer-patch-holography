"""Reproduce every native covariance column and the protected-load identities."""
import ast
import hashlib
import importlib.util
import json
from pathlib import Path

from tools import strict_json
from . import algebra

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results.json'


def bindings():
    names = {p.relative_to(ROOT).as_posix() for p in HERE.iterdir()
             if p.suffix in ('.py', '.md')}
    names.update(('evidence/native_geometric_source_20260925/spec.json',
                  'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
                  'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py',
                  'model/protected_spin_limit/README.md',
                  'model/native_fair_local_readout/PROOF.md'))
    pending = list(names)
    visited = set()

    def add(module):
        if not module.startswith(('model.', 'tools.')):
            return
        parts = module.split('.')
        stem = '/'.join(parts)
        for name in [stem+'.py', stem+'/__init__.py']+[
                '/'.join(parts[:i])+'/__init__.py' for i in range(1, len(parts))]:
            if (ROOT/name).is_file() and name not in names:
                names.add(name)
                pending.append(name)

    while pending:
        name = pending.pop()
        if name in visited or not name.endswith('.py'):
            continue
        visited.add(name)
        module = name[:-3].replace('/', '.')
        package = module.removesuffix('.__init__') if module.endswith('.__init__') else module.rpartition('.')[0]
        add(module)
        for node in ast.walk(ast.parse((ROOT/name).read_text())):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    add(alias.name)
            elif isinstance(node, ast.ImportFrom):
                target = '.'*node.level+(node.module or '')
                target = importlib.util.resolve_name(target, package) if node.level else target
                add(target)
                for alias in node.names:
                    add(target+'.'+alias.name)
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            for name in sorted(names)}


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def build():
    before = bindings()
    result = {'schema': 'oph.carrier_position_translation.v1',
              'source_bindings': before,
              'native_operator_covariance': algebra.branch_certificate(),
              'complete_adaptive_observer': algebra.complete_observer_certificate(),
              'uniform_native_clock': algebra.mean_certificate(),
              'protected_source_load_readout': algebra.projected_load_certificate(),
              'scope': {'actual_rank_three_carrier_record_displacement': True,
                        'arbitrary_input_operator_covariance_on_branch_support': True,
                        'complete_native_records_retained': True,
                        'full_finite_unitary_nonzero_shift_covariance_possible': False,
                        'carrier_total_load_is_protected_spin_position': False,
                        'field_translation_generator_joined_to_carrier_position': False,
                        'universal_particle_Heisenberg_bound_from_native_records': False,
                        'physical_length_or_hbar_calibration': False}}
    if bindings() != before:
        raise RuntimeError('source changed during covariance verification')
    return result


def verify(path=OUTPUT):
    if canonical(strict_json.load(Path(path))) != canonical(build()):
        raise ValueError('native position covariance evidence or scope differs')
    return True
