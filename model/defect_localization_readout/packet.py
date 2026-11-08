"""Replay complete creation experiments and exact source invariants."""
import ast
import hashlib
import importlib.util
import json
from pathlib import Path

from tools import strict_json
from . import algebra, bounds, execution

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results.json'


def bindings():
    names = {p.relative_to(ROOT).as_posix() for p in HERE.iterdir() if p.suffix in ('.py', '.md')}
    names.update(('model/native_ordered_ground/witness.npz',
                  'model/native_ordered_ground/README.md',
                  'model/protected_spin_limit/README.md',
                  'model/strong_coupling_refinement/README.md',
                  'model/native_spin_observer/README.md', 'model/requirements.txt',
                  'code/source_selection_model/response.json',
                  'code/m1_source_realization/receipt.json',
                  'evidence/native_geometric_source_20260925/spec.json',
                  'model/native_spectral_detector/contract.cpp',
                  'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
                  'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/__init__.py'))
    pending = list(names)
    visited = set()
    legacy = 'evidence/observer_dynamics_20260925/oph-physics-sim/'
    def add(module):
        if not module.startswith(('model.', 'tools.', 'oph_exact.', 'oph_fpe.')):
            return
        prefix = legacy if module.startswith(('oph_exact.', 'oph_fpe.')) else ''
        parts = module.split('.')
        stem = '/'.join(parts)
        paths = [stem+'.py', stem+'/__init__.py']+['/'.join(parts[:i])+'/__init__.py'
                                                for i in range(1, len(parts))]
        for path in paths:
            path = prefix+path
            if (ROOT/path).is_file() and path not in names:
                names.add(path)
                pending.append(path)
    while pending:
        name = pending.pop()
        if name in visited or not name.endswith('.py'):
            continue
        visited.add(name)
        module = name.removeprefix(legacy)[:-3].replace('/', '.')
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
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(names)}


def build():
    before = bindings()
    result = {'schema': 'oph.defect_localization_readout.v1',
              'source_bindings': before,
              'exact_branch_symmetry': algebra.branch_certificate(),
              'actual_cube_excitation': algebra.cube_certificate(),
              'complete_source_execution': execution.run(),
              'full_record_lifetime_bound': bounds.trajectory_example(),
              'scope': {
                  'localized_creation_preserves_initial_exterior': True,
                  'complete_C_E_repairs_preserve_the_symmetry_sector': True,
                  'same_clock_local_signal_has_certified_finite_lifetime': True,
                  'fixed_spatial_support_is_permanently_preserved': False,
                  'topological_defect_or_infinite_volume_superselection_identified': False,
                  'all_admitted_coherent_controls_preserve_parity': False,
                  'bare_repair_records_read_the_internal_spin_defect': False,
                  'coherent_preparation_and_spin_pointer_are_declared_services': True,
                  'relativistic_mass_shell_or_observed_particle_identified': False,
                  'physical_mass_table': []}}
    if bindings() != before:
        raise RuntimeError('source changed during complete defect verification')
    return result


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def verify(path):
    if canonical(strict_json.load(Path(path))) != canonical(build()):
        raise ValueError('defect receipt, source bindings, experiment or scope differs')
