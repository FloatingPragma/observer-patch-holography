"""Recompute the quartic defect source and all independent finite controls."""
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
from . import algebra,finite,thermodynamics

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = Path(__file__).resolve().parent
OUTPUT = PACKAGE/'results.json'
DATA = (
    'code/particles/hierarchy/certificates/R_P_source_audit_pixel_certificate.json',
    'code/particles/hierarchy/certificates/R_EW_tick_projection_certificate.json',
    'code/particles/leptons/derive_charged_mcpr_completion_conditional.py',
    'model/observer_particles/results/representation_identification_operators.json',
    'model/observer_particles/results/representation_identification.json',
    'model/observer_particles/results/unified_yukawa_operators.json',
    'model/observer_particles/results/unified_yukawa_action.json',
    'model/radial_defect_gap/PROOF.md',
    'model/defect_current_bridge/PROOF.md',
    'model/mcpr_native_memory_bridge/PROOF.md',
    'evidence/native_geometric_source_20260925/spec.json',
    'model/requirements.txt',
)


def canonical(value):
    return json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n'


def source_paths():
    """Python import closure and explicit consumed data; no research-note pins."""
    paths = set(DATA)
    archive = ROOT/'evidence/observer_dynamics_20260925/oph-physics-sim'
    paths.update(str(path.relative_to(ROOT)) for path in archive.rglob('*.py'))
    paths.update(str(path.relative_to(ROOT)) for path in PACKAGE.iterdir()
                 if path.is_file() and path.suffix in ('.py','.md'))
    pending = list(paths)
    visited = set()
    def add(module):
        if not module.startswith(('model.','tools.')):
            return
        parts = module.split('.')
        candidates = [module.replace('.','/')+'.py',module.replace('.','/')+'/__init__.py']
        candidates += ['/'.join(parts[:i])+'/__init__.py' for i in range(1,len(parts))]
        for path in candidates:
            if (ROOT/path).is_file() and path not in paths:
                paths.add(path)
                pending.append(path)
    while pending:
        path = pending.pop()
        if path in visited or not path.endswith('.py'):
            continue
        visited.add(path)
        module = path[:-3].replace('/','.')
        package = module.removesuffix('.__init__') if module.endswith('.__init__') else module.rpartition('.')[0]
        add(module)
        for node in ast.walk(ast.parse((ROOT/path).read_text())):
            if isinstance(node,ast.Import):
                for alias in node.names:
                    add(alias.name)
            elif isinstance(node,ast.ImportFrom):
                target = '.'*node.level+(node.module or '')
                if node.level:
                    target = importlib.util.resolve_name(target,package)
                add(target)
                for alias in node.names:
                    add(target+'.'+alias.name)
    return paths


def source_bindings():
    return {path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(source_paths())}


def scope():
    return {'same_native90_mode_core_and_memory_mass_operator':True,
            'genuine_nonquadratic_core_interaction':True,
            'all12_native_currents_commute_with_interaction':True,
            'complete_record_core_limit_uniform_in_cube_radius':True,
            'full_history_core_comparison_rate': 'O(n^-1/2)',
            'unread_core_comparison_rate': 'O(n^-1)',
            'initial_memory_correlations_retained_in_ideal_history_law':True,
            'rare_conditioned_history_error_uniform':False,
            'intervening_field_measurements_covered':False,
            'finite_faithful_Gibbs_reference_and_noncommuting_BKM_response':True,
            'interaction_strictly_site_local':False,
            'interaction_exponentially_localized_about_prescribed_origin':True,
            'native_CE_prepares_Gibbs_reference':False,
            'invariant_four_mode_slice_is_full_core_vacuum':False,
            'unchanged_MCPR_coefficients_are_interacting_physical_masses':False,
            'local_Gauss_constraint_or_chiral_relativistic_continuum':False,
            'G1_through_G6_joint_geometric_premises_realized':False,
            'Einstein_geometry_derived':False,
            'integrated_gate_accepted':False}


def produce():
    pins = source_bindings()
    result = {'schema':'oph.interacting_defect_response.v1','all_pass':True,
              'issues':[1140,1142,1145],
              'exact_identities':algebra.exact_certificate(),
              'onsite_negative_control':finite.onsite_control(),
              'source_blocks':finite.source_controls(),
              'independent_nongeneric_control':finite.spectrum_control(1.,.5),
              'finite_Gibbs_response':thermodynamics.control(),
              'complete_record_nonlinear_control':finite.memory_control(),
              'radius_uniform_core_limit':algebra.radius_uniform_bound(1,1,algebra.F(1,12),1,1024),
              'localization_bounds':[algebra.localization_bound(i) for i in (1,2,4,8)],
              'scope':scope(),'canonical_physical_premises_discharged':[],
              'source_bindings':pins}
    if pins != source_bindings():
        raise ValueError('source changed during calculation')
    return result


def report():
    pins = source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-interacting-defect-') as directory:
        path = Path(directory)/'worker.json'
        environment = dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
        child = subprocess.run([sys.executable,'-W','error','-m',
                                'model.interacting_defect_response.packet','--worker',str(path)],
                               cwd=ROOT,env=environment,text=True,capture_output=True,timeout=300)
        if child.returncode:
            raise RuntimeError('cold interacting-defect worker failed:\n'+child.stdout+child.stderr)
        value = strict_json.load(path)
    if value.get('source_bindings') != pins or pins != source_bindings():
        raise ValueError('source changed across calculation')
    return value


def validate(value,replay=None):
    if canonical(value) != canonical(report() if replay is None else replay):
        raise ValueError('interacting-defect receipt differs from complete replay')
    if value.get('source_bindings') != source_bindings():
        raise ValueError('interacting-defect source manifest is stale or incomplete')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    choices = parser.add_mutually_exclusive_group(required=True)
    choices.add_argument('--write',nargs='?',const=str(OUTPUT))
    choices.add_argument('--check',nargs='?',const=str(OUTPUT))
    choices.add_argument('--worker')
    args = parser.parse_args()
    if args.worker:
        Path(args.worker).write_text(canonical(produce()))
    elif args.write:
        Path(args.write).write_text(canonical(report()))
        print('wrote '+args.write)
    else:
        validate(strict_json.load(Path(args.check)))
        print('complete interacting-defect replay passed')


if __name__ == '__main__':
    main()
