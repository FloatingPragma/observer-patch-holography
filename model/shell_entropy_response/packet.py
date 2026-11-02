"""Reproduce the positive shell entropy-response certificate."""
import argparse
import hashlib
import json
from pathlib import Path
from tools import strict_json
from . import algebra, bounds, controls

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).parent
OUTPUT=HERE/'results.json'


def canonical(value):
    return json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n'


def source_bindings():
    paths={p.relative_to(ROOT).as_posix() for p in HERE.rglob('*.py')}
    paths.update({
        'model/__init__.py','tools/strict_json.py',
        'model/shell_entropy_response/README.md',
        'paper/tex_fragments/SHELL_ENTROPY_RESPONSE.tex',
        'model/shrinking_modular_entropy/__init__.py',
        'model/shrinking_modular_entropy/algebra.py',
        'model/shrinking_modular_entropy/README.md',
        'model/shrinking_modular_entropy/bounds.py',
        'model/weyl_gravity_response/README.md',
        'model/weyl_gravity_response/bounds.py',
        'model/lattice_modular_limit/README.md',
        'Lean/ObserverPatchHolography/EinsteinBranch/SmallBall.lean',
        'Lean/ObserverPatchHolography/EinsteinBranch/RegisterSurface.lean',
    })
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def theorem():
    return {'source_U':0,'source_valleys_retained':8,
            'same_natural_exterior_shell_states':True,'actual_cell_norm_used':True,
            'rank_one_noncommuting_entropy_Hessian':True,
            'strict_positive_shell_lower_bound':True,
            'finite_amplitude_cubic_remainder':True,
            'actual_Araki_relative_entropy_limit':True,
            'Hiai_Theorem_4_1_v_imported':True,
            'uniform_compact_amplitude_lattice_Hessian_limit':True,
            'positive_symbolic_joint_mesh':True,
            'uniform_covariance_gap_assumed':False,
            'absolute_continuum_entropy_computed':False,
            'global_affine_pure_vacuum_family_asserted':False,
            'nonzero_U_entropy_limit_proved':False,
            'area_dynamics_or_stationarity_derived':False,
            'Einstein_endpoint_closed':False,
            'native_entropy_circuit_executed':False,
            'astronomical_refinement_numerically_executed':False}


def custody():
    return {'current_narrow_mathematical_sources_checked':True,
            'predecessor_full_repository_receipts_replayed':False,
            'runtime_Git_history_required':False,
            'scope':'analytic theorem, exact finite controls, and separately identified floating diagnostics'}


def validate(value):
    if type(value) is not dict or set(value)!={'schema','finite_attachment','new_checks_pass',
            'full_repository_replay_performed','theorem','bounds','algebra','controls','custody','source_bindings'}:
        raise ArithmeticError('unexpected entropy-response root fields')
    if value.get('schema')!='oph.shell_entropy_response.v1' or value.get('finite_attachment')!='FA106':
        raise ArithmeticError('wrong entropy-response identity')
    if value.get('new_checks_pass') is not True or value.get('full_repository_replay_performed') is not False:
        raise ArithmeticError('unperformed check cannot be promoted')
    expected={'theorem':theorem(),'bounds':bounds.certificate(),'algebra':algebra.certificate(),
              'controls':controls.certificate(),'custody':custody(),'source_bindings':source_bindings()}
    for key,item in expected.items():
        if canonical(value.get(key))!=canonical(item):
            raise ArithmeticError('entropy-response certificate differs: '+key)


def build():
    pins=source_bindings()
    result={'schema':'oph.shell_entropy_response.v1','finite_attachment':'FA106','new_checks_pass':True,
            'full_repository_replay_performed':False,'theorem':theorem(),'bounds':bounds.certificate(),
            'algebra':algebra.certificate(),'controls':controls.certificate(),'custody':custody(),
            'source_bindings':pins}
    if source_bindings()!=pins:raise ArithmeticError('source changed during replay')
    validate(result)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',type=Path,nargs='?',const=OUTPUT)
    group.add_argument('--check',type=Path,nargs='?',const=OUTPUT)
    args=parser.parse_args();result=build()
    if args.write:args.write.write_text(canonical(result))
    else:
        retained=strict_json.load(args.check);validate(retained)
        if canonical(retained)!=canonical(result):raise ArithmeticError('cold replay differs')
    print(json.dumps({'new_checks_pass':True,'full_repository_replay_performed':False}))


if __name__=='__main__':main()
