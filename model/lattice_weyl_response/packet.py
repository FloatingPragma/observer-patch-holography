"""Replay the free-source lattice stress-response attachment."""
import argparse
import hashlib
import json
from pathlib import Path
from tools import strict_json
from . import algebra,bounds,controls

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).parent
PRIOR='model/weyl_gravity_response/results.json'
OUTPUT=HERE/'results.json'


def canonical(value):
    return json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n'


def source_bindings():
    paths={p.relative_to(ROOT).as_posix() for p in HERE.rglob('*.py')}
    paths.update({'model/lattice_weyl_response/README.md',
        'paper/tex_fragments/LATTICE_WEYL_RESPONSE.tex',PRIOR,
        'model/weyl_gravity_response/README.md','model/weyl_gravity_response/algebra.py',
        'model/weyl_gravity_response/bounds.py','model/weyl_geometric_modular/README.md',
        'model/weyl_geometric_modular/algebra.py','model/continuum_field_net/README.md',
        'model/continuum_field_net/cells.py','model/continuum_field_net/valleys.py',
        'model/continuum_free_limit/README.md',
        'Lean/ObserverPatchHolography/EinsteinBranch/SmallBall.lean',
        'Lean/ObserverPatchHolography/EinsteinBranch/RegisterSurface.lean',
        'tools/strict_json.py'})
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def prior_custody():
    value=strict_json.load(ROOT/PRIOR)
    if value.get('new_checks_pass') is not True:
        raise ArithmeticError('prior response scientific checks are absent')
    stale=[p for p,digest in value['source_bindings'].items()
           if not (ROOT/p).is_file() or hashlib.sha256((ROOT/p).read_bytes()).hexdigest()!=digest]
    return {'prior':PRIOR,'direct_prior_source_mismatches':stale,
            'prior_new_checks_pass':True,'prior_full_source_replay':False,
            'inherited_dependency_custody_pass':value['dependency_custody_pass'],
            'dependency_custody_pass':not stale and value['dependency_custody_pass'] is True}


def theorem():
    return {'source_interaction_U':0,'valleys_retained':8,'response_dimension':9,
            'actual_compact_lattice_states_and_quadratic_detectors':True,
            'gapless_projection_and_first_derivative_bounds':True,
            'explicit_joint_lattice_radius_limit':True,
            'finite_spacing_response_invertibility_proved':True,
            'arbitrary_all_valley_state_total_stress_claimed':False,
            'nonzero_U_stress_convergence_claimed':False,
            'lattice_regional_modular_log_convergence_claimed':False,
            'boosted_modular_caps_realized':False,
            'native_gate_program_for_shell_detector_supplied':False,
            'area_Newton_or_stationarity_supplied':False,
            'complete_Einstein_register_inhabited':False}


def validate(value):
    if value.get('schema')!='oph.lattice_weyl_response.v1' or value.get('new_checks_pass') is not True:
        raise ArithmeticError('invalid lattice response certificate')
    expected={'finite_attachment':'FA94','theorem':theorem(),'algebra':algebra.certificate(),
              'bounds':bounds.certificate(),'controls':controls.certificate(),'prior_custody':prior_custody()}
    for key,item in expected.items():
        if canonical(value.get(key))!=canonical(item):
            raise ArithmeticError('lattice response proof, scope or custody changed: '+key)
    passed=expected['prior_custody']['dependency_custody_pass']
    if value.get('dependency_custody_pass') is not passed or value.get('all_pass') is not passed:
        raise ArithmeticError('inherited custody must remain explicit')
    if value.get('source_bindings')!=source_bindings():
        raise ArithmeticError('lattice response current source bindings differ')


def build():
    pins=source_bindings();prior=prior_custody()
    value={'schema':'oph.lattice_weyl_response.v1','finite_attachment':'FA94',
           'new_checks_pass':True,'dependency_custody_pass':prior['dependency_custody_pass'],
           'all_pass':prior['dependency_custody_pass'],'theorem':theorem(),
           'algebra':algebra.certificate(),'bounds':bounds.certificate(),
           'controls':controls.certificate(),'prior_custody':prior,'source_bindings':pins}
    if pins!=source_bindings():raise ArithmeticError('sources changed during replay')
    validate(value)
    return value


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',type=Path,nargs='?',const=OUTPUT)
    group.add_argument('--check',type=Path,nargs='?',const=OUTPUT)
    args=parser.parse_args();value=build()
    if args.write:args.write.write_text(canonical(value))
    else:
        retained=strict_json.load(args.check);validate(retained)
        if canonical(retained)!=canonical(value):raise ArithmeticError('cold replay differs')
    print(json.dumps({k:value[k] for k in ('new_checks_pass','dependency_custody_pass','all_pass')}))


if __name__=='__main__':main()
