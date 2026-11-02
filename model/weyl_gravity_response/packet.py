"""Replay the explicit same-vacuum stress response attachment."""
import argparse
import hashlib
import json
from pathlib import Path
from tools import strict_json
from . import algebra,bounds

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).parent
PRIOR='model/weyl_geometric_modular/results.json'
OUTPUT=HERE/'results.json'


def canonical(value):
    return json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n'


def source_bindings():
    paths={p.relative_to(ROOT).as_posix() for p in HERE.rglob('*.py')}
    paths.update({'model/weyl_gravity_response/README.md',
                  'paper/tex_fragments/WEYL_GRAVITY_RESPONSE.tex',PRIOR,
                  'model/weyl_geometric_modular/README.md',
                  'model/weyl_geometric_modular/algebra.py',
                  'model/weyl_geometric_modular/entropy.py',
                  'model/weyl_geometric_modular/geometry.py',
                  'Lean/ObserverPatchHolography/EinsteinBranch/Tensor.lean',
                  'Lean/ObserverPatchHolography/EinsteinBranch/SmallBall.lean',
                  'Lean/ObserverPatchHolography/EinsteinBranch/RegisterSurface.lean',
                  'tools/strict_json.py'})
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def prior_custody():
    value=strict_json.load(ROOT/PRIOR)
    if value.get('new_checks_pass') is not True:
        raise ArithmeticError('prior modular scientific checks are not present')
    stale=[p for p,digest in value['source_bindings'].items()
           if not (ROOT/p).is_file() or hashlib.sha256((ROOT/p).read_bytes()).hexdigest()!=digest]
    return {'prior':PRIOR,'direct_prior_source_mismatches':stale,
            'prior_new_checks_pass':value['new_checks_pass'],
            'prior_dependency_custody_pass':value['dependency_custody_pass'],
            'prior_full_source_replay':False,
            'historical_broad_mismatch_count':value['custody']['current_prior_source_mismatch_count'],
            'dependency_custody_pass':not stale and value['dependency_custody_pass'] is True}


def theorem():
    return {'variation_dimension':9,'local_normal_two_sided_states_constructed':True,
            'all_tracefree_symmetric_stress_tangents_spanned':True,
            'finite_boosted_ball_modular_responses_invertible':True,
            'finite_relative_entropy_differences_invertible':True,
            'same_stress_and_vacuum_for_all_nine_observers':True,
            'Frechet_matter_first_law_on_common_variation_space':True,
            'explicit_common_continuum_radius_and_perturbation_limits':True,
            'all_eight_Weyl_summands_retained':True,
            'global_pure_vacuum_affine_two_sided_state_family_claimed':False,
            'ten_component_trace_direction_supplied':False,
            'native_shell_operations_implemented':False,
            'joint_lattice_collar_limit_supplied':False,
            'area_Newton_or_generalized_stationarity_supplied':False,
            'geometry_stress_null_balance_supplied':False,
            'complete_Einstein_register_inhabited':False,
            'interacting_or_massive_continuum_supplied':False,
            'measured_particle_species_identified':False}


def validate(value):
    if value.get('schema')!='oph.weyl_gravity_response.v1' or value.get('new_checks_pass') is not True:
        raise ArithmeticError('invalid new response certificate')
    for key,expected in [('theorem',theorem()),('algebra',algebra.certificate()),
                         ('bounds',bounds.certificate()),('prior_custody',prior_custody())]:
        if canonical(value.get(key))!=canonical(expected):
            raise ArithmeticError('response proof, scope or custody changed: '+key)
    passed=prior_custody()['dependency_custody_pass']
    if value.get('dependency_custody_pass') is not passed or value.get('all_pass') is not passed:
        raise ArithmeticError('inherited custody must remain explicit')
    if value.get('source_bindings')!=source_bindings():
        raise ArithmeticError('response source bindings changed')


def build():
    pins=source_bindings();prior=prior_custody()
    value={'schema':'oph.weyl_gravity_response.v1','finite_attachment':'FA91',
           'new_checks_pass':True,'dependency_custody_pass':prior['dependency_custody_pass'],
           'all_pass':prior['dependency_custody_pass'],'theorem':theorem(),
           'algebra':algebra.certificate(),'bounds':bounds.certificate(),
           'prior_custody':prior,'source_bindings':pins}
    if pins!=source_bindings():raise ArithmeticError('source changed during response replay')
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
        if canonical(retained)!=canonical(value):raise ArithmeticError('response receipt differs from cold replay')
    print(json.dumps({k:value[k] for k in ('new_checks_pass','dependency_custody_pass','all_pass')}))


if __name__=='__main__':main()
