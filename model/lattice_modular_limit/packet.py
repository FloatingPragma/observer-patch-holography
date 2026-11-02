"""Replay the fixed-region free-source modular limit."""
import argparse
import hashlib
import json
from pathlib import Path
from tools import strict_json
from . import algebra,bounds,controls

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).parent
OUTPUT=HERE/'results.json'
CHECKPOINT='ba312dbf0'
PRIOR='model/lattice_weyl_response/results.json'
HISTORICAL=HERE/'historical/FA94-ba312dbf0.json'
HISTORICAL_SHA256='ba89e41c2bc0068684f8fbe048941e76d046db33709cfc07b62f4bf11fc2b1df'


def canonical(value):return json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n'


def source_bindings():
    paths={p.relative_to(ROOT).as_posix() for p in HERE.rglob('*.py')}
    paths.update({'model/lattice_modular_limit/README.md',
        HISTORICAL.relative_to(ROOT).as_posix(),
        'paper/tex_fragments/LATTICE_MODULAR_LIMIT.tex',
        'model/lattice_weyl_response/README.md','model/weyl_geometric_modular/README.md',
        'model/continuum_field_net/README.md','model/continuum_field_net/cells.py',
        'model/continuum_free_limit/README.md',
        'Lean/ObserverPatchHolography/EinsteinBranch/SmallBall.lean',
        'Lean/ObserverPatchHolography/EinsteinBranch/RegisterSurface.lean','tools/strict_json.py'})
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def historical_context():
    data=HISTORICAL.read_bytes()
    if hashlib.sha256(data).hexdigest()!=HISTORICAL_SHA256:
        raise ArithmeticError('immutable historical source receipt changed')
    value=strict_json.loads(data.decode())
    if value.get('new_checks_pass') is not True:
        raise ArithmeticError('historical source theorem was not certified')
    return {'checkpoint':CHECKPOINT,'prior':PRIOR,'sha256':hashlib.sha256(data).hexdigest(),
            'historical_new_checks_pass':True,
            'historical_dependency_custody_pass':value['dependency_custody_pass'],
            'historical_full_source_replay_performed':False,
            'current_broad_CI_receipts_reinterpreted':False}


def theorem():
    return {'source_U':0,'region':'fixed bounded ball','valleys':8,
            'finite_regional_covariances_strictly_between_zero_and_one':True,
            'continuum_covariance_has_no_endpoint_eigenvectors':True,
            'actual_modular_flow_converges_for_all_fixed_probes':True,
            'convergence_uniform_on_compact_modular_time_intervals':True,
            'generator_convergence_on_matched_graph_core':True,
            'arbitrary_delta_a_tending_to_zero_allowed':True,
            'finite_CAR_words_and_correlations_converge':True,
            'explicit_unitary_relative_entropies_converge':True,
            'finite_valley_maps_are_exact_isometries':False,
            'uniform_covariance_gap_assumed':False,
            'bare_logs_converge_on_all_unmodified_probes':False,
            'operator_norm_modular_flow_convergence_claimed':False,
            'absolute_continuum_entropy_computed':False,
            'shrinking_region_area_or_geometry_limit_proved':False,
            'nonzero_U_modular_limit_proved':False,
            'native_modular_gate_execution_supplied':False,
            'Einstein_or_particle_identification_endpoint_closed':False}


def validate(value):
    if value.get('schema')!='oph.lattice_modular_limit.v1' or value.get('new_checks_pass') is not True:
        raise ArithmeticError('invalid modular-limit certificate')
    for key,item in {'finite_attachment':'FA96','theorem':theorem(),'algebra':algebra.certificate(),
                     'bounds':bounds.certificate(),'controls':controls.certificate(),
                     'historical_context':historical_context()}.items():
        if canonical(value.get(key))!=canonical(item):
            raise ArithmeticError('modular proof, scope or historical record changed: '+key)
    if value.get('dependency_custody_pass') is not False or value.get('all_pass') is not False:
        raise ArithmeticError('unreplayed historical custody cannot be promoted')
    if value.get('source_bindings')!=source_bindings():
        raise ArithmeticError('current direct mathematical source bindings changed')


def build():
    pins=source_bindings()
    value={'schema':'oph.lattice_modular_limit.v1','finite_attachment':'FA96','new_checks_pass':True,
           'dependency_custody_pass':False,'all_pass':False,'theorem':theorem(),
           'algebra':algebra.certificate(),'bounds':bounds.certificate(),'controls':controls.certificate(),
           'historical_context':historical_context(),'source_bindings':pins}
    if pins!=source_bindings():raise ArithmeticError('source changed during replay')
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
