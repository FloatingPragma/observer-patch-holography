"""Reproduce the exact small-ball entropy theorem and current narrow bindings."""
import argparse
import hashlib
import json
from pathlib import Path
from tools import strict_json
from . import algebra,bounds,controls

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).parent
OUTPUT=HERE/'results.json'


def canonical(value):
    return json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n'


def source_bindings():
    paths={p.relative_to(ROOT).as_posix() for p in HERE.rglob('*.py')}
    paths.update({'model/__init__.py','tools/strict_json.py',
        'model/small_ball_entropy_asymptotic/README.md',
        'paper/tex_fragments/SMALL_BALL_ENTROPY_ASYMPTOTIC.tex',
        'model/shell_entropy_response/__init__.py','model/shell_entropy_response/bounds.py',
        'model/shell_entropy_response/README.md',
        'model/shrinking_modular_entropy/__init__.py','model/shrinking_modular_entropy/bounds.py',
        'model/shrinking_modular_entropy/README.md','model/weyl_gravity_response/README.md',
        'model/weyl_geometric_modular/README.md','model/lattice_modular_limit/README.md',
        'Lean/ObserverPatchHolography/EinsteinBranch/SmallBall.lean',
        'Lean/ObserverPatchHolography/EinsteinBranch/RegisterSurface.lean'})
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def theorem():
    return {'source_U':0,'source_valleys_retained':8,
        'same_normalized_exterior_shell_preparations':True,
        'actual_cell_norm_used':True,
        'operator_exponential_domain_proved':True,
        'fractional_covariance_inverse_moments_proved':True,
        'fixed_amplitude_nonlinear_remainder_o_R4_at_fixed_L':True,
        'all81_rest_ball_preparation_responses_controlled':True,
        'exact_leading_relative_entropy_coefficient':True,
        'actual_Araki_relative_entropy':True,
        'finite_lattice_density_entropy_convergence':True,
        'uniform_covariance_gap_assumed':False,
        'arbitrary_L_R_schedule_o_R4_asserted':False,
        'sharp_boosted_lattice_algebra_constructed':False,
        'area_Newton_or_stationarity_derived':False,
        'interacting_Einstein_endpoint_closed':False,
        'nonzero_U_entropy_theorem':False,
        'native_entropy_readout_executed':False,
        'measured_particle_species_identified':False}


def custody():
    return {'current_narrow_sources_checked':True,'old_receipts_rewritten':False,
            'runtime_Git_history_required':False,'full_repository_replay_performed':False}


def validate(value):
    keys={'schema','finite_attachment','new_checks_pass','theorem','bounds','algebra','controls','custody','source_bindings'}
    if type(value) is not dict or set(value)!=keys:raise ArithmeticError('unexpected asymptotic receipt root fields')
    if value['schema']!='oph.small_ball_entropy_asymptotic.v1' or value['finite_attachment']!='FA107':
        raise ArithmeticError('wrong small-ball theorem identity')
    if value['new_checks_pass'] is not True:raise ArithmeticError('strict check flag required')
    expected={'theorem':theorem(),'bounds':bounds.certificate(),'algebra':algebra.certificate(),
              'controls':controls.certificate(),'custody':custody(),'source_bindings':source_bindings()}
    for key,item in expected.items():
        if canonical(value[key])!=canonical(item):raise ArithmeticError('asymptotic certificate differs: '+key)


def build():
    pins=source_bindings()
    value={'schema':'oph.small_ball_entropy_asymptotic.v1','finite_attachment':'FA107',
           'new_checks_pass':True,'theorem':theorem(),'bounds':bounds.certificate(),
           'algebra':algebra.certificate(),'controls':controls.certificate(),'custody':custody(),
           'source_bindings':pins}
    if source_bindings()!=pins:raise ArithmeticError('source changed during theorem replay')
    validate(value)
    return value


def main():
    p=argparse.ArgumentParser(description=__doc__);g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',type=Path,nargs='?',const=OUTPUT)
    g.add_argument('--check',type=Path,nargs='?',const=OUTPUT)
    a=p.parse_args();value=build()
    if a.write:a.write.write_text(canonical(value))
    else:
        old=strict_json.load(a.check);validate(old)
        if canonical(old)!=canonical(value):raise ArithmeticError('cold theorem replay differs')
    print(json.dumps({'new_checks_pass':True,'full_repository_replay_performed':False}))


if __name__=='__main__':main()
