"""Cold source-bound replay for the full centered bounded-contact free limit."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from tools import strict_json
from tools.build_completion_program import integration_source_paths
from . import algebra,bounds

ROOT=Path(__file__).resolve().parents[2]
OUTPUT=Path(__file__).resolve().parent/'results.json'
PRIORS=('centered_source_ground','continuum_field_net',
        'continuum_composite_channel','continuum_ground_constraints')


def source_bindings():
    paths=integration_source_paths(ROOT,{
        'domain_and_proof':'model/continuum_free_limit/README.md',
        'source_receipts':[f'model/{name}/results.json' for name in PRIORS],
        'independent_checks':[],
    })
    paths.add('paper/tex_fragments/CONTINUUM_FREE_LIMIT.tex')
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def prior_proofs():
    for name in PRIORS:
        row=strict_json.loads((ROOT/f'model/{name}/results.json').read_text())
        if row.get('all_pass') is not True or not row.get('source_bindings'):
            raise ValueError('prior theorem missing source custody')
        for path,digest in row['source_bindings'].items():
            if hashlib.sha256((ROOT/path).read_bytes()).hexdigest()!=digest:
                raise ValueError('prior proof source changed: '+path)
    return {'all_prior_source_hashes_verified':True,
            'prior_expensive_experiments_reexecuted':False}


def theorem():
    return {'reference':'new particle-hole invariant centered ground of K=H-(U/2)N',
            'refinement':'a->0, c(a)=v/a, fixed v>0, 0<=U(a)<=Umax<infinity',
            'uniform_composite_finite_energy_weight_vanishes':True,
            'both_vacuum_charge_branches_controlled':True,
            'initial_state_is_pure_free_filled_sea':True,
            'outside_initial_CAR_vacuum_cyclic':True,
            'graded_locality_promotes_vacuum_to_operator_identity':True,
            'complete_centered_field_hierarchy_is_free':True,
            'massless_Weyl_field_count':8,'Weyl_chiralities_each_sign':4,
            'particle_and_hole_energy':'v |p|','rest_masses_in_this_branch':'all zero',
            'nonzero_continuum_interaction_survives':False,
            'massive_sector_supplied':False,
            'equal_time_composite_CAR_field_set_to_zero':False,
            'original_uncentered_ground_reidentified':False,
            'unbounded_coupling_trajectory_decided':False,
            'external_renormalization_theorem_imported':False,
            'physical_species_or_unit_dictionary_identified':False,
            'Einstein_requirements_completed':False}


def validate(value):
    if value.get('schema')!='oph.continuum_free_limit.v1' or value.get('all_pass') is not True:
        raise ValueError('invalid centered free-limit receipt')
    row=value['theorem']
    for key in ('uniform_composite_finite_energy_weight_vanishes',
                'both_vacuum_charge_branches_controlled','initial_state_is_pure_free_filled_sea',
                'outside_initial_CAR_vacuum_cyclic',
                'graded_locality_promotes_vacuum_to_operator_identity',
                'complete_centered_field_hierarchy_is_free'):
        if row[key] is not True:raise ValueError('free hierarchy proof premise lost')
    for key in ('nonzero_continuum_interaction_survives','massive_sector_supplied',
                'equal_time_composite_CAR_field_set_to_zero','original_uncentered_ground_reidentified',
                'unbounded_coupling_trajectory_decided','external_renormalization_theorem_imported',
                'physical_species_or_unit_dictionary_identified','Einstein_requirements_completed'):
        if row[key] is not False:raise ValueError('centered bounded-contact scope was promoted')
    if row['massless_Weyl_field_count']!=8 or row['Weyl_chiralities_each_sign']!=4:
        raise ValueError('a Weyl valley was deleted')
    if row['rest_masses_in_this_branch']!='all zero':raise ValueError('free mode mass changed')
    for key in ('uniform_finite_energy_composite_escape',
                'free_composite_has_no_linear_wick_component','all_eight_weyl_nodes_retained',
                'vacuum_duhamel_requires_uniform_input_mass_bound',
                'vacuum_duhamel_applies_to_both_adjoint_branches'):
        if value['bounds'][key] is not True:raise ValueError('uniform spectral bound premise lost')
    if value['controls']['vacuum_only_countercontrol']['vacuum_motion_alone_proves_full_free_dynamics'] is not False:
        raise ValueError('vacuum-only countercontrol lost')


def canonical(value):return json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n'


def _produce():
    pins=source_bindings()
    value={'schema':'oph.continuum_free_limit.v1','all_pass':True,
           'prior_custody':prior_proofs(),'theorem':theorem(),
           'bounds':bounds.certificate(),'controls':algebra.certificate(),
           'source_bindings':pins,
           'proof_scope':'The self-contained analytic proof determines the infinite hierarchy. Exact bounds and finite full-Fock countercontrols verify its constants and indispensable distinctions; finite tests are not substituted for the proof.'}
    validate(value)
    if source_bindings()!=pins:raise ValueError('source changed during free-limit replay')
    return value


def report():
    pins=source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-continuum-free-') as temp:
        path=Path(temp)/'worker.json';env=dict(os.environ)
        env['OPENBLAS_NUM_THREADS']=env['OMP_NUM_THREADS']='1'
        subprocess.run([sys.executable,'-W','error','-X','pycache_prefix='+str(Path(temp)/'pycache'),
                        '-m','model.continuum_free_limit.packet','--worker',str(path)],
                       cwd=ROOT,env=env,check=True,capture_output=True,text=True)
        value=strict_json.loads(path.read_text())
    if source_bindings()!=pins or value['source_bindings']!=pins:
        raise ValueError('source changed across cold free-limit replay')
    validate(value);return value


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',nargs='?',const=OUTPUT,type=Path)
    group.add_argument('--check',nargs='?',const=OUTPUT,type=Path)
    group.add_argument('--worker',type=Path)
    args=parser.parse_args()
    if args.worker:args.worker.write_text(canonical(_produce()));return
    value=report()
    if args.write:args.write.write_text(canonical(value))
    else:
        retained=strict_json.loads(args.check.read_text());validate(retained)
        if canonical(value)!=canonical(retained):raise ValueError('free-limit replay changed')
    print(json.dumps({'all_pass':True,'source_bindings':len(value['source_bindings']),
                      'centered_bounded_contact_limit':'free eight-valley Weyl'}))


if __name__=='__main__':main()
