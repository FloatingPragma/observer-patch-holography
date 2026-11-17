"""Cold, source-bound positive energy readout and spectral-tail certificates."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from fractions import Fraction as F

from tools import strict_json
from tools.build_completion_program import integration_source_paths
from . import energy,finite,kernels


ROOT=Path(__file__).resolve().parents[2]
OUTPUT=Path(__file__).resolve().parent/'results.json'


def source_bindings():
    paths=integration_source_paths(ROOT,{
        'domain_and_proof':'model/spectral_readout/README.md',
        'source_receipts':['model/retarded_matter_response/results.json',
                           'model/ground_spectral_source/results.json',
                           'model/retarded_spatial_locality/results.json'],
        'independent_checks':[]})
    paths.add('paper/tex_fragments/SPECTRAL_READOUT.tex')
    return {path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(paths)}


def validate(value):
    if value.get('schema')!='oph.spectral_readout.v1' or value.get('all_pass') is not True:
        raise ValueError('invalid positive spectral readout packet')
    if value['spectral_measure']['second_operator_moment']!='h^2+2U r h+U^2 r I':
        raise ValueError('contact second moment changed')
    if value['finite_volume_error_contract']['total_infinite_ground_bin_error'] is not None:
        raise ValueError('unknown state discrepancy silently discharged')
    if value['finite_volume_error_contract']['state_discrepancy'] is not None:
        raise ValueError('retained packet has no measured state discrepancy')
    for row,m in zip(value['finite_source_readouts'],(8,16,64),strict=True):
        if (row['outcome_count']!=m or len(row['rows'])!=m or
                row['source_lines_used']!=78 or row['all_outcomes_retained'] is not True or
                row['direct_positive_and_Fourier_enclosures_agree'] is not True):
            raise ValueError('complete positive spectral output changed')
        total_low=total_high=F(0)
        for index,item in enumerate(row['rows']):
            low,high=map(F,item['probability_interval'])
            if item['index']!=index or not 0<=low<=high<=1:
                raise ValueError('invalid complete probability interval')
            total_low+=low;total_high+=high
        if not total_low<=1<=total_high:raise ValueError('complete readout lost probability')
        for key in ('native_full_energy_POVM_executed','finite_shot_records_generated',
                    'infinite_ground_values_inferred','energy_centers_are_relativistic_masses'):
            if row[key] is not False:raise ValueError('finite spectral calculation scope promoted')
    measure=value['spectral_measure']
    for key in ('common_Hilbert_space_strong_resolvent_convergence_of_dilations_claimed',
                'infinite_ground_spectral_values_computed','real_axis_density_or_isolated_pole_inferred',
                'relativistic_rest_mass_inferred'):
        if measure[key] is not False:raise ValueError('spectral convergence scope promoted')
    if value['native_first_harmonic']['all_higher_harmonics_natively_executed'] is not False:
        raise ValueError('one native harmonic promoted to a complete energy measurement')


def canonical(value):return json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n'


def _produce():
    pins=source_bindings()
    examples=[]
    for m in (1024,4096,16384):
        values=kernels.real_energy_leakage(m,F(1,4),F(1,2),F(67,128))
        examples.append({'outcome_count':m,'time_step':'1/4','resolution':'1/2',
                         **{name:str(value) for name,value in values.items()}})
    value={'schema':'oph.spectral_readout.v1','all_pass':True,
           'spectral_measure':energy.certificate(),
           'finite_source_readouts':[finite.readout(m) for m in (8,16,64)],
           'infinite_ground_uniform_resolution_examples':examples,
           'native_first_harmonic':finite.native_first_harmonic(),
           'finite_volume_error_contract':energy.comparison_contract(
               8,F(1,8),1,half_side=4096,buffer_radius=64,support_radius=0,support_sites=1),
           'source_bindings':pins,
           'proof_scope':'Complete positive finite-resolution energy POVM; exact finite-source probabilities and energy tails; same-subsequence weak spectral and strong compressed-resolvent convergence. No numerical infinite ground, pole, native higher-harmonic execution or relativistic mass claim.'}
    if source_bindings()!=pins:raise ValueError('source changed during spectral replay')
    validate(value)
    return value


def report():
    pins=source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-spectral-readout-') as temp:
        target=Path(temp)/'worker.json';environment=dict(os.environ)
        environment['OPENBLAS_NUM_THREADS']=environment['OMP_NUM_THREADS']='1'
        worker=subprocess.run([sys.executable,'-W','error','-X','pycache_prefix='+str(Path(temp)/'pycache'),
                               '-m','model.spectral_readout.packet','--worker',str(target)],
                              cwd=ROOT,env=environment,capture_output=True,text=True)
        if worker.returncode:raise RuntimeError('cold spectral readout failed: '+worker.stderr)
        value=strict_json.load(target)
    if source_bindings()!=pins or value['source_bindings']!=pins:
        raise ValueError('source changed across spectral replay')
    validate(value)
    return value


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    modes=parser.add_mutually_exclusive_group()
    modes.add_argument('--write',type=Path,nargs='?',const=OUTPUT)
    modes.add_argument('--check',type=Path,nargs='?',const=OUTPUT)
    modes.add_argument('--worker',type=Path,help=argparse.SUPPRESS)
    args=parser.parse_args();value=_produce() if args.worker else report()
    if args.check:
        previous=strict_json.load(args.check);validate(previous)
        if canonical(previous)!=canonical(value):raise ValueError('spectral readout receipt changed')
    elif args.worker or args.write:
        (args.worker or args.write).write_text(canonical(value))
    if not args.worker:print('Verified complete positive spectral readouts: M=8,16,64; no physical mass assignment.')


if __name__=='__main__':main()
