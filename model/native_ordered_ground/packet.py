"""Cold full-native replay and strict source custody for FA98."""
import argparse
import ast
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import numpy as np
from tools import strict_json
from . import experiment,source,angles

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).parent
OUTPUT=HERE/'results.json'
REPLAY_ATOL=2e-11


def canonical(x):return json.dumps(x,indent=2,sort_keys=True,allow_nan=False)+'\n'


def compare(expected,actual,path='receipt'):
    if type(expected) is not type(actual):raise ValueError(path+': type changed')
    if type(actual) is dict:
        if expected.keys()!=actual.keys():raise ValueError(path+': keys changed')
        for k in actual:compare(expected[k],actual[k],path+'.'+k)
    elif type(actual) is list:
        if len(expected)!=len(actual):raise ValueError(path+': length changed')
        for i,(a,b) in enumerate(zip(expected,actual,strict=True)):compare(a,b,path+str(i))
    elif type(actual) is float:
        if not math.isfinite(actual) or not math.isfinite(expected):raise ValueError('nonfinite diagnostic')
        # All floating fields are norms, masses or probabilities; exact signed
        # quantities are rational strings and never use this tolerance.
        if min(actual,expected)<0:raise ValueError(path+': negative diagnostic')
        if path.rsplit('.',1)[-1] in ('complete_first_failure_probability','rounded_negative_mass'):
            if abs(math.sqrt(actual)-math.sqrt(expected))>REPLAY_ATOL:
                raise ValueError(path+': failure amplitude changed')
        elif abs(actual-expected)>REPLAY_ATOL:raise ValueError(path+': diagnostic changed')
    elif type(actual) not in (str,int,bool,type(None)) or expected!=actual:raise ValueError(path+': value changed')


def compare_arrays(retained,current,key):
    """Replay agreement only; validate independently checks the source target.

    The failure budget is global in square-root weights, not per outcome. Thus
    many individually tiny forged outcomes cannot evade the total error bound.
    """
    if retained.shape!=current.shape or retained.dtype!=current.dtype:
        raise ValueError('native replay array carrier changed')
    if not np.isfinite(retained).all() or not np.isfinite(current).all():
        raise ValueError('nonfinite native replay array')
    if key=='failure_weights':
        if np.min(retained,initial=0)<0 or np.min(current,initial=0)<0:
            raise ValueError('negative native replay failure weight')
        if np.linalg.norm(np.sqrt(retained)-np.sqrt(current))>REPLAY_ATOL:
            raise ValueError('cold total failure amplitude changed')
    elif not np.allclose(retained,current,atol=REPLAY_ATOL,rtol=0):
        raise ValueError('cold native array changed')


def source_bindings():
    paths=set(HERE.glob('*.py'));queue=list(paths)
    legacy=ROOT/'evidence/observer_dynamics_20260925/oph-physics-sim'
    def include(path):
        if path.is_file() and path not in paths:paths.add(path);queue.append(path)
    for path in (ROOT/'model/__init__.py',HERE/'__init__.py'):include(path)
    # slow_geometry.source loads this declared archive module dynamically.
    include(legacy/'oph_exact/carrier.py');include(legacy/'oph_exact/__init__.py')
    while queue:
        path=queue.pop();base=legacy if path.is_relative_to(legacy) else ROOT
        parts=path.relative_to(base).with_suffix('').parts
        package='.'.join(parts[:-1]);tree=ast.parse(path.read_text())
        for node in ast.walk(tree):
            names=[]
            if isinstance(node,ast.Import):names=[a.name for a in node.names]
            if isinstance(node,ast.ImportFrom):
                name=('.'*node.level)+(node.module or '')
                base=importlib.util.resolve_name(name,package) if node.level else name
                names=[base]+[base+'.'+a.name for a in node.names]
            for name in names:
                root=legacy if name.split('.')[0] in ('oph_fpe','oph_exact') else ROOT
                segments=name.split('.')
                include(root.joinpath(*segments).with_suffix('.py'))
                for j in range(1,len(segments)+1):include(root.joinpath(*segments[:j])/'__init__.py')
    paths.update(ROOT/p for p in ('code/source_selection_model/response.json',
       'evidence/native_geometric_source_20260925/spec.json',
       'code/m1_source_realization/receipt.json','model/native_spectral_detector/contract.cpp',
       'model/native_ordered_ground/README.md','model/native_ordered_ground/witness.npz',
       'paper/tex_fragments/NATIVE_ORDERED_GROUND.tex','model/spin_goldstone_sector/README.md',
       'model/protected_spin_limit/README.md','model/native_spin_observer/README.md'))
    return {p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}


def theorem():
    return {'finite_attachment':'FA98','finite_bipartite_positive_field_ground_preparation':True,
      'exact_rational_ground_algorithm_and_gap_lower_bound':True,
      'arbitrary_precision_rational_angle_construction':True,
      'complete_native_finite_ground_and_pointer_executed':True,
      'all_256_finite_eigenvalues_checked':True,'supplied_ground_amplitude_oracle':False,
      'ordered_thermodynamic_then_zero_field_limits':True,
      'effective_infinite_ground_marginal_modulus':False,
      'finite_field_state_stationary_for_zero_field_source':False,
      'continuing_autonomous_reference_realized':False,'relativistic_rest_mass_identified':False,
      'physical_species_or_Einstein_dynamics_derived':False,
      'arbitrary_precision_interface':'FA93 exact dyadic integrated port commands',
      'local_source_ground_defect':'B*s+24*J*s*epsilon+6*J/E',
      'historical_context_commit':'ba312dbf0674abfb8ed08366953fe5e6d3e13202',
      'historical_broad_receipts_replayed_here':False}


def exact_controls():
    bound=source.perron_bound(8,source.cube()[0],__import__('fractions').Fraction(1),__import__('fractions').Fraction(1,4))
    nodes=[]
    for a,b,bits in ((1,1,40),(3,7,80),(1,10**30,120),(0,1,40),(1,0,40)):
        row=angles.node(a,b,bits)
        nodes.append({'a':a,'b':b,'bits':bits,**{k:str(v) for k,v in row.items()}})
    return {'Perron_bound':bound,'rational_node_controls':nodes}


def produce():
    pins=source_bindings();checked,execution,arrays=experiment.run()
    value={'schema':'oph.native_ordered_ground.v1','new_checks_pass':True,'all_pass':True,
       'theorem':theorem(),'finite_spectrum':checked,'exact_controls':exact_controls(),
       'execution':execution,'source_bindings':pins}
    if pins!=source_bindings():raise ValueError('source changed during native replay')
    validate(value,arrays);return value,arrays


def validate(value,arrays):
    keys={'schema','new_checks_pass','all_pass','theorem','finite_spectrum','exact_controls','execution','source_bindings'}
    if set(value)!=keys or value['schema']!='oph.native_ordered_ground.v1' or value['new_checks_pass'] is not True or value['all_pass'] is not True:
        raise ValueError('invalid complete native ground receipt')
    word,initial,target,checked,plan=experiment.plan()
    for key,expected in [('theorem',theorem()),('finite_spectrum',checked),('exact_controls',exact_controls()),('source_bindings',source_bindings())]:
        if canonical(value[key])!=canonical(expected):raise ValueError('exact proof/source/type mismatch: '+key)
    if set(arrays)!={'success','failure_weights'}:raise ValueError('complete native arrays required')
    success=arrays['success'];weights=arrays['failure_weights']
    if success.shape!=initial.shape or success.dtype!=np.complex128 or weights.shape!=(1090,4) or weights.dtype!=np.float64:
        raise ValueError('complete native array carrier changed')
    if not np.isfinite(success).all() or not np.isfinite(weights).all() or np.min(weights,initial=0)<0:
        raise ValueError('invalid native probabilities')
    pair=np.array([e['kind']=='J' for e in word])
    if np.any(weights[pair,2:]!=0):raise ValueError('nonphysical padded outcome')
    residual=abs(float(np.vdot(success,success).real)+float(weights.sum())-1)
    distance=float(np.sqrt(np.linalg.norm(success-target)**2+weights.sum()))
    if residual>2e-8 or distance>float(__import__('fractions').Fraction(plan['complete_ground_isometry_error']))+2e-8:
        raise ValueError('complete output misses certified source target')
    row=value['execution'];count=row.get('rounded_negative_count');mass=row.get('rounded_negative_mass')
    if type(count) is not int or not 0<=count<=4360 or type(mass) is not float or not 0<=mass<1e-35:
        raise ValueError('invalid roundoff diagnostic')
    expected=dict(plan,complete_output_distance_diagnostic=distance,
      complete_first_failure_probability=float(weights.sum()),normalization_residual=residual,
      native_first_failure_slots=4360-2*int(pair.sum()),rounded_negative_count=count,rounded_negative_mass=mass,
      pointer_bit_zero_success_probability=float(np.linalg.norm(success[:1<<16])**2),
      pointer_bit_one_success_probability=float(np.linalg.norm(success[1<<16:])**2))
    compare(expected,row)


def write(path,value,arrays):
    archive=path.with_suffix('.npz');np.savez_compressed(archive,**arrays)
    entry={'filename':archive.name,'sha256':hashlib.sha256(archive.read_bytes()).hexdigest()}
    path.write_text(canonical(dict(value,complete_record_archive=entry)))


def read(path):
    value=strict_json.load(path);entry=value.pop('complete_record_archive')
    if set(entry)!={'filename','sha256'} or entry['filename']!=path.with_suffix('.npz').name:raise ValueError('local paired archive required')
    archive=path.with_suffix('.npz')
    if hashlib.sha256(archive.read_bytes()).hexdigest()!=entry['sha256']:raise ValueError('native archive digest mismatch')
    with np.load(archive,allow_pickle=False) as data:arrays={k:data[k] for k in data.files}
    validate(value,arrays);return value,arrays


def report():
    pins=source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-native-ground-') as tmp:
        output=Path(tmp)/'worker.json';env=dict(os.environ)
        for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):env[key]='1'
        p=subprocess.run([sys.executable,'-W','error','-X','pycache_prefix='+str(Path(tmp)/'cache'),'-m',
           'model.native_ordered_ground','--worker',str(output)],cwd=ROOT,env=env,capture_output=True,text=True)
        if p.returncode:raise RuntimeError(p.stdout+p.stderr)
        value,arrays=read(output)
    if pins!=source_bindings():raise ValueError('source changed during cold replay')
    return value,arrays


def main():
    parser=argparse.ArgumentParser(description=__doc__);group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',type=Path,nargs='?',const=OUTPUT);group.add_argument('--check',type=Path,nargs='?',const=OUTPUT)
    group.add_argument('--worker',type=Path);args=parser.parse_args()
    value,arrays=produce() if args.worker else report()
    if args.check:
        old,retained=read(args.check);compare(old,value)
        for k in arrays:
            compare_arrays(retained[k],arrays[k],k)
    else:write(args.worker or args.write,value,arrays)
    if not args.worker:print(canonical({'new_checks_pass':True,'all_pass':True}),end='')
