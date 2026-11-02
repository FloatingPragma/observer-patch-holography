"""Narrow-input cold receipt for the spatial and normalized scaling theorem."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from tools import strict_json
from . import algebra, bounds, experiment

ROOT=Path(__file__).resolve().parents[2]
OUTPUT=Path(__file__).with_name('results.json')
PRIORS=('model/spin_goldstone_sector/results.json','model/protected_spin_limit/results.json')
CHECKPOINT='ba312dbf0'
SNAPSHOTS={
    PRIORS[0]:('model/spin_particle_scaling/historical/FA92.json','509181dfd80bec6863c2f8fcfa637677a8b3f73ebc3fb5afb869be6c1ba3b575'),
    PRIORS[1]:('model/spin_particle_scaling/historical/FA90.json','0d8a144cab3c774166ddacf421c663e2908e155159395446b76bbcfe2552dbe6'),
}
INPUTS=('model/spin_goldstone_sector/README.md','model/spin_goldstone_sector/experiment.py',
        'model/protected_spin_limit/README.md','model/cubic_particle_channels/README.md',
        'model/field_symmetry/symmetry.py','tools/strict_json.py',
        'paper/tex_fragments/SPIN_PARTICLE_SCALING.tex')


def canonical(value):
    return json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n'


def source_bindings():
    paths=set(INPUTS)
    paths.update(path for path,_ in SNAPSHOTS.values())
    paths.update(str(p.relative_to(ROOT)) for p in Path(__file__).parent.iterdir() if p.suffix in ('.py','.md'))
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def inherited():
    out={}
    for name in PRIORS:
        snapshot,expected=SNAPSHOTS[name]
        data=(ROOT/snapshot).read_bytes()
        if hashlib.sha256(data).hexdigest()!=expected:raise ValueError('historical receipt snapshot changed')
        row=strict_json.loads(data.decode())
        if row.get('new_checks_pass') is not True:raise ValueError('missing inherited mathematical checks')
        out[name]={'checkpoint':CHECKPOINT,'snapshot':snapshot,'sha256':expected,
                   'historical_new_checks_pass':True,'historical_all_pass':row.get('all_pass'),
                   'parent_full_source_closure_replayed':False,'expensive_parent_experiments_reexecuted':False}
    return out


def theorem():
    return {'source':'FA90 c=J E^4,U=J E^8, fixed-spacing cubic AF spin limit',
            'ground':'FA92 selected ordered zero-field state, internal n=e1',
            'two_orthogonal_equal_weight_transverse_channels':True,
            'source_spatial_quantum_numbers_resolved':True,
            'nonzero_normalized_energy_momentum_scaling_measures':True,
            'natural_p_inverse_spectral_measures_have_positive_finite_mass_limits':True,
            'scaling_limit_unique':False,
            'full_continuum_field_algebra_constructed':False,
            'surviving_scattering_or_nonGaussianity_proved':False,
            'sharp_particle_pole_or_Lorentz_helicity_proved':False,
            'measured_particle_species_or_rest_mass_identified':False,
            'source_dynamical_diagonal_preserves_normalized_readouts':True,
            'effective_thermodynamic_ground_preparation_cost_proved':False,
            'finite_source_correlators_assumed_stationary':False,
            'native_pulses_executed_here':False,
            'all_source_repair_records_retained':True}


def direct_checks():
    return {'spatial':algebra.certificate(),'scaling':bounds.certificate(),'finite_control':experiment.check()}


def produce():
    pins=source_bindings()
    value={'schema':'oph.spin_particle_scaling.v1','theorem':theorem(),
           'direct_checks':direct_checks(),'inherited_custody':inherited(),
           'source_bindings':pins,'source_binding_scope':'explicit narrow mathematical inputs; not parent transitive closures',
           'new_checks_pass':True,'dependency_custody_pass':False,'all_pass':False}
    if source_bindings()!=pins:raise ValueError('source changed during replay')
    return value


def validate(value):
    if type(value) is not dict or value.get('schema')!='oph.spin_particle_scaling.v1':raise ValueError('invalid schema')
    if canonical(value)!=canonical(produce()):raise ValueError('payload, source, scope, or JSON type changed')


def cold():
    with tempfile.TemporaryDirectory(prefix='oph-spin-scaling-') as tmp:
        out=Path(tmp)/'receipt.json'
        env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
        subprocess.run([sys.executable,'-W','error','-m','model.spin_particle_scaling','--worker',str(out)],
                       cwd=ROOT,env=env,check=True,timeout=180)
        return strict_json.loads(out.read_text())


def main():
    parser=argparse.ArgumentParser();group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true');group.add_argument('--worker')
    args=parser.parse_args()
    if args.worker:
        Path(args.worker).write_text(canonical(produce()));return
    result=cold()
    if args.write:OUTPUT.write_text(canonical(result))
    elif canonical(strict_json.loads(OUTPUT.read_text()))!=canonical(result):raise ValueError('cold receipt changed')
    print(canonical({key:result[key] for key in ('new_checks_pass','dependency_custody_pass','all_pass')}).strip())
