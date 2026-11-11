"""Reproduce finite blocking records and conditional matrix-running controls."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import sympy as sp

from . import threshold, running

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def bindings():
    paths=set(HERE.glob('*.py'))|set(HERE.glob('*.md'))
    base=ROOT/'model/interacting_extension/higgs_join'
    paths.update(base/name for name in ('interacting_carrier.py','run_joined_response.py',
        'joined_operators.npz','joined_exact_operators.json','joined_basis.json','joined_operator_receipt.json'))
    paths.update(ROOT/name for name in ('model/interacting_extension/dependencies/source_reference.py',
        'model/interacting_extension/dependencies/exact_family_operators.json',
        'model/observer_particles/results/unified_yukawa_operators.json'))
    return {p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}


def report():
    from ..interacting_extension.dependencies.source_reference import exact_source_packet
    from ..interacting_extension.higgs_join.interacting_carrier import primary_family_matrices
    source=exact_source_packet(ROOT)
    y,sigma=primary_family_matrices()
    if np.max(abs(y-np.asarray(source['Y'],complex))) > 1e-14:
        raise ArithmeticError('primary family coefficient differs from complete source replay')
    exact={'literal_CAR':threshold.exact_local_car(),
           'complete_source_operators':threshold.exact_source_commutators(),
           'source_spectral_bounds':threshold.source_spectral_bounds(source['Y'])}
    before=bindings()
    arrays={'source_Y':y,'source_sigma':sigma}
    cases={}
    catalog=[('native_no_hopping',y,sigma,1.,0.,0.),
             ('native_half_Higgs',y,sigma,.5,0.,0.),
             ('native_opposite_Higgs',y,sigma,-1.,0.,0.),
             ('native_reference_shift',y,sigma,1.,0.,.125),
             ('zero_vertex',np.zeros((3,3)),np.zeros((3,3)),1.,0.,0.),
             ('symmetric_family',np.eye(3),np.eye(3),1.,0.,0.),
             ('native_quarter_hopping',y,sigma,1.,.25,0.),
             ('native_half_hopping',y,sigma,1.,.5,0.),
             ('native_full_hopping',y,sigma,1.,1.,0.),
             ('native_full_hopping_quarter_Higgs',y,sigma,.25,1.,0.)]
    for name,yy,ss,g,h,r in catalog:
        result,aa=threshold.complete_control(yy,ss,coupling=g,hopping=h,reference=r)
        cases[name]=result;arrays.update({name+'__'+k:v for k,v in aa.items()})
        print('Computed complete finite record:',name,flush=True)
    # A complex positive matrix controls all-left conjugation and arbitrary
    # family frames independently of the real source example.
    u=np.array([[1,1j,0],[1j,1,0],[0,0,np.sqrt(2)]])/np.sqrt(2)
    complex_y=u@np.diag([.125,.25,.5])@u.conj().T
    complex_result,aa=threshold.complete_control(complex_y,np.eye(3))
    cases['complex_family_control']=complex_result
    arrays.update({'complex_family_control__'+k:v for k,v in aa.items()})
    rg,ra=running.evolve(np.asarray([y,sigma,y]))
    arrays.update({'running__'+k:v for k,v in ra.items()})
    if bindings()!=before:
        raise ArithmeticError('consumed implementation changed during calculation')
    result={'schema':'oph.yukawa_source_running.v1','implementation_sha256':before,
        'source_sha256':source['source_hashes'],'exact':exact,'finite_cases':cases,
        'conditional_running':rg,'running_controls':running.algebra_controls(),
        'source_family_singular_values':np.linalg.svd(y,compute_uv=False)[::-1].tolist(),
        'status':{'empirical_particle_inputs':[], 'source_Y_reconstructed_from_complete_records':True,
          'source_selected_Higgs_coupling_scalar_energy_or_hopping':False,
          'ground_preparation_and_reference_energy_are_declared':True,
          'no_hopping_threshold_formula_is_exact_on_declared_finite_source':True,
          'nonzero_hopping_controls_are_complete_numeric_diagonalizations':True,
          'numeric_errors_are_rigorous_interval_bounds':False,
          'conditional_minimal_SM_one_loop_flow':True,
          'source_to_running_Yukawa_matching_derived':False,
          'physical_mass_or_continuum_prediction':False}}
    return result,arrays


def _close(a,b,path='root'):
    if isinstance(a,dict):
        if not isinstance(b,dict) or a.keys()!=b.keys():raise ArithmeticError('receipt keys '+path)
        for key in a:_close(a[key],b[key],path+'/'+key)
    elif isinstance(a,list):
        if not isinstance(b,list) or len(a)!=len(b):raise ArithmeticError('receipt list '+path)
        for i,(x,y) in enumerate(zip(a,b)):_close(x,y,path+'/'+str(i))
    elif isinstance(a,float):
        if not np.isclose(a,b,rtol=1e-7,atol=3e-10):raise ArithmeticError('receipt number '+path)
    elif a!=b:raise ArithmeticError('receipt exact '+path)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true')
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result,arrays=report()
    if args.write:
        np.savez_compressed(HERE/'spectral_records.npz',**arrays)
        (HERE/'results.json').write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
    if args.check:
        _close(result,json.loads((HERE/'results.json').read_text()))
        with np.load(HERE/'spectral_records.npz',allow_pickle=False) as saved:
            if set(saved.files)!=set(arrays):raise ArithmeticError('complete array inventory')
            for key,value in arrays.items():
                # Amplitudes and created vectors depend on arbitrary numerical
                # eigenvector phases. Compare invariant residues and record laws.
                if 'amplitudes' in key or 'created_vectors' in key:continue
                if saved[key].shape!=value.shape or not np.allclose(saved[key],value,atol=3e-9,rtol=1e-7):
                    raise ArithmeticError('complete array '+key)
    print('Verified finite Higgs blocking and conditional running.')


if __name__=='__main__':main()
