"""Exact-input high-precision refinement of the lowest addition lines.

Sparse matrix coefficients are rebuilt from integer CAR/Gauss numerators,
exact basis norms, exact native family expressions and log(5/4).  A double
precision bordered inverse is used only as a preconditioner; every residual,
normalization, Rayleigh quotient and update is evaluated in mpmath.  Agreement
between 50 and 80 digits is a precision audit, not interval arithmetic.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import numpy as np
import sympy as sp
from scipy.linalg import lu_factor,lu_solve
from .interacting_carrier import HERE,build_join


def mp_scalar(text,dps):
    value=sp.sympify(text)
    a,b=sp.expand_complex(value).as_real_imag()
    return mp.mpc(str(sp.N(a,dps)),str(sp.N(b,dps)))


def rebuild_entries(root,reference_file,q,dps,g=1):
    exact=json.loads((root/'joined_exact_operators.json').read_text())['operators']
    basis=json.loads((root/'joined_basis.json').read_text())[str(q)]
    refs=json.loads(reference_file.read_text())
    family=lambda key:[[mp_scalar(v,dps) for v in row] for row in refs[key]]
    y,d,j=family('Y'),family('D'),family('J');phi=mp.log(mp.mpf(5)/4)
    sigma=[[y[a][b]+(1-mp.cos(phi))*d[a][b]+mp.sin(phi)*j[a][b] for b in range(3)] for a in range(3)]
    norms=[mp_scalar(x,dps).real for x in basis['raw_norms']]
    entries={}
    def add(i,k,value):
        entries[i,k]=entries.get((i,k),mp.mpc(0))+value
    diagonal=exact[f'q{q}_diagonal']
    for i,(electric,nh) in enumerate(zip(diagonal['electric_native_trace'],diagonal['N_H'],strict=True)):
        add(i,i,mp_scalar(electric,dps)+mp_scalar(nh,dps))
    def hermitian_term(name,coefficient):
        if not coefficient:return
        for i,k,value in exact[name]:
            a=coefficient*mp_scalar(value,dps)/mp.sqrt(norms[i]*norms[k])
            add(i,k,a);add(k,i,mp.conj(a))
    hermitian_term(f'q{q}_loop',-1)
    for a in range(3):
        for b in range(3):
            hermitian_term(f'q{q}_D0_{a}_{b}',-y[a][b])
            hermitian_term(f'q{q}_D1_{a}_{b}',-sigma[a][b])
            hermitian_term(f'q{q}_W_{a}_{b}',g*y[a][b])
    return {key:value for key,value in entries.items() if value},basis


def matrix_rounding_error(entries,indices,float_h):
    """Frobenius difference on the union of exact and floating supports."""
    inverse={int(old):new for new,old in enumerate(indices)}
    if float_h.shape!=(len(indices),len(indices)):
        raise ValueError('matrix and index dimensions differ')
    exact={(inverse[a],inverse[b]):value for (a,b),value in entries.items()
           if a in inverse and b in inverse}
    support=set(exact)|set(zip(*np.nonzero(float_h)))
    return mp.sqrt(mp.fsum(abs(exact.get((i,j),0)-
        mp.mpc(float(float_h[i,j].real),float(float_h[i,j].imag)))**2
        for i,j in support))


def refine(entries,indices,float_h,float_energy,float_vector,dps):
    inverse={old:new for new,old in enumerate(indices)};n=len(indices)
    rows=[[] for _ in range(n)]
    for (a,b),value in entries.items():
        if a in inverse:
            if b not in inverse:raise ArithmeticError('exact spin block leakage')
            rows[inverse[a]].append((inverse[b],value))
    def matvec(x):return [mp.fsum(value*x[j] for j,value in row) for row in rows]
    def normalize(x):
        length=mp.sqrt(mp.fsum(abs(z)**2 for z in x))
        return [z/length for z in x]
    x=normalize([mp.mpc(float(z.real),float(z.imag)) for z in float_vector])
    border=np.zeros((n+1,n+1),complex)
    border[:n,:n]=float_h-float_energy*np.eye(n)
    border[:n,n]=float_vector;border[n,:n]=float_vector.conj()
    lu=lu_factor(border)
    ladder=[]
    for iteration in range(12):
        hx=matvec(x);rho=mp.re(mp.fsum(mp.conj(z)*w for z,w in zip(x,hx,strict=True)))
        residual=[w-rho*z for z,w in zip(x,hx,strict=True)]
        residual_norm=mp.sqrt(mp.fsum(abs(z)**2 for z in residual))
        ladder.append(mp.nstr(residual_norm,8))
        if residual_norm<mp.power(10,-dps+12):break
        scale=max(abs(z) for z in residual)
        rhs=np.zeros(n+1,complex);rhs[:n]=[-complex(z/scale) for z in residual]
        delta=lu_solve(lu,rhs)[:n]
        x=normalize([z+scale*mp.mpc(float(w.real),float(w.imag)) for z,w in zip(x,delta,strict=True)])
    else:raise ArithmeticError('high precision eigenvector refinement did not converge')
    # Full norm of exact-H minus float-H on this block supplies a simple
    # matrix-rounding diagnostic independent of the eigensolver residual.
    return rho,x,{'working_decimal_digits':dps,'dimension':n,
                  'iterations':len(ladder)-1,'residual_norm_ladder':ladder,
                  'final_eigen_residual_norm':mp.nstr(residual_norm,25),
                  'Rayleigh_energy':mp.nstr(rho,dps-10),
                  'exact_to_float_matrix_Frobenius_difference':mp.nstr(matrix_rounding_error(entries,indices,float_h),20)}


def _run_at_precision(root,reference_file,response_root=None,output=None):
    start=time.monotonic();root=Path(root);reference_file=Path(reference_file)
    response_root=Path(response_root) if response_root is not None else root
    output=Path(output) if output is not None else response_root
    output.mkdir(parents=True,exist_ok=True)
    with np.load(response_root/'joined_response_arrays.npz',allow_pickle=False) as a:
        stored={k:a[k].copy() for k in a.files if k.startswith('dressed_q') or k in ('source_Y','source_sigma')}
    joined=build_join(root,family0=stored['source_Y'],family1=stored['source_sigma'],
                      yukawa=stored['source_Y'],fermion_energy=0.,scalar_energy=1.)
    targets=[(0,0),(1,0),(1,2),(1,4)]
    all_tables=[];final_vectors={};summaries=[]
    for dps in (50,80):
        with mp.workdps(dps):
            table={};entries={q:rebuild_entries(root,reference_file,q,dps) for q in (0,1)}
            for q,col in targets:
                h=joined['sectors'][q]['H'];obs=joined['sectors'][q]['observables']
                weight=int(stored[f'dressed_q{q}_twice_Sz'][col])
                indices=np.flatnonzero(obs['twice_Sz'].diagonal().real.astype(int)==weight)
                energy=float(stored[f'dressed_q{q}_energies'][col]);vector=stored[f'dressed_q{q}_vectors'][indices,col]
                rho,x,receipt=refine(entries[q][0],indices,h[indices][:,indices].toarray(),energy,vector,dps)
                table[q,col]=rho;receipt.update(Qstar=q,eigenstate_index=col,twice_Sz=weight)
                summaries.append(receipt)
                if dps==80:
                    full=np.zeros(h.shape[0],complex);full[indices]=[complex(z) for z in x]
                    final_vectors[f'q{q}_state{col}']=full
                print(json.dumps({'stage':'precision_refinement','dps':dps,'Qstar':q,'column':col,
                                  'residual':receipt['final_eigen_residual_norm'],
                                  'seconds':time.monotonic()-start}),flush=True)
            all_tables.append(table)
    disagreement=max(abs(all_tables[0][k]-all_tables[1][k]) for k in targets)
    if disagreement>mp.mpf('1e-40'):raise ArithmeticError('50 to80 digit precision ladder')
    gaps=[all_tables[1][1,k]-all_tables[1][0,0] for k in (0,2,4)]
    float_gaps=[float(stored['dressed_q1_energies'][k]-stored['dressed_q0_energies'][0]) for k in (0,2,4)]
    if not all(x>0 for x in gaps):raise ArithmeticError('lowest dressed addition lines changed sign')
    result={'schema':'oph.interacting.lepton_higgs_precision.v1',
            'exact_native_family_expressions':str(reference_file.resolve()),
            'exact_native_family_file_sha256':hashlib.sha256(reference_file.read_bytes()).hexdigest(),
            'integer_CAR_and_exact_Gauss_normalization_rebuilt':True,
            'history_angle_exact':'log(5/4)','controls':'same epsilon_F=0, epsilon_H=1 and all dynamical strengths1',
            'refinement_method':'mpmath sparse residual with float64 bordered preconditioner',
            'interval_arithmetic_or_rigorous_floating_enclosure_claimed':False,
            'all_target_energies_at50_and80_digits':summaries,
            'maximum50_to80_digit_energy_difference':mp.nstr(disagreement,25),
            'lowest_three_distinct_addition_lines_40_digits':[mp.nstr(x,40) for x in gaps],
            'corresponding_float64_subtractions':float_gaps,
            'absolute_float64_addition_line_errors':[mp.nstr(abs(x-mp.mpf(y)),20) for x,y in zip(gaps,float_gaps,strict=True)],
            'all_other_eigenstates_retained_in_main_response_packet':True,
            'particle_mass_claimed':False,'elapsed_seconds':time.monotonic()-start}
    (output/'joined_precision_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    np.savez_compressed(output/'joined_refined_low_vectors.npz',**final_vectors,
                        addition_lines=np.asarray([float(x) for x in gaps]))
    return result


def run(root,reference_file,response_root=None,output=None):
    """Refine at controlled precision without changing the caller's context."""
    with mp.workdps(80):
        return _run_at_precision(root,reference_file,response_root,output)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=HERE)
    p.add_argument('--references',type=Path,default=HERE.parent/'dependencies'/'exact_family_operators.json')
    p.add_argument('--response-root',type=Path,required=True)
    p.add_argument('--output',type=Path)
    args=p.parse_args();run(args.root,args.references,args.response_root,args.output)
