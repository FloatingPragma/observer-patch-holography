"""Actual weak current of the new uncalibrated colour current-square source.

Independent numerical vectors are enclosed using exact rational inertia and
radical residual bounds. No calibrated hadron arrays are reused.
"""
from fractions import Fraction as F
from pathlib import Path
import importlib.util
import hashlib
import json
import numpy as np

HERE=Path(__file__).resolve().parent
HAD=HERE.parent/'hadrons'
spec=importlib.util.spec_from_file_location('oph_repair_hadron_audit',HAD/'audit.py')
had=importlib.util.module_from_spec(spec);spec.loader.exec_module(had)
src=had.source


def add(a,b):return a[0]+b[0],a[1]+b[1]
def scale(a,q):return (a[0]*q,a[1]*q) if q>=0 else (a[1]*q,a[0]*q)
def sup(a):return max(abs(a[0]),abs(a[1]))


def kth(A,G,k):
    def count(x):
        B=[row[:] for row in A]
        for i,g in enumerate(G):B[i][i]-=x*g
        return had.inertia(B)
    lo,hi=F(0),F(16)
    assert count(lo)<k<=count(hi)
    while hi-lo>F(1,1<<36):
        mid=(lo+hi)/2
        if count(mid)<k:lo=mid
        else:hi=mid
    return lo,hi


def weak_operator(dn, target_counts, cell, mu):
    """sigma_mu times weak raising: down -> up, full occupation amplitudes."""
    idx={x:i for i,x in enumerate(target_counts)}
    P=(np.eye(2),np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1,-1]))[mu]
    M=np.zeros((len(idx),dn['dimension']),complex)
    exact={}
    for col,occ in enumerate(dn['counts']):
        for spin_out in range(2):
            for spin_in in range(2):
                coef=P[spin_out,spin_in]
                if coef==0:continue
                target,z=src.act(occ,4*cell+2*spin_out,4*cell+2*spin_in+1)
                if target is not None and target in idx:
                    row=idx[target];M[row,col]+=coef*src.value(z)
                    if mu in (0,3):exact[row,col]=src.add(exact.get((row,col),{}),z,F(int(coef.real)))
    return M,exact


def spin_raise(d):
    target_counts=[]
    for occ in d['counts']:
        for x in range(d['L']):
            for w in range(2):
                t,z=src.act(occ,4*x+w,4*x+2+w)
                if t is not None and t not in target_counts:target_counts.append(t)
    idx={x:i for i,x in enumerate(target_counts)};op={}
    for col,occ in enumerate(d['counts']):
        for x in range(d['L']):
            for w in range(2):
                t,z=src.act(occ,4*x+w,4*x+2+w)
                if t is not None:op[idx[t],col]=src.add(op.get((idx[t],col),{}),z)
    return target_counts,op


def vector(weakup):
    d=src.data(2,F(1),weakup);vals=had.extension_values(d,'current_square')
    H=src.dense(d)[:d['dimension'],:d['dimension']]+np.diag([float(v) for v in vals[:d['dimension']]])
    values,U=np.linalg.eigh(H);q=[F(f'{v:.14f}') for v in U[:,0]]
    A,G=had.rational_matrix(d,vals);ground=kth(A,G,1);second=kth(A,G,2)
    energy=sum(ground,F())/2;gap=second[0]-ground[1];assert gap>0
    residual=F()
    for i in range(d['dimension']):
        r=(F(),F())
        for j in range(d['dimension']):
            entry=src.interval(d['entries'].get((i,j),{}),80)
            if i==j:entry=add(entry,(vals[i]-energy,vals[i]-energy))
            r=add(r,scale(entry,q[j]))
        residual+=sup(r)
    norm2=sum(x*x for x in q);assert F(999,1000)**2<norm2<F(1001,1000)**2
    # normalized trial residual, plus the ground interval's half width
    delta=residual/F(999,1000)
    distance=2*delta/gap+abs(norm2-1)
    assert distance<F(1,100000)
    targets,raising=spin_raise(d)
    splus=F()
    for i in range(len(targets)):
        r=(F(),F())
        for j in range(d['dimension']):r=add(r,scale(src.interval(raising.get((i,j),{}),80),q[j]))
        splus+=sup(r)
    assert splus+3*distance<F(1,10)
    lookup={c:i for i,c in enumerate(d['counts'])}
    reflected=[lookup[c[4:]+c[:4]] for c in d['counts']]
    reflection_error=sum(abs(q[reflected[i]]-q[i]) for i in range(len(q)))
    assert reflection_error+2*distance<F(1,10)
    # [H,S²]=0 follows from each source term being spin scalar. The unique
    # ground cannot have S=3/2: then ||S+ psi||=sqrt3, contradicting this bound.
    return d,H,np.array(list(map(float,q))),q,dict(weakup=weakup,
        ground_base_interval=list(map(str,ground)),second_base_interval=list(map(str,second)),
        spectral_gap_lower=str(gap),rational_trial_norm2=str(norm2),
        exact_radical_residual_l1_upper=str(residual),ground_vector_distance_upper=str(distance),
        spin_raising_norm_upper=str(splus+3*distance),ground_spin='1/2',
        cell_reflection_even_error_upper=str(reflection_error+2*distance),
        exact_cell_reflection_parity=1),distance


def sandwich(op,qp,qn):
    val=(F(),F())
    for (i,j),z in op.items():val=add(val,scale(src.interval(z,80),qp[i]*qn[j]))
    return val


def main():
    inputs=[HAD/'audit.py',HAD/'inputs/0_source.py',HAD/'RESULT.json',HAD/'inputs/3_results.json']
    pins={str(p.relative_to(HERE.parents[3])):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    dn,Hn,n,qn,rn,en=vector(1);dp,Hp,p,qp,rp,ep=vector(2)
    exact_v={};exact_a={}
    for x in range(2):
        _,vv=weak_operator(dn,dp['counts'],x,0);_,aa=weak_operator(dn,dp['counts'],x,3)
        for key,z in vv.items():exact_v[key]=src.add(exact_v.get(key,{}),z)
        for key,z in aa.items():exact_a[key]=src.add(exact_a.get(key,{}),z)
    v=sandwich(exact_v,qp,qn);a=sandwich(exact_a,qp,qn)
    if sum(v)<0:
        p=-p;qp=[-q for q in qp];v=(-v[1],-v[0]);a=(-a[1],-a[0])
    error=3*(ep+F(1001,1000)*en)
    v=(v[0]-error,v[1]+error);a=(a[0]-error,a[1]+error)
    assert v[0]>0
    ga=(a[0]/v[1],a[1]/v[0]) if a[0]>=0 else (a[0]/v[0],a[1]/v[1])
    # Build the proton spin partner by its actual total lowering operator.
    down=[tuple(occ[i^2] for i in range(8)) for occ in dp['counts']];idx={c:i for i,c in enumerate(down)}
    lower=np.zeros((len(down),dp['dimension']))
    for j,occ in enumerate(dp['counts']):
        for x in range(2):
            for w in range(2):
                target,z=src.act(occ,4*x+2+w,4*x+w)
                if target is not None:lower[idx[target],j]+=src.value(z)
    pdown=lower@p
    local=np.zeros((2,4,2),complex)
    for x in range(2):
        for mu in range(4):
            local[x,mu,0]=np.vdot(p,weak_operator(dn,dp['counts'],x,mu)[0]@n)
            local[x,mu,1]=np.vdot(pdown,weak_operator(dn,down,x,mu)[0]@n)
    fv=local[:,0,0].sum().real;fa=local[:,3,0].sum().real
    P=[np.eye(2),np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1,-1])]
    expected=np.array([[.5*(fv if mu==0 else fa)*P[mu][:,0] for mu in range(4)] for _ in range(2)])
    matrix_error=float(np.max(abs(local-expected)))
    assert matrix_error<2e-10 and abs(np.vdot(pdown,pdown)-1)<2e-10
    gap=json.loads((HAD/'RESULT.json').read_text())['gap_intervals_exact']['current_square']
    out=dict(schema='oph.current_square_weak_quartet.v1',input_sha256=pins,
        source='H_colour + vectorlike number + sum_x(Bc-Qc)^2; U=1, kappa=1, L=2',
        calibrated_hadron_arrays_reused=False,source_gap_exact=gap,
        ground_certificates=[rn,rp],vector_interval=list(map(str,v)),axial_interval=list(map(str,a)),
        axial_over_vector_interval=list(map(str,ga)),vector_float=fv,axial_float=fa,ratio_float=fa/fv,
        local_current_real=local.real.tolist(),local_current_imag=local.imag.tolist(),
        local_quartet_numeric_error=matrix_error,
        exact_scope='Ground-vector enclosures and vector/axial matrix intervals are exact rational/radical bounds. Spin-half follows from spin invariance, uniqueness and raising-norm exclusion. Local quartet factors also follow from spin covariance and cell reflection; matrices are independently checked numerically.',
        remaining=['full native amended ground preparation/filter','common lepton/color clock normalization','3D weak evolution','physical hadron masses and splitting','physical lifetime'])
    for p in inputs:assert hashlib.sha256(p.read_bytes()).hexdigest()==pins[str(p.relative_to(HERE.parents[3]))]
    (HERE/'AMENDED_CURRENT.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ('vector_float','axial_float','ratio_float','local_quartet_numeric_error','source_gap_exact')},indent=2))


if __name__=='__main__':main()
