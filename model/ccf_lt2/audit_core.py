#!/usr/bin/env python3
"""Independent finite charge/current/field bridge audits.

Inherited electromagnetic algebra and preparation checks for CCF-LT-2.
Algebraic residuals are checked exactly with integers or SymPy rationals.
Generic graph examples are algebra controls, not the one-edge transport run.
"""
from __future__ import annotations
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import argparse, itertools, json, math, random, hashlib
from collections import deque, Counter
from fractions import Fraction
from pathlib import Path
import numpy as np
import scipy.sparse as sp
import sympy as sy

HERE=Path(__file__).resolve().parent
# Exact committed oriented face table, LocalFaceMaxwellAction.lean, pinned in SOURCES.md.
ICO_FACES=[(0,1,2),(0,3,1),(0,2,4),(0,6,3),(0,4,6),(1,5,2),(1,3,7),(1,7,5),(2,8,4),(2,5,8),(3,6,9),(3,9,7),(4,10,6),(4,8,10),(5,7,11),(5,11,8),(6,10,9),(7,9,11),(8,11,10),(9,10,11)]

def complex_data(cone=False):
    edges=sorted({tuple(sorted((f[i],f[(i+1)%3]))) for f in ICO_FACES for i in range(3)})
    faces=ICO_FACES.copy(); nv=12
    if cone:
        faces += [(12,a,b) for a,b in edges]
        edges=sorted(edges+[(a,12) for a in range(12)])
        nv=13
    B=np.zeros((nv,len(edges)),dtype=np.int64)
    C=np.zeros((len(faces),len(edges)),dtype=np.int64)
    ei={e:i for i,e in enumerate(edges)}
    for k,(a,b) in enumerate(edges): B[a,k]=1;B[b,k]=-1
    for k,f in enumerate(faces):
        for a,b in zip(f,f[1:]+f[:1]): C[k,ei[tuple(sorted((a,b)))]]+=1 if a<b else -1
    assert not np.any(B@C.T)
    return nv,edges,faces,B,C

def charge_audit():
    """Derive charges from Lambda^0 + Lambda^2 + Lambda^4 C^5, not from a species lookup."""
    Y6=(-2,-2,-2,3,3); T32=(0,0,0,1,-1)
    rows=[]
    for degree in (0,2,4):
        for subset in itertools.combinations(range(5),degree):
            y6=sum(Y6[a] for a in subset); t32=sum(T32[a] for a in subset)
            q3=(3*t32+y6)//2
            assert (3*t32+y6)%2==0
            triality=sum(a<3 for a in subset)%3
            assert (q3+triality)%3==0
            rows.append(dict(subset=list(subset),degree=degree,Y6=y6,T3_twice=t32,Q3=q3,triality=triality))
    assert len(rows)==16
    counts=Counter(r['Q3'] for r in rows)
    assert counts==Counter({0:2,-2:3,2:3,-1:3,1:3,-3:1,3:1})
    # Full spinless 16-coordinate Fock census. Triality zero is necessary, not
    # sufficient for a color singlet. It is enough for the charge-integrality implication.
    zero_triality=0; no_half=True
    for mask in range(1<<16):
        q=sum(rows[i]['Q3'] for i in range(16) if mask>>i&1)
        t=sum(rows[i]['triality'] for i in range(16) if mask>>i&1)%3
        assert (q+t)%3==0
        if t==0:
            zero_triality+=1; assert q%3==0
        no_half &= (2*q != 3 and 2*q != -3)
    assert no_half
    return {'module':'Lambda^0 plus Lambda^2 plus Lambda^4 of C^5',
            'rows':rows,'Q3_multiplicities':dict(sorted(counts.items())),
            'fock_occupations_checked':1<<16,'triality_zero_occupations':zero_triality,
            'charge_integral_on_triality_zero':True,
            'scope':'Representation and Fock charge identity, not source selection of physical particles.'}

def bfs_tree(n,edges,root=0):
    adj=[[] for _ in range(n)]
    for e,(a,b) in enumerate(edges): adj[a].append((b,e));adj[b].append((a,e))
    parent=[-1]*n; pe=[-1]*n; depth=[0]*n;parent[root]=root;order=[root]
    for v in order:
        for u,e in sorted(adj[v]):
            if parent[u]<0: parent[u]=v;pe[u]=e;depth[u]=depth[v]+1;order.append(u)
    if len(order)!=n: raise ValueError('Graph must be connected')
    return parent,pe,depth,order

def gauss_repair(B,edges,q,electric,root=0,order=None):
    """Postorder local batched repair; q is protected. Stores each removed syndrome.
    B uses outward divergence. Units are integer native-third-charge units.
    """
    n,m=B.shape;parent,pe,depth,bfs=bfs_tree(n,edges,root)
    e=np.array([int(x) for x in electric],dtype=object);q=np.array([int(x) for x in q],dtype=object)
    g=B@e-q
    if int(sum(g))!=0: raise ValueError('Unbalanced total charge: retain a reference/boundary')
    agenda=list(reversed(bfs[1:])) if order is None else list(order)
    if sorted(agenda)!=[v for v in range(n) if v!=root]: raise ValueError('Agenda must visit every nonroot vertex exactly once')
    budget=(n-1)*sum(depth)
    history=[]; attempted=0; cost=lambda z:sum(depth[i]*abs(int(z[i])) for i in range(n))
    # Audit-only early termination; a fixed circuit can run the full proved budget.
    while np.any(g):
        changed=False
        for v in agenda:
            attempted+=1
            x=int(g[v]); old=cost(g)
            if not x: continue
            delta=-int(B[v,pe[v]])*x
            e[pe[v]]+=delta; g[v]=0;g[parent[v]]+=x
            assert cost(g)<=old-abs(x)
            assert np.array_equal(B@e-q,g)
            history.append((v,pe[v],x,delta));changed=True
        if not changed: raise AssertionError('False terminal state')
        if attempted>budget: raise AssertionError('Batched token bound violated')
    assert np.array_equal(B@e,q)
    return e,history,attempted

def topology_repair_audit():
    rng=random.Random(20261117); reports=[]; schedules_checked=0
    for cone in (False,True):
        n,edges,faces,B,C=complex_data(cone)
        Bs=sy.Matrix(B);Cs=sy.Matrix(C)
        rankB=Bs.rank();rankC=Cs.rank()
        assert rankB==n-1
        assert rankC==len(edges)-rankB
        L=Bs*Bs.T; J=sy.ones(n,n)/n;G=(L+J).inv()-J
        assert L*G==sy.eye(n)-J
        max_events=max_attempts=0;charge_kept=True;endpoint_same=True
        examples=[]
        for trial in range(128):
            q=np.array([rng.randint(-6,6) for _ in range(n)],dtype=np.int64)
            # Explicit compensating reference at the root, not deletion of charge.
            matter=q.copy();reference=-int(sum(q));q[0]+=reference
            e0=np.array([rng.randint(-4,4) for _ in edges],dtype=np.int64)
            _,_,_,bfs=bfs_tree(n,edges)
            finals=[]
            for agenda in (list(reversed(bfs[1:])),bfs[1:],rng.sample(list(range(1,n)),n-1)):
                e,h,attempts=gauss_repair(B,edges,q,e0,order=agenda)
                finals.append(e)
                assert len(h)<=sum(bfs_tree(n,edges)[2])
                max_events=max(max_events,len(h));max_attempts=max(max_attempts,attempts)
                schedules_checked+=1
            assert all(np.array_equal(finals[0],e) for e in finals[1:])
            # All chord fluxes are protected; tree edge solution is unique.
            tree_edges=set(bfs_tree(n,edges)[1][1:]); chords=[k for k in range(len(edges)) if k not in tree_edges]
            assert all(finals[0][k]==e0[k] for k in chords)
            if trial==0:
                qS=sy.Matrix(q); phi=G*qS;EL=Bs.T*phi
                assert Bs*EL==qS;assert Cs*EL==sy.zeros(len(faces),1)
                assert (EL.T*EL)[0]==(qS.T*G*qS)[0]
                examples.append({'matter_charge':matter.tolist(),'reference_charge':reference,
                    'physical_charge':q.tolist(),'final_flux':finals[0].tolist(),
                    'minimum_real_field_energy':str((EL.T*EL)[0]/2),
                    'minimum_real_field_is_integer':all(x.q==1 for x in EL)})
        # Exact Gauss charge/flux equality for every region of one chosen charged state.
        e=finals[0]
        for mask in range(1<<n):
            region=np.array([(mask>>i)&1 for i in range(n)],dtype=np.int64)
            assert int(region@q)==int((region@B)@e)
        reports.append({'carrier':'icosahedral_cone' if cone else 'icosahedral_screen',
            'vertices':n,'edges':len(edges),'faces':len(faces),'rank_divergence':rankB,'rank_curl':rankC,
            'random_inputs':128,'schedules_per_input':3,'maximum_accepted_repairs':max_events,
            'maximum_attempts':max_attempts,'regions_checked':1<<n,'examples':examples,
            'scope':'Cone is a declared extension; real Coulomb solution is not an integer compact-rotor basis state.'})
    # Charge-erasure integrity control: repair must preserve native charge.
    q=np.array([3,-3,0]);e=np.array([0,0]);B=np.array([[1,0],[-1,1],[0,-1]])
    fixed,_,_=gauss_repair(B,[(0,1),(1,2)],q,e)
    assert np.array_equal(q,np.array([3,-3,0])) and np.array_equal(fixed,np.array([3,0]))
    try:gauss_repair(B,[(0,1),(1,2)],[1,0,0],[0,0])
    except ValueError:unbalanced_rejected=True
    else:unbalanced_rejected=False
    assert unbalanced_rejected
    # Binary-encoded large charges use the same number of batched operations.
    big=1<<256
    bigE,big_history,big_attempts=gauss_repair(B,[(0,1),(1,2)],[big,-big,0],[0,0])
    assert list(bigE)==[big,0] and len(big_history)==1
    return {'reports':reports,'total_schedule_runs':schedules_checked,
            'charged_pair_survives_zero_residual':{'q':q.tolist(),'E':fixed.tolist(),'G':[0,0,0]},
            'unbalanced_sector_rejected':unbalanced_rejected,
            'large_integer_charge_bits':257,'large_integer_accepted_repairs':len(big_history)}

def csr_int(rows,cols,vals,d):
    x=sp.coo_matrix((np.array(vals,dtype=np.int64),(rows,cols)),shape=(d,d)).tocsr();x.sum_duplicates();x.eliminate_zeros();return x

def maxabs(M):
    M=M.tocsr();M.eliminate_zeros()
    if not np.issubdtype(M.dtype,np.integer): raise TypeError('Exact residual checker requires integer matrices')
    return 0 if M.nnz==0 else int(max(abs(M.data)))

def charge_field_kinematic(n=3,cap=3,q_unit=3,broken=False,wrap=False):
    """Complete kinematic Fock x hard-wall link space of a triangle. Exact integer H*16."""
    edges=[(0,1),(0,2),(1,2)];C=np.array([[1,-1,1]],dtype=np.int64)
    B=np.array([[1,1,0],[-1,0,1],[0,-1,-1]],dtype=np.int64)
    fluxes=list(itertools.product(range(-cap,cap+1),repeat=3))
    fidx={f:i for i,f in enumerate(fluxes)};nf=len(fluxes);dim=(1<<n)*nf
    Ediag=np.tile(np.array(fluxes,dtype=np.int64),(1<<n,1)).T
    background=np.array([0,1,0]); qdiag=np.zeros((n,dim),dtype=np.int64)
    for mask in range(1<<n):
        qdiag[:,mask*nf:(mask+1)*nf]=q_unit*(np.array([(mask>>v)&1 for v in range(n)])-background)[:,None]
    Gdiag=B@Ediag-qdiag
    def shiftflux(flux,deltas):
        out=tuple(int(x+y) for x,y in zip(flux,deltas))
        if wrap: return tuple(((x+cap)%(2*cap+1))-cap for x in out)
        return out if all(-cap<=x<=cap for x in out) else None
    T=[]
    for edge,(a,b) in enumerate(edges):
        rr=[];cc=[];vv=[]
        for mask in range(1<<n):
            if not(mask>>a&1) or mask>>b&1:continue
            after=mask^(1<<a);sgn=(-1)**((mask&((1<<a)-1)).bit_count()+ (after&((1<<b)-1)).bit_count())
            after|=1<<b
            for k,f in enumerate(fluxes):
                ds=[0]*3;ds[edge]=0 if broken else -q_unit
                fp=shiftflux(f,ds)
                if fp is None:continue
                rr.append(after*nf+fidx[fp]);cc.append(mask*nf+k);vv.append(sgn)
        T.append(csr_int(rr,cc,vv,dim))
    rr=[];cc=[];vv=[]
    for mask in range(1<<n):
        for k,f in enumerate(fluxes):
            fp=shiftflux(f,C[0]);
            if fp is not None:rr.append(mask*nf+fidx[fp]);cc.append(mask*nf+k);vv.append(1)
    P=csr_int(rr,cc,vv,dim)
    # H=mass + (1/16) E^2 + (1/4)(1-Re P) -(1/8)sum(T+T*)
    diag=np.sum(Ediag**2,axis=0)+4
    for mask in range(1<<n):
        diag[mask*nf:(mask+1)*nf]+=16*sum((1 if v!=1 else -1)*((mask>>v)&1) for v in range(n))
    H=sp.diags(diag,format='csr',dtype=np.int64)-2*(P+P.T)
    for t in T:H=H-2*(t+t.T)
    H=H.astype(np.int64)
    assert maxabs(H-H.T)==0
    # J = i/16 Jint, Jint=2q(T-T*)
    Jint=[2*q_unit*(t-t.T) for t in T]
    comm_diagonal=lambda X,d: X@sp.diags(d,dtype=np.int64)-sp.diags(d,dtype=np.int64)@X
    gauss=[maxabs(comm_diagonal(H,g)) for g in Gdiag]
    continuity=[];ampere=[]
    for v in range(n):
        rhs=sp.csr_matrix(H.shape,dtype=np.int64)
        for k in range(3):rhs=rhs-int(B[v,k])*Jint[k]
        continuity.append(maxabs(comm_diagonal(H,qdiag[v])-rhs))
    for k in range(3):
        # magnetic H intensity /i = 2*(P-P*)/16; outward E convention.
        rhs=2*int(C[0,k])*(P-P.T)-Jint[k]
        ampere.append(maxabs(comm_diagonal(H,Ediag[k])-rhs))
    Qsum=np.sum(qdiag,axis=0)
    physical=np.all(Gdiag==0,axis=0)
    # Complete constraint and charge measurements commute with the constraint;
    # duplicated orthogonal labels do not copy unknown private amplitudes.
    phys_count=int(physical.sum())
    leakage=maxabs(H[np.flatnonzero(~physical)][:,np.flatnonzero(physical)])
    return {'dimension':dim,'cap':cap,'q_unit':q_unit,'H_integer_denominator':16,
            'nonzero_H_entries':int(H.nnz),'physical_dimension':phys_count,'physical_kernel_leakage_integer':leakage,
            'max_Gauss_commutator':max(gauss),'continuity_residual':max(continuity),
            'Ampere_residual':max(ampere),'global_charge_residual':maxabs(comm_diagonal(H,Qsum)),
            'mutation':'omit_matter_flux_shift' if broken else ('cyclic_flux_wrap' if wrap else 'none')}

def quantum_audit():
    good=[charge_field_kinematic(cap=3,q_unit=q) for q in (1,3)]
    for r in good:
        assert r['max_Gauss_commutator']==r['continuity_residual']==r['Ampere_residual']==r['global_charge_residual']==r['physical_kernel_leakage_integer']==0
    bad1=charge_field_kinematic(cap=3,q_unit=3,broken=True)
    bad2=charge_field_kinematic(cap=3,q_unit=3,wrap=True)
    assert bad1['max_Gauss_commutator']>0 and bad2['max_Gauss_commutator']>0
    assert bad1['physical_kernel_leakage_integer']>0 and bad2['physical_kernel_leakage_integer']>0
    return {'complete_kinematic_checks':good,'negative_controls':[bad1,bad2],
            'scope':'Diagnostic local charged coordinates. Not the full native 192-mode numerical Hamiltonian.'}

def annihilator(n,j):
    A=sy.zeros(1<<n)
    for mask in range(1<<n):
        if mask>>j&1:A[mask^(1<<j),mask]=(-1)**((mask&((1<<j)-1)).bit_count())
    return A

def fock_audit():
    a=[annihilator(4,j) for j in range(4)];ad=[x.T for x in a]
    c=[(a[0]-ad[1])/sy.sqrt(2),(a[1]+ad[0])/sy.sqrt(2),
       (a[2]-ad[3])/sy.sqrt(2),(a[3]+ad[2])/sy.sqrt(2)]
    for i in range(4):
        for j in range(4):
            assert sy.simplify(c[i]*c[j].T+c[j].T*c[i])==(sy.eye(16) if i==j else sy.zeros(16))
    q=(-3,3,-3,3)
    Q=sum((q[i]*ad[i]*a[i] for i in range(4)),sy.zeros(16))
    Qc=sum((q[i]*c[i].T*c[i] for i in range(4)),sy.zeros(16))
    assert sy.simplify(Qc-Q)==sy.zeros(16)
    # Literal native bilinear hop/current transformed before any number restriction.
    hop=c[2].T*c[0]; J=sy.I*q[0]*(hop-hop.T)
    Ndiag=[m.bit_count() for m in range(16)]
    pair=sy.zeros(16)
    for i in range(16):
        for j in range(16):
            if Ndiag[i]!=Ndiag[j]:pair[i,j]=J[i,j]
    norm2=sy.simplify(sy.trace(pair.conjugate().T*pair))
    assert norm2>0;assert sy.simplify(Q*J-J*Q)==sy.zeros(16)
    # Native source field-dependent frame warning, independent 4-state control at theta=0.
    # U(theta) even block [[c,s exp(iq theta)],[-s exp(-iq theta),c]], odd=I.
    cc=sy.Rational(4,5);ss=sy.Rational(3,5);qq=3
    U=sy.eye(4);Up=sy.zeros(4)
    U[0,0]=U[3,3]=cc;U[0,3]=ss;U[3,0]=-ss
    Up[0,3]=sy.I*qq*ss;Up[3,0]=sy.I*qq*ss
    conn=sy.simplify(-sy.I*U.T*Up)
    assert conn==conn.conjugate().T
    pair_amp2=sy.simplify(abs(conn[3,0])**2)
    assert pair_amp2==sy.Rational(1296,625)
    # A simple complete physical-code charge meter on two logical labels,
    # demonstrating preservation of coherence before, but not after charge readout.
    V=sy.zeros(6,2);V[1,0]=1;V[4,1]=1
    assert V.T*V==sy.eye(2)
    rho=sy.Matrix([[1,-sy.I],[sy.I,1]])/2
    assert V.T*(V*rho*V.T)*V==rho
    return {'CAR_pairs_checked':16,'charge_transform_exact':True,
            'full_Fock_dimension':16,'pair_current_squared_Frobenius_norm':str(norm2),
            'field_connection_pair_amplitude_squared':str(pair_amp2),
            'scope':'Exact finite operator checks; source pair coefficients retained, not discarded.'}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('task',choices=['charge','topology','quantum','fock','all'])
    parser.add_argument('--output',type=Path);args=parser.parse_args()
    jobs={'charge':charge_audit,'topology':topology_repair_audit,'quantum':quantum_audit,'fock':fock_audit}
    if args.task=='all':out={k:f() for k,f in jobs.items()}
    else:out=jobs[args.task]()
    text=json.dumps(out,indent=2,sort_keys=True)
    if args.output:args.output.write_text(text+'\n')
    print(text)
if __name__=='__main__':main()
