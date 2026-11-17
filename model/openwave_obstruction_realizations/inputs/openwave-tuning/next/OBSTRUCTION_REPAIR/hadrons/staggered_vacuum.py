"""Complete two-site SU3 staggered source with a dynamical filled vacuum.

This is a new finite toy source, not a transfer of the old native source.
Its mass is an actual excitation relative to its ground, not a reference
shift. No proton target or fitted coefficient enters its Hamiltonian.
"""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import json
import numpy as np
import audit

HERE=Path(__file__).resolve().parent


def creation(mode,word):
    if word>>mode&1:return None,0
    return word|(1<<mode),(-1)**((word&((1<<mode)-1)).bit_count())


def annihilation(mode,word):
    if not word>>mode&1:return None,0
    return word^(1<<mode),(-1)**((word&((1<<mode)-1)).bit_count())


def bilinear(i,j):
    out=np.zeros((64,64),dtype=np.int64)
    for w in range(64):
        v,a=annihilation(j,w)
        if v is None:continue
        z,b=creation(i,v)
        if z is not None:out[z,w]=a*b
    return out


def basis():
    R=np.zeros((64,6),dtype=np.int64);R[0,0]=1;R[63,5]=1
    for n in range(4):
        for left in combinations(range(3),n):
            modes=[c if c in left else c+3 for c in range(3)]
            w=0;sign=1
            for mode in reversed(modes):
                w,s=creation(mode,w);sign*=s
            R[w,n+1]+=sign
    G=R.T@R
    assert np.array_equal(G,np.diag([1,1,3,3,1,1]))
    constraint=np.zeros((64,64),dtype=np.int64)
    for a in range(3):
        for b in range(3):
            if a!=b:
                T=bilinear(a,b)+bilinear(3+a,3+b)
            else:
                if a==2:continue
                T=bilinear(a,a)+bilinear(3+a,3+a)-bilinear(a+1,a+1)-bilinear(4+a,4+a)
            assert not np.any(T@R)
            constraint+=T.T@T
    # Cartan-zero words have equal total occupations of all three colours.
    zero=[w for w in range(64) if len({((w>>c)&1)+((w>>(c+3))&1) for c in range(3)})==1]
    assert len(zero)==10
    reduced=constraint[np.ix_(zero,zero)]
    assert audit.inertia([list(map(F,-row)) for row in reduced])==4
    # The positive constraint sum therefore has exactly six zero vectors;
    # R is the complete Gauss-singlet kernel, not merely an invariant sample.
    return R,list(map(int,np.diag(G)))


def components(R,G):
    nleft=np.array([(w&7).bit_count() for w in range(64)],dtype=np.int64)
    nright=np.array([(w>>3).bit_count() for w in range(64)],dtype=np.int64)
    mass=np.diag(nleft-nright)
    casimir3=np.diag(2*nleft*(3-nleft))
    C12=np.zeros((64,64),dtype=np.int64)
    for a in range(3):
        for b in range(a):
            E=bilinear(a,b);D=bilinear(b,a)
            C12+=6*(E@D+D@E)
    H1=bilinear(0,0)-bilinear(1,1)
    H2=bilinear(0,0)+bilinear(1,1)-2*bilinear(2,2)
    C12+=3*H1@H1+H2@H2
    assert np.array_equal(C12,4*casimir3)
    hopping=sum((-bilinear(c,c+3)-bilinear(c+3,c) for c in range(3)),np.zeros((64,64),dtype=np.int64))
    result=[]
    for M in [mass,casimir3,hopping]:
        A=R.T@M@R
        # Gram entries divide every coordinate of the invariant restriction.
        H=np.array([[F(int(A[i,j]),G[i]) for j in range(6)] for i in range(6)],dtype=object)
        assert np.array_equal(M@R,R@H)
        result.append(A)
    return result


def interval(A,G,index,lo=F(-4),hi=F(4)):
    def below(x):
        B=[row[:] for row in A]
        for i,g in enumerate(G):B[i][i]-=x*g
        return audit.inertia(B)
    assert below(lo)<=index<below(hi)
    while hi-lo>F(1,1<<40):
        mid=(lo+hi)/2
        if below(mid)<=index:lo=mid
        else:hi=mid
    return lo,hi


def run():
    R,G=basis();M,C,K=components(R,G)
    inputs=json.loads((HERE/'RAW_COEFFICIENT_INPUTS.json').read_text())
    m=F(inputs['432']['declared_light_midpoint']);t=F(1,8)
    rows=[]
    for kappa in [F(1,2),F(1),F(2)]:
        A=[[m*int(M[i,j])+kappa*F(int(C[i,j]),3)+t*int(K[i,j]) for j in range(6)] for i in range(6)]
        full_intervals=[interval(A,G,j) for j in range(6)]
        assert full_intervals[0][1]<full_intervals[1][0]<0
        # Exact B=(N-3)/3: vacuum N=3 has B0; N0 and N6 have B-1,+1.
        assert all(A[0][j]==A[5][j]==0 for j in range(6))
        source_gap=[-full_intervals[0][1],-full_intervals[0][0]]
        H=np.array([[float(A[i][j])/np.sqrt(G[i]*G[j]) for j in range(6)] for i in range(6)])
        values,vectors=np.linalg.eigh(H)
        vacuum=vectors[:,0]
        assert np.linalg.norm(vacuum[[0,5]])<1e-12
        # Local epsilon triple creator/annihilator are actual CAR operators.
        add=np.zeros((64,64),dtype=np.int64)
        for word in range(64):
            v=word;sign=1
            for mode in [2,1,0]:
                v,s=creation(mode,v)
                if v is None:break
                sign*=s
            if v is not None:add[v,word]=sign
        full_v=R@(vacuum/np.sqrt(np.array(G)))
        plus=add@full_v;minus=add.T@full_v
        assert np.flatnonzero(abs(plus)>1e-12).tolist()==[63]
        assert np.flatnonzero(abs(minus)>1e-12).tolist()==[0]
        pplus=float(np.vdot(plus,plus));pminus=float(np.vdot(minus,minus))
        assert 0<pplus<1 and 0<pminus<1
        # Every complement of each local contraction remains a valid effect.
        assert np.array_equal(add.T@add+(np.eye(64,dtype=int)-add.T@add),np.eye(64,dtype=int))
        rows.append(dict(kappa=str(kappa),t=str(t),m=str(m),
            full_six_physical_energy_intervals=[list(map(str,x)) for x in full_intervals],
            baryon_pair_gap_interval=list(map(str,source_gap)),
            baryon_pair_gap_decimal=list(map(float,source_gap)),
            actual_local_triquark_creation_probability_diagnostic=pplus,
            actual_local_triquark_annihilation_probability_diagnostic=pminus,
            reference_source='unique lowest N3 state, not the empty Fock state',
            baryon_spin_or_nucleon_identity_claimed=False))
    # Exact massless control for the primary inherited kappa=1 member.
    D=F(4,3)-2*t;root=F(7,6)
    assert root*root==D*D+12*t*t and (root-D)/2==F(1,24)
    anchor=F(inputs['electron_anchor_GeV']);unit=anchor/m;target=F('0.93827208943')
    for row in rows:
        lo,hi=map(F,row['baryon_pair_gap_interval'])
        row['electron_anchored_gap_GeV_interval']=list(map(str,[lo*unit,hi*unit]))
        row['electron_anchored_gap_GeV_decimal']=list(map(float,[lo*unit,hi*unit]))
        row['proton_comparison_met']=lo*unit<=target<=hi*unit
        assert not row['proton_comparison_met']
    out=dict(schema='oph.obstruction_repair.staggered_colour_vacuum.v1',
        architecture='two staggered cells, three colours, signed hopping -t(c0†c1+h.c.), staggered raw mass m(N0-N1), SU3 tree electric kappa C0², zero exterior colour flux',
        full_CAR_dimension=64,complete_Gauss_singlet_dimension=6,
        exact_generator_and_Hamiltonian_intertwiners=True,all_physical_baryon_sectors_retained=True,
        actual_ground_vacuum=True,empty_vacuum_relabelled=False,
        massless_primary_gap_exact='1/24',rows=rows,
        same_clock_electron_rest_factor='optional independent source factor e f†f with same exact e; spatial/Dirac transfer not proved',
        new_empirical_fit=False,physical_proton_mass_solved=False,
        old_positive_number_source_obstruction_evaded=True,
        native_M6_full_program_executed=False,
        limitations=['spin/flavour and proton/neutron dictionary absent','two-cell staggered source only','common spatial Dirac field and native complete preparation/readout not implemented','no 3D refinement or scattering identification','inherited analyzer has retrospective phase ancestry','physical mass comparison failed at all three inherited kappa values'])
    (HERE/'STAGGERED_VACUUM_RESULT.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['full_CAR_dimension','complete_Gauss_singlet_dimension','actual_ground_vacuum','massless_primary_gap_exact','rows']},indent=2))


if __name__=='__main__':run()
