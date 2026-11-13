"""Two native colour families on the declared two-cell strong source.

The exact colour-singlet reductions retain all 816 three-quark and 256
quark-antiquark states. Light-only blocks equal the seven-input action.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import json
from pathlib import Path
import numpy as np
import scipy.sparse as sp
from scipy.linalg import eigh
from model.openwave_nuclear_distance.spatial import carrier
from model.openwave_common_hyperfine_spectrum import source as parent
from model.openwave_native_composites.spectrum import native, one_body

ORBITALS=tuple((x,family,s,weak) for x,family,s,weak in product(range(2),repeat=4))
INDEX={v:i for i,v in enumerate(ORBITALS)}


def lifted_generator(d, matrix, species=None):
    out=sp.csr_matrix((len(d['labels']),len(d['labels'])),dtype=complex)
    for i,j in zip(*matrix.nonzero()):
        if d['kind']=='baryon':out+=matrix[i,j]*d['E'](i,j)
        elif species is None:
            out+=matrix[i,j]*(d['E'](0,i,j)+d['E'](1,i,j))
        else:out+=matrix[i,j]*d['E'](species,i,j)
    return out


def angular(d,cells,spin,species=None):
    R=np.zeros((16,16),dtype=int);z=np.zeros(16,dtype=int)
    for i,(x,family,s,w) in enumerate(ORBITALS):
        if x not in cells or (not spin and family!=0):continue
        axis=s if spin else w
        z[i]=1-2*axis
        if axis==0:
            j=INDEX[x,family,1,w] if spin else INDEX[x,family,s,1]
            R[i,j]=1
    # The inherited strong weak Casimir acts on quarks only, including mesons.
    chosen=species if spin or species is not None else 0
    Rp=lifted_generator(d,R,chosen)
    Z=lifted_generator(d,np.diag(z),chosen)
    return ((Z@Z+2*Z)/4+Rp.conj().T@Rp).tocsr()


@lru_cache(None)
def reduced_carrier(kind):
    if kind=='baryon':
        base=carrier(4)
        assert all(ORBITALS[i]==(x//2,x%2,s,w)
                   for i,(x,s,w) in enumerate(base['orbitals']))
        return {**base,'kind':kind,'counts_q':base['counts'],
                'counts_a':np.zeros_like(base['counts'])}
    assert kind=='meson'
    labels=list(product(range(16),repeat=2));size=len(labels)
    countsq=np.zeros((size,16),dtype=int);countsa=countsq.copy()
    for k,(q,a) in enumerate(labels):countsq[k,q]=1;countsa[k,a]=1
    def E(species,i,j):
        unit=sp.csr_matrix(([1.],([i],[j])),shape=(16,16))
        return sp.kron(unit,sp.eye(16)) if species==0 else sp.kron(sp.eye(16),unit)
    # Literal native q-antiquark colour contraction, including every flavour.
    modes=tuple((species,*o,c) for species,o,c in product(range(2),ORBITALS,range(3)))
    mode_index={v:i for i,v in enumerate(modes)}
    basis=tuple((1<<i)|(1<<j) for i in range(48) for j in range(48,96))
    lookup={v:i for i,v in enumerate(basis)}
    joins=native()[2];rows=[];cols=[];values=[]
    for k,(q,a) in enumerate(labels):
        for c,cc in product(range(3),repeat=2):
            value=joins[ORBITALS[a][3]][c,cc]/np.sqrt(3)
            if value:
                mask=(1<<mode_index[(0,*ORBITALS[q],c)])|(1<<mode_index[(1,*ORBITALS[a],cc)])
                rows.append(lookup[mask]);cols.append(k);values.append(value)
    B=sp.csr_matrix((values,(rows,cols)),shape=(2304,256))
    assert np.max(abs((B.conj().T@B-sp.eye(256)).data),initial=0)<1e-14
    return {'kind':kind,'labels':labels,'counts_q':countsq,'counts_a':countsa,
            'E':E,'B':B,'full':basis,'modes':modes}


@lru_cache(None)
def source(kind):
    d=reduced_carrier(kind);p=parent.parameters();cq,ca=d['counts_q'],d['counts_a']
    counts=cq+ca;size=len(counts)
    localq=np.array([cq[:,[i for i,o in enumerate(ORBITALS) if o[0]==x]].sum(axis=1) for x in range(2)]).T
    locala=np.array([ca[:,[i for i,o in enumerate(ORBITALS) if o[0]==x]].sum(axis=1) for x in range(2)]).T
    number=3 if kind=='baryon' else 2
    difference=counts@np.array([2*o[3]-1 for o in ORBITALS])
    strange=counts[:,[i for i,o in enumerate(ORBITALS) if o[1:2]==(1,) and o[3]==1]].sum(axis=1)
    charm=counts[:,[i for i,o in enumerate(ORBITALS) if o[1:2]==(1,) and o[3]==0]].sum(axis=1)
    S=sum((angular(d,(x,),True) for x in range(2)),sp.csr_matrix((size,size)))
    I=sum((angular(d,(x,),False) for x in range(2)),sp.csr_matrix((size,size)))
    Sg=angular(d,(0,1),True);Ig=angular(d,(0,1),False)
    gradient=2*S-Sg
    hop=np.eye(16)
    for i,(x,f,s,w) in enumerate(ORBITALS):hop[i,INDEX[1-x,f,s,w]]=-1
    K=lifted_generator(d,hop)
    if kind=='baryon':
        colour=F(2,3)*np.sum(localq*(3-localq),axis=1)
        pair=sp.csr_matrix((size,size))
    else:
        colour=np.array([F(0) if ORBITALS[q][0]==ORBITALS[a][0] else F(8,3) for q,a in d['labels']])
        pair=Sg
    diag=number*float(p['mbar'])+float(p['delta']/2)*difference+float(p['A']/9)*np.sum((localq-locala)**2,axis=1)+float(p['U'])*np.asarray(colour,float)
    H=sp.diags(diag)+float(p['t'])*K+float(p['Js'])*S+float(p['Jw'])*I+float(p['eta'])*Ig+float(p['spin_gradient'])*gradient+float(p['pair_triplet'])*pair
    # Declared C-symmetric extension. In a light one-pair block both separate
    # isospin Casimirs are3/4, so the subtraction exactly preserves that full
    # block. In pure-quark blocks both added terms vanish. On the complete
    # carrier its lower bound is -3(Jw+eta)N/8.
    nsq=cq[:,[i for i,o in enumerate(ORBITALS) if o[1]==1 and o[3]==1]].sum(axis=1)
    if kind=='meson':
        Ia=sum((angular(d,(x,),False,1) for x in range(2)),sp.csr_matrix((size,size)))
        Iag=angular(d,(0,1),False,1)
        H+=float(p['Jw'])*Ia+float(p['eta'])*Iag
        H-=float(F(3,4)*(p['Jw']+p['eta']))*sp.eye(size)
    assert np.max(abs((H-H.conj().T).data),initial=0)<1e-12
    # Exact calibrated-sector transfer, in the inherited native colour basis.
    light_ids=np.flatnonzero(strange+charm==0)
    old=parent.sector(kind)
    old_labels=old['labels']
    lookup={tuple(lab):i for i,lab in enumerate(d['labels'])}
    old_to_new=[]
    for lab in old_labels:
        target=tuple(INDEX[x,0,s,w] for x,s,w in lab)
        old_to_new.append(lookup[target])
    assert set(light_ids)==set(old_to_new)
    error=float(np.max(abs(H[old_to_new][:,old_to_new].toarray()-old['h'])))
    assert error<5e-10,error
    nq=cq.sum(axis=1);na=ca.sum(axis=1)
    nsa=ca[:,[i for i,o in enumerate(ORBITALS) if o[1]==1 and o[3]==1]].sum(axis=1)
    meson_strange_pair=nsq*na+nsa*nq
    return {**d,'H0':H.tocsr(),'S':Sg,'I':Ig,'Nstrange':strange,'Ncharm':charm,
            'strange_pair':meson_strange_pair,'light_transfer_error':error,
            'sz':counts@np.array([1-2*o[2] for o in ORBITALS]),
            'flavour_q':np.array([cq[:,[i for i,o in enumerate(ORBITALS) if 2*o[1]+o[3]==f]].sum(axis=1) for f in range(4)]).T,
            'flavour_a':np.array([ca[:,[i for i,o in enumerate(ORBITALS) if 2*o[1]+o[3]==f]].sum(axis=1) for f in range(4)]).T}


def select(kind,quarks,antiquarks,spin,isospin,dm=0.,pair_shift=0.):
    d=source(kind)
    ids=np.flatnonzero(np.all(d['flavour_q']==quarks,axis=1)&np.all(d['flavour_a']==antiquarks,axis=1)&(d['sz']==round(2*spin)))
    s=d['S'][ids][:,ids].toarray();i=d['I'][ids][:,ids].toarray();eye=np.eye(len(ids))
    selector=(s-spin*(spin+1)*eye)@(s-spin*(spin+1)*eye)+(i-isospin*(isospin+1)*eye)@(i-isospin*(isospin+1)*eye)
    z,P=eigh(selector);P=P[:,z<1e-9]
    assert P.shape[1],(kind,quarks,antiquarks,spin,isospin)
    H=d['H0']+sp.diags(dm*d['Nstrange']+pair_shift*d['strange_pair']+10000*d['Ncharm'])
    h=H[ids][:,ids].toarray();ev,U=eigh(P.conj().T@h@P)
    vector=np.zeros(len(d['labels']),complex);vector[ids]=P@U[:,0]
    seed=int(np.flatnonzero(abs(vector)>np.max(abs(vector))/2)[0])
    vector/=vector[seed]/abs(vector[seed])
    residual=float(np.linalg.norm(H@vector-ev[0]*vector));assert residual<1e-9
    return {'energy':float(ev[0]),'all_levels':ev.tolist(),'dimension':P.shape[1],
            'residual':residual,'state':vector}


if __name__=='__main__':
    unit=float(parent.INPUTS['electron_unit_GeV'])
    L=select('baryon',(1,1,0,1),(0,0,0,0),.5,0)
    K=select('meson',(1,0,0,0),(0,0,0,1),0,.5)
    dm=1.115683/unit-L['energy'];pair=.493677/unit-K['energy']-dm
    print(json.dumps({'Lambda_base_GeV':L['energy']*unit,'K_base_GeV':K['energy']*unit,'dm':dm,'strange_pair':pair}),flush=True)
    assert dm>0 and pair>0
    for name,kind,q,a,s,i in [('Lambda','baryon',(1,1,0,1),(0,0,0,0),.5,0),('Sigma0','baryon',(1,1,0,1),(0,0,0,0),.5,1),('Xi0','baryon',(1,0,0,2),(0,0,0,0),.5,.5),('K+','meson',(1,0,0,0),(0,0,0,1),0,.5),('K0','meson',(0,1,0,0),(0,0,0,1),0,.5),('K-','meson',(0,0,0,1),(1,0,0,0),0,0)]:
        z=select(kind,q,a,s,i,dm,pair)
        print(json.dumps({'name':name,'GeV':z['energy']*unit,'dimension':z['dimension'],'residual':z['residual']}),flush=True)
