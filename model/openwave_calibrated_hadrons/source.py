"""One positive native-current action and its complete exact colour sectors."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import json
from pathlib import Path
import numpy as np
import scipy.sparse as sp
from scipy.linalg import eigh
from scipy.optimize import brentq
import sympy as sy
from model.openwave_native_composites.spectrum import sector as old_sector,one_body
from model.openwave_baryon_ordering.source import balanced_baryon
from model.coherent_defect_wave.physical_attachment.chiral import representation

HERE=Path(__file__).resolve().parent
INPUTS=json.loads((HERE/'inputs.json').read_text())
DECLARED={'U':F(2048),'Js':F(1,4),'Jw':F(1,4),'t':F(1,8),'eta':F(1)}


def integral(matrix):
    a=matrix.tocsr();assert np.max(abs(a.data.imag),initial=0)<1e-12
    assert np.max(abs(a.data.real-np.rint(a.data.real)),initial=0)<2e-11
    out=sp.csr_matrix((np.rint(a.data.real).astype(np.int64),a.indices,a.indptr),shape=a.shape)
    out.eliminate_zeros();return out


def native_colour48(modes,basis):
    """Actual native integer root currents and orthogonal Cartan metric."""
    subsets,_,g,_=representation()
    qids=[[subsets.index((c,w))for c in range(3)]for w in(3,4)]
    aids=[[subsets.index(tuple(c for c in range(3)if c!=missing))for missing in range(3)],
          [subsets.index(tuple(c for c in range(5)if c!=missing))for missing in range(3)]]
    index={x:i for i,x in enumerate(modes)};n=len(basis);out=sp.csr_matrix((n,n),dtype=np.int64)
    generators=[2*a for a in g[:7]]+[2*(g[6]+2*g[7])]
    for cell in range(2):
        for axis,gen in enumerate(generators):
            imaginary=any(sy.im(x)!=0 for x in gen)
            raw=gen/sy.I if imaginary else gen
            A=np.zeros((len(modes),len(modes)),dtype=np.int64)
            for species,_,spin,f,c in modes:
                ids=(qids if species==0 else aids)[f]
                for cc in range(3):
                    value=raw[ids[c],ids[cc]];assert value.is_Integer
                    A[index[species,cell,spin,f,c],index[species,cell,spin,f,cc]]=int(value)
            J=integral(one_body(A,basis));coefficient=(4 if axis==7 else 12)*(-1 if imaginary else 1)
            out+=coefficient*(J@J)
    return out


def angular4(modes,basis,cells,spin):
    index={x:i for i,x in enumerate(modes)};A=np.zeros((len(modes),len(modes)),dtype=np.int64)
    for species,cell,ss,f,c in modes:
        if cell not in cells or (not spin and species!=0):continue
        if spin:
            if ss==0:A[index[species,cell,0,f,c],index[species,cell,1,f,c]]=1
        elif f==0:A[index[0,cell,ss,0,c],index[0,cell,ss,1,c]]=1
    R=integral(one_body(A,basis))
    z=np.array([sum((1 if (ss if spin else f)==0 else -1)
         for i,(species,cell,ss,f,c)in enumerate(modes)
         if cell in cells and(spin or species==0)and(mask>>i)&1)for mask in basis],dtype=np.int64)
    return (sp.diags(z*z+2*z,dtype=np.int64)+4*(R.T@R)).tocsr()


@lru_cache(None)
def components(kind):
    old=old_sector(kind);modes=old['modes'];basis=old['basis'];n=len(basis);number=3 if kind=='baryon'else 2
    C=native_colour48(modes,basis);assert np.max(abs((C-48*old['C']).data),initial=0)<3e-11
    K=integral(8*old['kinetic']);S=angular4(modes,basis,(0,),True)+angular4(modes,basis,(1,),True)
    I=angular4(modes,basis,(0,),False)+angular4(modes,basis,(1,),False);Ig=angular4(modes,basis,(0,1),False)
    B=[];difference=[]
    for mask in basis:
        local=[sum((1 if species==0 else -1)for i,(species,x,ss,f,c)in enumerate(modes)if x==cell and(mask>>i)&1)for cell in range(2)]
        B.append(sum(x*x for x in local));difference.append(sum((1 if f else -1)for i,(_,x,ss,f,c)in enumerate(modes)if(mask>>i)&1))
    parts={'N':number*sp.eye(n,dtype=np.int64,format='csr'),'D':sp.diags(difference,dtype=np.int64,format='csr'),
           'K':K,'C48':C,'S4':S,'I4':I,'B9':sp.diags(B,dtype=np.int64,format='csr'),'Ig4':Ig}
    if kind=='baryon':
        ref=balanced_baryon();R=ref['Binteger'];gram=ref['Gram']
    else:gram=np.full(len(old['labels']),3,dtype=np.int64);R=integral(old['B']@sp.diags(np.sqrt(gram)))
    assert (R.T@R-sp.diags(gram,dtype=float)).nnz==0
    reduced={name:(R.T@a@R).toarray()for name,a in parts.items()}
    denominator=np.sqrt(np.outer(gram,gram))
    hparts={name:a/denominator for name,a in reduced.items()}
    for name,a in parts.items():assert np.max(abs(a@old['B']-old['B']@hparts[name]))<2e-11
    # Fixed-charge, fixed-spin finite subspaces, constructed from their full operators.
    filt=(old['q']-np.eye(len(gram)))@(old['q']-np.eye(len(gram)))
    target=F(3,4)if kind=='baryon'else F(0)
    filt+=(old['spin']-float(target)*np.eye(len(gram)))@(old['spin']-float(target)*np.eye(len(gram)))
    e,v=eigh(filt);P=v[:,e<1e-8]
    return {**old,'parts':parts,'integer_reduced_parts':reduced,'hparts':hparts,'Binteger':R,'Gram':gram,'P_calibration':P}


def weights(p):
    return {'N':p['mbar'],'D':p['delta']/2,'K':p['t'],'C48':p['U']/48,
            'S4':p['Js']/4,'I4':p['Jw']/4,'B9':p['A']/9,'Ig4':p['eta']/4}


def matrix(kind,p):
    d=components(kind)
    return sum(float(v)*d['hparts'][k]for k,v in weights(p).items())


@lru_cache(None)
def parameters():
    unit=F(INPUTS['electron_unit_GeV']);targets={k:F(v)/unit for k,v in INPUTS['calibration_targets_GeV'].items()}
    delta=targets['neutron']-targets['proton'];targetbar=(targets['proton']+targets['neutron'])/2
    p={**DECLARED,'mbar':F(0),'delta':F(0),'A':F(0)}
    def energy(kind,A):
        H=matrix(kind,{**p,'A':A});P=components(kind)['P_calibration'];return float(eigh(P.conj().T@H@P,eigvals_only=True)[0])
    def objective(A):return energy('baryon',A)-1.5*energy('meson',A)-float(targetbar-F(3,2)*targets['charged_pion'])
    assert objective(0)<0<objective(4096)
    fit=brentq(objective,0.,4096.,xtol=3e-12);A=F(round(fit*2**36),2**36)
    mbar=F(round((float(targets['charged_pion'])-energy('meson',A))/2*2**36),2**36)
    assert A>0 and mbar>delta/2>0
    return {**DECLARED,'mbar':mbar,'delta':delta,'A':A}


@lru_cache(None)
def sector(kind):
    p=parameters();d=components(kind);H=sum(float(v)*d['parts'][k]for k,v in weights(p).items());h=matrix(kind,p)
    assert np.linalg.norm(H@d['B']-d['B']@h)<2e-8
    return {**d,'H':H,'h':h,'parameters':p}


def baryon_states():
    """Root beta API: full colour sector; no truncation to nucleons or one lepton pair."""
    return sector('baryon')


def exact_matrix(kind):
    p=parameters();d=components(kind);g=d['Gram'];parts=d['integer_reduced_parts'];w=weights(p);n=len(g)
    A=sy.zeros(n)
    for name,weight in w.items():A+=sy.Rational(weight.numerator,weight.denominator)*sy.Matrix(parts[name].tolist())
    return sy.Matrix(n,n,lambda i,j:A[i,j]/sy.sqrt(int(g[i])*int(g[j])))
