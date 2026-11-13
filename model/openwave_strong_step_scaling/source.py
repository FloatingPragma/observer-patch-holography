"""Complete native Gauss chains and exact energy-dependent coarse actions."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import factorial
import hashlib
import json
from pathlib import Path
import numpy as np
import scipy.sparse as sp
from scipy.linalg import eigh
from model.openwave_nuclear_distance.spatial import carrier
from model.openwave_common_hyperfine_spectrum import source as native

from model.openwave_string_positive_cone.source import sqrt_interval


@lru_cache(None)
def shared_length_attachment():
    # Typed, authenticated theorem output: the parent independently certifies
    # this exact scale. This consumer reexecutes its own full native chains.
    root=Path(__file__).resolve().parents[2]
    parent=json.loads((root/'model/openwave_common_hyperfine_spectrum/results.json').read_text())
    for path,digest in parent['source_sha256'].items():
        assert hashlib.sha256((root/path).read_bytes()).hexdigest()==digest,path
    assert parent['parameters_rational']=={k:str(v)for k,v in native.parameters().items()}
    return parent['common_length_attachment']


@lru_cache(None)
def chain(length):
    d, p = carrier(length), native.parameters()
    kappa = F(shared_length_attachment()['kappa'])
    n, counts, idx, E, angular = d['n'], d['counts'], d['index'], d['E'], d['angular']
    N = len(n)
    K = sp.csr_matrix((N,N)); S=K.copy(); I=K.copy(); edgeI=K.copy(); gradient=K.copy()
    for x in range(length):
        S += 4*angular((x,), True)
        I += 4*angular((x,), False)
    for x in range(length-1):
        for spin, f in product(range(2), repeat=2):
            a,b = idx[x,spin,f],idx[x+1,spin,f]
            K += E(a,a)+E(b,b)-E(a,b)-E(b,a)
        edgeI += 4*angular((x,x+1),False)
        gradient += 8*(angular((x,),True)+angular((x+1,),True))-4*angular((x,x+1),True)
    prefix = np.cumsum(n,axis=1)[:,:-1]
    localC3=2*np.sum(n*(3-n),axis=1); electric3=2*np.sum(prefix*(3-prefix),axis=1)
    spans=np.array([max(d['orbitals'][j][0] for j in lab)-min(d['orbitals'][j][0] for j in lab) for lab in d['labels']])
    assert np.array_equal(electric3,4*spans)
    nu=counts[:,[j for j,o in enumerate(d['orbitals']) if o[2]==0]].sum(axis=1)
    ns=counts[:,[j for j,o in enumerate(d['orbitals']) if o[1]==0]].sum(axis=1)
    parts={'N':3*sp.eye(N,dtype=np.int64),'D':sp.diags(3-2*nu,dtype=np.int64),
           'B9':sp.diags(np.sum(n*n,axis=1),dtype=np.int64),'C3':sp.diags(localC3,dtype=np.int64),
           'Electric3':sp.diags(electric3,dtype=np.int64),'K':K,'S4':S,'I4':I,'EdgeI4':edgeI,'Gradient4':gradient}
    weights={'N':p['mbar'],'D':p['delta']/2,'B9':p['A']/9,'C3':(p['U']-kappa/2)/3,
             'Electric3':kappa/3,'K':p['t'],'S4':p['Js']/4,'I4':p['Jw']/4,
             'EdgeI4':p['eta']/4,'Gradient4':p['spin_gradient']/4}
    H=sum(float(weights[k])*a for k,a in parts.items()).tocsr()
    row,col=H.nonzero()
    assert np.all(nu[row]==nu[col]) and np.all(ns[row]==ns[col])
    gram=np.array([6//np.prod([factorial(int(v)) for v in a]) for a in counts],dtype=np.int64)
    Jp=sp.csr_matrix((N,N)); Ip=Jp.copy()
    for x,spin in product(range(length),range(2)):
        term=E(idx[x,spin,0],idx[x,spin,1]);Ip+=term;Jp+=(1-2*spin)*term
    if length==2:
        parent=native.sector('baryon')
        overlap=parent['B'].T@d['B']
        assert np.max(abs(parent['h']@overlap-overlap@H))<4e-10
    return {**d,'H':H,'parts':parts,'weights':weights,'spans':spans,'nu':nu,'ns':ns,
            'gram':gram,'Ip':Ip,'Jp':Jp,'kappa':kappa,'parameters':p}


def block(d, nu):
    ids=np.flatnonzero((d['nu']==nu)&(d['ns']==2))
    h=d['H'][ids][:,ids].toarray()
    p=np.flatnonzero(d['spans'][ids]<=1);q=np.flatnonzero(d['spans'][ids]>1)
    gram=d['gram'][ids]
    # Each exact part is integer in the inherited unnormalized symmetric CAR basis.
    exact=np.full(h.shape,F(0),dtype=object)
    for name,w in d['weights'].items():
        a=d['parts'][name].tocsr()[ids][:,ids].toarray()*np.sqrt(np.outer(gram,gram))
        rounded=np.rint(a).astype(np.int64)
        assert np.max(abs(a-rounded))<2e-10
        for i,j in zip(*rounded.nonzero()):exact[i,j]+=w*int(rounded[i,j])
    assert np.max(abs(np.asarray(exact,float)/np.sqrt(np.outer(gram,gram))-h))<4e-10
    return {'ids':ids,'H':h,'P':p,'Q':q,'gram':gram,'exact_numerator':exact}


def complement_bounds(b, energy_upper):
    q,p,g,R=b['Q'],b['P'],b['gram'],b['exact_numerator']
    if not len(q):return {'floor':None,'gap':None,'coupling_squared':F(0),'tau':F(0)}
    inverse_root={int(a):sqrt_interval(F(1,int(a)))[1] for a in np.unique(np.outer(g,g))}
    lower=min(R[i,i]/int(g[i])-sum(abs(R[i,j])*inverse_root[int(g[i]*g[j])] for j in q if j!=i) for i in q)
    gap=lower-F(energy_upper)
    assert gap>0
    v2=sum(R[i,j]**2/F(int(g[i]*g[j])) for i in p for j in q)
    tau=sqrt_interval(v2)[1]/gap
    return {'floor':lower,'gap':gap,'coupling_squared':v2,'tau':tau}


def effective(b, energy):
    h,p,q=b['H'],b['P'],b['Q'];A=h[np.ix_(p,p)]
    R=np.zeros((len(h),len(p)));R[p]=np.eye(len(p))
    if len(q):
        D=h[np.ix_(q,q)];V=h[np.ix_(p,q)]
        R[q]=np.linalg.solve(energy*np.eye(len(q))-D,V.T)
        Sigma=V@R[q]
    else:Sigma=np.zeros_like(A)
    metric=R.T@R
    e,U=eigh(metric);Z=R@(U@np.diag(1/np.sqrt(e))@U.T)
    assert np.max(abs(Z.T@Z-np.eye(len(p))))<3e-12
    return {'R':R,'Z':Z,'metric':metric,'Sigma':Sigma,'H_effective':A+Sigma}
