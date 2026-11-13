"""Five-band compact-cluster symbol extracted from the full native chain."""
from functools import lru_cache
import numpy as np
from scipy.linalg import eigh
from .source import chain


def canonical_range(projector,rank):
    columns=[]
    for j in range(len(projector)):
        v=projector[:,j].copy()
        for q in columns:v-=q*(q@v)
        norm=np.linalg.norm(v)
        if norm>1e-9:
            v/=norm
            if v[np.flatnonzero(abs(v)>1e-9)[0]]<0:v=-v
            columns.append(v)
            if len(columns)==rank:break
    B=np.array(columns).T
    assert np.max(abs(B.T@B-np.eye(rank)))<1e-12
    return B


@lru_cache(None)
def compact_frame(L=4):
    d=chain(L);N=len(d['labels']);S=d['angular'](tuple(range(L)),True);I=d['angular'](tuple(range(L)),False)
    columns=[];groups=[]
    occupations=[]
    for x in range(L):
        a=np.zeros(L,dtype=int);a[x]=3;occupations.append(('site',x,a))
    for x in range(L-1):
        for nl in (2,1):
            a=np.zeros(L,dtype=int);a[x]=nl;a[x+1]=3-nl;occupations.append(('edge',x,a))
    for kind,x,occupation in occupations:
        ids=np.flatnonzero(np.all(d['n']==occupation,axis=1)&(d['nu']==1)&(d['ns']==2))
        s=S[ids][:,ids].toarray();i=I[ids][:,ids].toarray();eye=np.eye(len(ids))
        filt=(s-.75*eye)@(s-.75*eye)+(i-.75*eye)@(i-.75*eye)
        e,V=eigh(filt);V=V[:,e<1e-10];B=canonical_range(V@V.T,V.shape[1]);rank=len(B.T)
        assert rank==(1 if kind=='site' else 2)
        indices=list(range(len(columns),len(columns)+rank))
        for j in range(rank):
            v=np.zeros(N);v[ids]=B[:,j];columns.append(v)
        groups.append((kind,x,tuple(occupation),indices))
    P=np.array(columns).T
    assert np.max(abs(P.T@P-np.eye(len(columns))))<1e-12
    H=P.T@d['H']@P
    J=P.T@d['Ip'].T@d['Jp']@P
    return {'d':d,'P':P,'H':H,'J':J,'groups':groups}


@lru_cache(None)
def symbol_data():
    z=compact_frame(4);H=z['H'];J=z['J'];groups=z['groups']
    sites={x:indices for kind,x,n,indices in groups if kind=='site'}
    edges={x:sum([inds for kind,xx,n,inds in groups if kind=='edge' and xx==x],[]) for x in range(3)}
    a=H[sites[0][0],sites[0][0]];b=H[np.ix_(edges[0],edges[0])]
    vl=H[sites[0]][:,edges[0]].ravel();vr=H[sites[1]][:,edges[0]].ravel()
    jr=J[np.ix_(edges[0],edges[0])];ja=J[sites[0][0],sites[0][0]]
    # All finite boundary and interior matrix elements obey the same symbol.
    rebuilt=np.zeros_like(H);jrebuilt=np.zeros_like(J)
    for x in sites:
        rebuilt[np.ix_(sites[x],sites[x])]=a;jrebuilt[np.ix_(sites[x],sites[x])]=ja
    for x,ids in edges.items():
        rebuilt[np.ix_(ids,ids)]=b;jrebuilt[np.ix_(ids,ids)]=jr
        for site,v in ((x,vl),(x+1,vr)):
            rebuilt[np.ix_(sites[site],ids)]=v
            rebuilt[np.ix_(ids,sites[site])]=v[:,None]
    assert np.max(abs(rebuilt-H))<2e-11,np.max(abs(rebuilt-H))
    assert np.max(abs(jrebuilt-J))<2e-12
    residual=z['d']['H']@z['P']-z['P']@H
    square=residual.T@residual
    T0=square[np.ix_(edges[1],edges[1])]
    T1=square[np.ix_(edges[1],edges[2])]
    assert np.max(abs(square[np.ix_(edges[0],edges[2])]))<1e-12
    assert np.max(abs(square[np.ix_(edges[0],edges[1])]-T1))<1e-12
    assert np.max(abs(square[:4]))<1e-18
    return {'a':a,'b':b,'vl':vl,'vr':vr,'ja':ja,'j':jr,'T0':T0,'T1':T1,
            'finite_intertwiner_error':float(np.max(abs(rebuilt-H)))}


def symbol(k):
    d=symbol_data();H=np.zeros((5,5),complex);J=H.copy()
    H[0,0]=d['a'];H[1:,1:]=d['b'];H[0,1:]=d['vl']+np.exp(-1j*k)*d['vr'];H[1:,0]=H[0,1:].conj()
    J[0,0]=d['ja'];J[1:,1:]=d['j']
    return H,J


def coupling_square(k):
    d=symbol_data();T=d['T0']+np.exp(1j*k)*d['T1']+np.exp(-1j*k)*d['T1'].T
    assert np.linalg.eigvalsh(T)[0]>-1e-12
    out=np.zeros((5,5),complex);out[1:,1:]=T
    return out


if __name__=='__main__':
    import json
    for k in (0,np.pi/4,np.pi/2,np.pi):
        H,J=symbol(k);e,V=eigh(H);v=V[:,0]
        T=coupling_square(k)
        print(json.dumps({'k':k,'all_five_energies_GeV':(e*.00051099895069).tolist(),'gA':float(np.vdot(v,J@v).real),'gap_source':float(e[1]-e[0]),'coupling_squared_norm':float(np.linalg.eigvalsh(T)[-1]),'ground_coupling_squared':float(np.vdot(v,T@v).real)}),flush=True)
