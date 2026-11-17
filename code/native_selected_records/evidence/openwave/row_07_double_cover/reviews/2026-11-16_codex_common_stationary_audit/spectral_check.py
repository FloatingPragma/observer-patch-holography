"""Independent four-mode spectral counter and analytic native-transfer check.
No producer imports; only NumPy and standard library. Verification, not a run gate.
"""
import json, math
import numpy as np

I=np.eye(3); edges=((0,1),(1,2),(0,2))
def lap(edge):
    v=np.zeros(3);v[edge[0]]=1;v[edge[1]]=-1
    return np.outer(v,v)
L=2*lap((0,1))+2*lap((0,2))+lap((1,2));P=I-L/10
f=np.array([0.,1.,-1.])/math.sqrt(2);g=np.array([-2.,1.,1.])/math.sqrt(6)
E=np.column_stack((f,g));perm=I[[1,2,0]]
D=np.c_[np.eye(2),np.zeros((2,2))]
B=np.zeros((4,4),complex)
for m in (0,1):
    B[m,2*m]=B[m,2*m+1]=1/math.sqrt(2)
    B[2+m,2*m]=-1j/math.sqrt(2);B[2+m,2*m+1]=1j/math.sqrt(2)

def counter(k):
    eps=2.**(-k);Q=256*2**k
    phases=np.array([t for lam in (4,6) for t in (2*math.asin(eps*math.sqrt(lam/20)),-2*math.asin(eps*math.sqrt(lam/20)))])
    delta=phases[None,:]-2*math.pi*np.arange(Q)[:,None]/Q
    values=np.exp(.5j*(Q-1)*delta)*np.sinc(Q*delta/(2*math.pi))/np.sinc(delta/(2*math.pi))
    F=np.einsum('ai,ri,bi->rab',B,values,B.conj())
    assert np.linalg.norm(np.einsum('rxa,rxb->ab',F.conj(),F)-np.eye(4))<1e-9
    plus=F[1:Q//2]
    return np.einsum('rab,rcd->acbd',plus,plus.conj())

def cp(C,rho):return np.einsum('acbd,bd->ac',C,rho)
def prep(C):
    v=(I-P)@np.array([0.,1.,0.])/(2*math.sqrt(3))
    for _ in range(8):v=P@v/math.sqrt(3)
    p=float(v@v);v/=math.sqrt(p)
    optical=E.T@np.outer(v,v)@E
    rho=D.T@optical@D;rho=cp(C,rho);q=float(np.trace(rho).real);p*=q;rho/=q
    position=D@rho@D.T;q=float(np.trace(position).real);p*=q;position/=q
    return position,math.log10(p)
def mean(slot,m):
    v,k=divmod(slot,m);w=(v+1)%3
    star=lambda center:I-sum(lap(e) for e in edges if center in e)/4
    return (1-k/m)*star(v)+(k/m)*star(w)
def normalize(rho):
    p=float(np.trace(rho).real);assert p>0
    return rho/p,math.log10(p)

rows=[]
for k in range(7):
    C=counter(k);position,preplog=prep(C)
    for m in (4,8):
        means=[mean(slot,m) for slot in range(3*m)]
        for reference in ('matched','stationary'):
            for law in ('native','swap'):
                for turns in (1,2):
                    for direction in (1,-1):
                        rho=np.kron(np.ones((2,2))/2,position);reach=preplog
                        frames=3*m*turns
                        reach+=frames*(2 if reference=='matched' else 1)*math.log10(1-2**-8)
                        for step in range(1,frames+1):
                            q=step%m;rot=(2*m+direction*step)%(3*m)
                            ref=(2*m+direction*min(q,m-q))%(3*m) if reference=='matched' else 2*m
                            G=perm@E
                            a,b=G.T@means[rot]@G,G.T@means[ref]@G
                            if law=='swap':a,b=2*a-np.eye(2),2*b-np.eye(2)
                            op=np.zeros((4,4));op[:2,:2]=a/math.sqrt(3);op[2:,2:]=b/math.sqrt(3)
                            for _ in range(4):rho,dlog=normalize(op@rho@op.T);reach+=dlog
                            # Full plane filter is scalar on this invariant optical input.
                            reach+=math.log10(1/(48 if law=='native' else 12))
                        blocks=rho.reshape(2,2,2,2);out=np.empty((2,4,2,4),complex)
                        for a in (0,1):
                            for b in (0,1):
                                pos=blocks[a,:,b,:]
                                out[a,:,b,:]=cp(C,D.T@pos@D)
                        out,dlog=normalize(out.reshape(8,8));reach+=dlog
                        blocks=out.reshape(2,4,2,4)
                        cross=float(np.trace(blocks[0,:,1,:]).real)
                        visibility=2*cross;sign=-1 if law=='native' and turns==1 else 1
                        path=(1+sign*visibility)/2
                        weights=[float(np.trace(blocks[a,:,a,:]).real) for a in (0,1)]
                        normalized=cross/math.sqrt(weights[0]*weights[1])
                        assert visibility<0 if sign==-1 else visibility>0
                        s5=1-(3/8)**8;s3=1-(1/4)**8
                        reach+=11*math.log10(s5)+frames*math.log10(s3)
                        rows.append(dict(k=k,M=m,reference=reference,law=law,turns=turns,direction=direction,
                            visibility=visibility,normalized_cross=normalized,path_probability=path,
                            log_interference=reach,log_clock=reach+math.log10(path)))
print(json.dumps({'method':'Four-mode analytic spectral counter, native Laplacian transfers; no producer imports',
                  'cases':rows,'count':len(rows),'passed':True,'external_G6':False},indent=2))
