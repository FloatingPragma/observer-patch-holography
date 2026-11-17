import json
import numpy as np
from scipy.linalg import expm,logm
from scipy.special import xlogy
from bound import low_rank_log_Z,finite_entropy_interval

def annihilators(n):
    out=[]
    for j in range(n):
        a=np.zeros((2**n,2**n),complex)
        for state in range(2**n):
            if (state>>j)&1:a[state^(1<<j),state]=(-1)**((state&((1<<j)-1)).bit_count())
        out.append(a)
    return out

def density(C,ops):
    e,U=np.linalg.eigh(C);h=(U*np.log((1-e)/e))@U.conj().T
    H=sum(h[i,j]*ops[i].conj().T@ops[j] for i in range(len(ops)) for j in range(len(ops)))
    rho=expm(-H);return rho/np.trace(rho),h

def entropy(rho):
    e=np.clip(np.linalg.eigvalsh(rho),0,1);return float(np.sum(-xlogy(e,e)))

rng=np.random.default_rng(7103);ops=annihilators(4);rows=[]
for sample in range(8):
    U,_=np.linalg.qr(rng.normal(size=(4,4))+1j*rng.normal(size=(4,4)))
    V,_=np.linalg.qr(rng.normal(size=(4,4))+1j*rng.normal(size=(4,4)))
    C0=(U*np.array([.08,.22,.67,.91]))@U.conj().T
    C1=(V*np.array([.13,.31,.58,.86]))@V.conj().T
    e,W=np.linalg.eigh(C1-C0);X=W*np.sqrt(abs(e));signs=np.sign(e)
    result=low_rank_log_Z(C0,X,signs)
    rho0,h=density(C0,ops);rho1,_=density(C1,ops)
    Z=float(np.trace(rho1@rho1@np.linalg.inv(rho0)).real)
    J=float(np.trace(h@(C1-C0)).real)
    tests=[]
    for epsilon in [.001,.03,.2]:
        interval=finite_entropy_interval(J,epsilon,result['log_Z'])
        actual=entropy((1-epsilon)*rho0+epsilon*rho1)-entropy(rho0)
        interval.update(actual_non_gaussian_entropy_change=actual,inside_bound=interval['entropy_lower']-1e-12<=actual<=interval['entropy_upper']+1e-12)
        tests.append(interval)
    rows.append(dict(sample=sample,full_Fock_Z=Z,low_rank_log_Z=result['log_Z'],log_Z_error=abs(np.log(Z)-result['log_Z']),phase=result['determinant_phase'],mixtures=tests))
out=dict(status='EXPLICIT_NONCOMMUTING_FOCK_MIXTURE_CONTROL',regional_modes=4,Fock_dimension=16,rows=rows,max_log_Z_error=max(r['log_Z_error'] for r in rows),all_actual_entropy_brackets_pass=all(t['inside_bound'] for r in rows for t in r['mixtures']))
with open('FOCK_CONTROL.json','w') as f:json.dump(out,f,indent=2)
print(json.dumps({k:v for k,v in out.items() if k!='rows'},indent=2),flush=True)
assert out['max_log_Z_error']<1e-10 and out['all_actual_entropy_brackets_pass']
