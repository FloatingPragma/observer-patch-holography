#!/usr/bin/env python3
"""Post-run independent artifact audit. No native source or compiler calls."""
from fractions import Fraction as F
import gzip
import hashlib
import json
from pathlib import Path
import time
import numpy as np

HERE=Path(__file__).resolve().parent

def distribution(rho,effects):
    result=[]
    for label in range(256):
        effect=np.array([[1.]],complex)
        for bit in range(7,-1,-1):effect=np.kron(effect,effects[(label>>bit)&1])
        result.append(np.sum(rho*effect.T).real)
    return np.array(result)

def main():
    started=time.perf_counter()
    receipt=json.loads((HERE/'RECEIPT.json').read_text())
    arrays=np.load(HERE/'evidence/raw_source_arrays.npz')
    laws=np.load(HERE/'evidence/complete_outcome_laws.npz')
    # Recompute every law from stored source matrices and pointer effects.
    vectors={name:arrays['U_'+name][:,0] for name in ('f','g')}
    pure={name:np.outer(v,v.conj()) for name,v in vectors.items()}
    mixed={name:np.diag(np.diag(rho)) for name,rho in pure.items()}
    plant=vectors['f'].copy();plant[4]*=-1
    plant={'f':np.outer(plant,plant.conj()),'g':pure['g']}
    max_error=0.
    for law in ('C','E'):
        effects=arrays['effects_'+law]
        for arm,states,energy in (
            ('energy_coherent',pure,True),('energy_dephased',mixed,True),
            ('energy_sign_plant',plant,True),('position_coherent',pure,False),
            ('position_dephased',mixed,False)):
            p={}
            for name,rho in states.items():
                if energy:
                    U=arrays['inverse_'+name];rho=U@rho@U.conj().T
                p[name]=distribution(rho,effects)
            recomputed=np.outer(p['f'],p['g'])
            max_error=max(max_error,float(np.max(abs(recomputed-laws[law+'_'+arm]))))
    assert max_error<5e-12
    # Independently retain every native auxiliary tape terminal outcome.
    with gzip.open(HERE/'evidence/auxiliary_reset_tapes.json.gz','rt') as f:
        tapes=json.load(f)
    checks=0
    for row in tapes:
        state=list(row['initial']);execution=row['execution']
        for (a,b),coin,label in execution['records']:
            d=state[a]-state[b];j=(d-(2*coin-1)*(abs(d)%2))//2
            assert label==(j if row['law']=='C' else d)
            state[a]-=j;state[b]+=j;checks+=1
        assert state==execution['output']
        if execution['outcome']['success']:
            assert execution['outcome']['initial']==row['initial']
            assert state==[1,0,1,0,0,0,0,0,0,0,0,0]
    # Exact logarithm series, including its positive tail enclosure.
    for lam in (6,12):
        x=F(lam,1260);lo=1260*sum((x**k/k for k in range(1,33)),F())
        hi=lo+1260*x**33/(33*(1-x))
        assert [str(lo),str(hi)]==receipt['energy']['exact_log_intervals'][str(lam)]
    assert F(209,210)**2/F(104,105)==F(43681,43680)
    # Compose reset, pointer, conditional-X, preparations and inverse words.
    # e bounds each complete source block on the declared proper-code service.
    e=F(1,10**11);t=F(35,2**31);blocks=32+16+4
    gain=(1+e)**(2*blocks)
    complete=(gain-1)/2+gain*32*t
    assert complete<F(6,10**7)
    X=np.eye(6,dtype=complex);X[:2,:2]=np.array([[0,1j],[1j,0]])
    assert np.linalg.norm(arrays['full_M6_X']-X,'fro')<=float(e)
    result={'schema':'oph.localized_signed_rest.acceptance_audit.v1','all_pass':True,
        'native_calls':0,'recomputed_complete_laws':10,'labels_per_law':65536,
        'probability_array_max_error':max_error,'auxiliary_tape_endpoint_checks':checks,
        'exact_finite_clock_excess':'1260*log(43681/43680)',
        'complete_prepare_and_energy_read_error_upper':'6/10000000',
        'complete_error_exact_composition':str(complete),
        'proper_code_services_and_precision_supplied':True,
        'physical_hardware_error_authenticated':False,
        'elapsed_seconds':time.perf_counter()-started}
    (HERE/'ACCEPTANCE_AUDIT.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='complete_error_exact_composition'},indent=2,sort_keys=True))

if __name__=='__main__':main()
