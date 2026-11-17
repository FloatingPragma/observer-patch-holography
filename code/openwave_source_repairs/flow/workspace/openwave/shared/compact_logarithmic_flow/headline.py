#!/usr/bin/env python3
"""Small independent arithmetic/record check; no producer/provider imports.

This does not replay millions of source pulses. It reconstructs source
transfer, rational certificate and ideal flow, then verifies all saved
complete detector arrays against their unnormalized source outputs.
"""
import argparse,hashlib,json,time
from fractions import Fraction as F
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def check(attempt='stage2_001'):
    started=time.monotonic();output=HERE/'attempts'/attempt
    receipt=json.loads((output/'RECEIPT.json').read_text());science=receipt['scientific']
    for rel,digest in receipt['artifacts'].items():assert sha(output/rel)==digest,rel
    for pin in ('SOURCE_PINS.json','PREREG_PROPOSAL.json','SOURCE_PINS_FLOW.json','PREREG_FLOW.json'):
        for rel,digest in json.loads((HERE/pin).read_text())['files'].items():assert sha(ROOT/rel)==digest,rel
    source=np.load(HERE/'source_columns.npz');proposal=json.loads((HERE/'proposal.json').read_text())
    L=np.zeros((240,240),dtype=np.int64)
    for a,b in source['seams']:L[a,a]+=1;L[b,b]+=1;L[a,b]-=1;L[b,a]-=1
    C=1260*np.eye(240,dtype=np.int64)-L
    assert np.array_equal(C,source['counts']) and int(np.trace(C))==301140
    f=np.zeros(240);g=np.zeros(240);f[[4,5,6,7]]=[.5,-.5,-.5,.5];g[[16,17,18,19]]=[.5,-.5,-.5,.5]
    assert np.array_equal(L@f,6*f) and np.array_equal(L@g,6*g)
    changed=L.copy();a,b=proposal['held_seam'];changed[a,a]-=1;changed[b,b]-=1;changed[a,b]+=1;changed[b,a]+=1
    assert np.array_equal(1260*np.eye(240,dtype=np.int64)-changed,source['changed_counts'])
    independent={}
    for law,matrix in [('native',C),('held_seam',source['changed_counts'])]:
        Q=proposal['denominator'];rows=proposal['laws'][law]
        V=np.array(rows['columns'],dtype=object).T;lam=np.array(rows['eigenvalues'],dtype=object)
        gram=V.T@V-Q*Q*np.eye(240,dtype=object)
        residual=Q*np.array(matrix,dtype=object)@V-1260*V*lam
        gbound=F(max(sum(abs(int(z)) for z in row) for row in gram),Q*Q)
        rbound=F(max(max(sum(abs(int(z)) for z in row) for row in residual),max(sum(abs(int(z)) for z in row) for row in residual.T)),1260*Q*Q)
        floor=F(1260-2*max(1260-int(matrix[i,i]) for i in range(240)),1260)
        assert gbound<=F(1,10**9) and rbound<=F(1,10**9)
        assert all(floor<=F(int(v),Q)<=1 for v in lam)
        cert=science['spectral_certificates'][law]
        assert F(cert['gram_operator_upper'])==gbound and F(cert['residual_operator_upper'])==rbound
        assert F(cert['energy_operator_upper'])==1260*(rbound+2*gbound)/floor
        independent[law]={'gram':float(gbound),'residual':float(rbound)}
    arrays=np.load(output/'evidence/flow_records.npz');initial=arrays['initial_source_columns']
    assert np.max(abs(initial[:240,0]-1j*f))<1e-12 and np.max(abs(initial[:240,1]-1j*g))<1e-12
    for law,matrix,times in [('native',L,[1,4]),('held_seam',changed,[4])]:
        eigen,vectors=np.linalg.eigh(matrix.astype(float));energies=-1260*np.log1p(-eigen/1260)
        for i,t in enumerate(times):
            state=arrays[law+'_terminal_source'][:,2*i:2*i+2]
            expected=vectors@((np.exp(-1j*t*energies)[:,None])*(vectors.T@initial[:240]))
            residual=float(np.linalg.norm(state[:240]-expected,'fro'))
            assert residual<1e-7,(law,t,residual)
            position=np.einsum('zi,zj->zij',state.conj(),state)
            key=law+'_t'+str(t)
            for reader in ('C','E'):
                assert np.max(abs(arrays[key+'_position_'+reader]-position))<1e-9
            outside=float(np.sum(abs(state[12:,0])**2));click=float(abs(np.vdot(initial[:,0],state[:,0]))**2)
            if law=='native':assert outside<1e-7 and click>1-1e-7
            else:assert outside>1e-3
            recorded=science['complete_pointer_records'][key+'_C']
            assert abs(recorded['outside_original_cell']-outside)<1e-9 and abs(recorded['energy_ray_click']-click)<1e-9
            independent[key]={'outside':outside,'f_ray':click,'source_vs_independent_flow_error':residual}
    controls=np.load(output/'evidence/coherence_controls.npz');contrast={}
    for name in ('coherent','which_cell_dephased','reset_to_f'):
        v=[]
        for setting in ('plus','minus'):
            record=controls[name+'_'+setting];assert record.shape==(256,2) and np.min(record)>-1e-12
            assert abs(float(record.sum())-1)<1e-7
            v.append(float(record[0,0]-record[0,1]))
        contrast[name]=v[0]-v[1]
    assert contrast['coherent']>1-1e-7
    assert abs(contrast['which_cell_dephased'])<1e-7 and abs(contrast['reset_to_f'])<1e-7
    assert all(abs(v-science['coherence_controls']['source_reference_contrast'][k])<1e-12 for k,v in contrast.items())
    reference=json.loads((output/'evidence/REFERENCE_PREPARATION.json').read_text())
    assert F(reference['ideal_plus_given_selectors'])==F(239,240) and F(reference['ideal_minus_given_selectors'])==F(1,240)
    assert reference['full_clean_permutation_columns']==630*512
    assert F(reference['complete_preparation_halfdiamond_upper'])<F(1,10**4)
    elapsed=time.monotonic()-started;assert elapsed<60,elapsed
    result={'independent_arithmetic':independent,'coherence_contrast':contrast,'seconds':elapsed,
            'full_source_word_replay_claimed':False,'official_row_award':False}
    print(json.dumps(result,indent=2,allow_nan=False));return result

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--attempt',default='stage2_001');args=parser.parse_args();check(args.attempt)
