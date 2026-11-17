#!/usr/bin/env python3
"""Certified first-cell reuse of the full-source native word implementation.

This is an operator-certificate transfer to another input subspace. It does
not count its interval matrix algebra as additional native event executions.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,math,sys,time
import numpy as np
HERE=Path(__file__).resolve().parent;OW=HERE.parents[2];FLOW=OW/'shared/compact_logarithmic_flow'
sys.path.insert(0,str(FLOW));sys.path.insert(0,str(FLOW/'inputs'))
import certify,intervals as iv,native_words as nw,reference

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rotations(angle):
    s,c=iv.sincos(F(angle)/16)
    for _ in range(4):c,s=iv.add(iv.mul(c,c),iv.neg(iv.mul(s,s))),iv.mul(iv.exact(2),iv.mul(s,c))
    return c,s
def main():
    started=time.monotonic();reg=json.loads((HERE/'PREREG_FIRST_CELL.json').read_text())
    for name,digest in reg['files'].items():assert sha(HERE/name)==digest
    for name,digest in reg['flow_files'].items():assert sha(FLOW/name)==digest
    proposal,source=certify.load();cert=certify.certify('native');rows=proposal['laws']['native']
    V=np.asarray(rows['columns'],dtype=object).T;den=proposal['denominator']**2*iv.Q
    phases=[rotations(certify.phase(e,4)[0]) for e in rows['eigenvalues']]
    Q=np.zeros((12,12),complex);entry=F()
    for a in range(12):
        for b in range(12):
            sums=[[0,0],[0,0]]
            for j,parts in enumerate(phases):
                coeff=int(V[a,j])*int(V[b,j])
                for k,interval in enumerate(parts):
                    vals=[coeff*x for x in interval];sums[k][0]+=min(vals);sums[k][1]+=max(vals)
            values=[]
            for lo,hi in sums:
                x=float(F(lo+hi,2*den));values.append(x)
                entry=max(entry,abs(F(x)-F(lo,den)),abs(F(x)-F(hi,den)))
            Q[a,b]=complex(*values)
    g=F(cert['gram_operator_upper'])
    delta=24*entry+3*g+4*F(cert['energy_operator_upper'])+F(2,10**12)
    D=np.zeros((240,2),dtype=np.int64);D[[4,5,6,7],0]=[1,-1,-1,1];D[[5,6,8,11],1]=[-1,-1,1,1]
    L=1260*np.eye(240,dtype=np.int64)-source['counts'];assert np.array_equal(L@D,6*D)
    local=D[:12];P=local@np.array([[2,-1],[-1,2]])@local.T/6
    assert np.max(abs(P[:,0]))==0
    tail=np.eye(12)-P;deficit=np.eye(12)-Q.conj().T@Q
    for _ in range(64):
        tail=Q.conj().T@tail@Q+np.trace(tail)/12*deficit
        assert np.max(abs(tail))<2
    kappa=4*delta+2*delta*delta;u=F(1,2**53);gamma128=128*u/(1-128*u)
    roundoff=64*12*4096*gamma128/(1-64*kappa)+F(1,10**12)
    frob=float(np.linalg.norm(tail,'fro'));upper=F(math.ceil(frob*10**12),10**12)+roundoff+64*kappa/(1-64*kappa)
    assert upper<F(1,4)
    # The first region predicate is0<=address<12. Each component uses
    # the same native CCX/flip routine already executed for the firstcell.
    ops=[]
    def gate(name,t):ops.append((name,tuple(t)))
    def controlled(ctrl,zeros):
        for z in zeros:gate('flip',(z,))
        scratch=list(range(9,13));chain=[(ctrl[0],ctrl[1],scratch[0])]
        for j in range(2,len(ctrl)-1):chain.append((scratch[j-2],ctrl[j],scratch[j-1]))
        for t in chain:gate('ccx',t)
        gate('ccx',(scratch[len(ctrl)-3],ctrl[-1],8))
        for t in reversed(chain):gate('ccx',t)
        for z in reversed(zeros):gate('flip',(z,))
    controlled([3,4,5,6,7],[3,4,5,6,7]);controlled([2,3,4,5,6,7],[2,4,5,6,7])
    initial=np.arange(512);truth=initial.copy()
    for name,t in ops:
        if name=='flip':truth^=1<<t[0]
        else:truth^=(((truth>>t[0])&1)&((truth>>t[1])&1))<<t[2]
    want=initial^(((initial&255)<12).astype(int)<<8)
    assert np.array_equal(truth,want)
    words,raw,failed,metadata=reference.routines()
    # Three complete13qubit raw witnesses execute the changedaddress words.
    state=np.zeros((8192,3),complex);state[:256,0]=1/16;state[255,1]=1
    state[[4,5,6,7],2]=[.5,-.5,-.5,.5];original=state.copy();failure=np.zeros((3,3),complex)
    for name,t in ops:
        ix=nw.indices(tuple(reversed(t)),13);block=state[ix]
        failure+=np.einsum('iac,ij,jad->cd',block.conj(),failed[name],block,optimize=True)
        state[ix]=(raw[name]@block.reshape(len(raw[name]),-1)).reshape(block.shape)
    wanted=np.zeros_like(state);wanted[want]=original[:512]
    operator=float(np.linalg.norm(state-wanted,2));assert operator<1e-8
    routes=[]
    for p in range(12):
        word=[q for bit in range(8) if p>>bit&1 for q in nw.c.remap(words['flip'],[bit])]
        blank=np.zeros(256);blank[0]=1;v,bad=nw.apply_word(word,blank,8);routes.append(v[:,0])
    route=np.asarray(routes).T;target=np.zeros((256,256));target[:12,:12]=np.eye(12)/12
    assert np.linalg.norm(route@route.conj().T/12-target,2)<1e-9
    pulses=sum(len(words[name]) for name,t in ops)
    reset=F(35,2147483648);pointer=F(200000000001,20000000000000000000000)
    app=21*(reset+pointer)+nw.complete_bound(pulses)['half_diamond']+pointer+8*nw.PULSE_GATE+F(1,4)**16
    result={'candidate':'A-native-poisson-compact-energy-v1','initial_source_address':0,
        'native_flow_program_reused':'the same full240word, with its arbitrary-input/operatorcertificate',
        'fresh_raw_fullflow_columns_executed':False,'certified_subblock_error_upper':str(delta),
        '64epoch_tail_frobenius':frob,'64epoch_tail_operator_upper':str(upper),
        '512epoch_ideal_bad_weight_upper':'1/65536','selected_density':'P_D0/2 in the ideal limit',
        'initial_dark_overlap_exact_zero':True,'512epoch_apparatus_error_upper':str(512*app),
        'all512_region_truth_columns_checked':True,'raw13qubit_witness_error':operator,
        'region_word_pulses':pulses,'native_processor_counts':nw.COUNTS,
        'firstcell_selector_coin_tapes_reused':'sameuniform12selector, sourceaddresses0through11inexecutedroutes',
        'full_row_solved':False}
    out=HERE/'first_cell_certified';out.mkdir(exist_ok=False)
    np.savez_compressed(out/'maps.npz',certified_Q_proposal=Q,dark_integer_basis=D,
        tail64=tail,region_states=state,region_processor_failure=failure,reinjection_routes=route)
    payload={'scientific':result,'artifacts':{'maps.npz':sha(out/'maps.npz')},
        'prereg_sha256':sha(HERE/'PREREG_FIRST_CELL.json'),'elapsed_seconds':time.monotonic()-started}
    (out/'RECEIPT.json').write_text(json.dumps(payload,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'firstcell_pass':True,'tail_upper':float(upper),'apparatus512':float(512*app),
                      'seconds':time.monotonic()-started}))

if __name__=='__main__':main()
