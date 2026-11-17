#!/usr/bin/env python3
"""Preregistered unified-controller replay and valid-current biased-law control."""
import gzip,hashlib,json,math,random
from pathlib import Path
import numpy as np
from instrument import Instrument,Settings,native,INTERNAL

HERE=Path(__file__).resolve().parent
STEPS=32768
BIAS_SEED=2026111435

def enc(x):return (json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode()

class OuterCoinZero:
    """A different random-input law; not an allowed preparation of the source."""
    def __init__(self,base):self.base=base;self.draw=0
    def randrange(self,*args):
        self.draw+=1;value=self.base.randrange(*args)
        return 0 if self.draw==2 else value
    def choice(self,x):return self.base.choice(x)

def main():
    oldpath=HERE/'inputs/PRODUCTION_RECORDS.json.gz'
    assert hashlib.sha256(oldpath.read_bytes()).hexdigest()=='f36f7e8e7b1c179632fb70fdc17de0571431ee423842e4e5f9705fc49121a001'
    old=json.loads(gzip.decompress(oldpath.read_bytes()))
    ins=Instrument(Settings());N=12*len(ins.outer['cells']);W=ins.proposal_weight
    assert N==972 and W==6318
    census_calls=full_calls=0;fair_moves=np.zeros((N,N),dtype=np.int64)
    for name,g in [('outer',ins.outer),('inner',ins.inner)]:
        assert [list(x) for x in g['cells']]==old[name+'_graph']['cells']
        assert [list(x) for x in g['edges']]==old[name+'_graph']['edges']
        for eid,slot,coin,current,target in old['branches'][name]:
            x=[0,0,0];x[slot]=1;y,rec=native(x,(0,1),coin);census_calls+=1
            assert rec==(coin,current) and y.index(1)==target
            a,b,_=g['edges'][eid]
            if name=='outer' and slot<2 and target!=slot:fair_moves[(a,b)[slot],(a,b)[target]]+=1
            if eid<64:
                spectator=next(i for i in range(3) if i not in (a,b));mapping=(a,b,spectator)
                xx=[0]*(12*len(g['cells']));xx[mapping[slot]]=1
                yy,rr=native(xx,(a,b),coin);full_calls+=1
                assert rr==rec and yy.index(1)==mapping[target]
    assert (census_calls,full_calls)==(82008,768)
    reproduced=0;native_comparison_calls=census_calls+full_calls;fair_downward=0;band_error=0.
    lf=np.diag(fair_moves.sum(axis=1))-fair_moves
    lint=lf[:12,:12]-8*np.eye(12,dtype=np.int64)
    ev=np.linalg.eigvalsh(lint);bands=[]
    for x in ev:
        if not bands or abs(x-bands[-1])>1e-8:bands.append(float(x))
    projectors=[]
    for x in bands:
        p=np.eye(12)
        for y in bands:
            if abs(x-y)>1e-8:p=p@(lint-y*np.eye(12))/(x-y)
        assert np.max(abs(p@p-p))<1e-12
        projectors.append(p)
    assert np.max(abs(sum(projectors)-np.eye(12)))<1e-12
    for run in old['runs']:
        rng=random.Random(run['seed']);coupled=not run['disabled']
        for channel in run['channels']:
            start=rng.randrange(N);assert start==channel['start'];memory=ins.prepare(start,coupled)
            native_before=memory.native_calls;offset=0;counts=[0]*N;clock=[]
            previous=start
            for step,expected in enumerate(channel['macros'],1):
                out=ins.step(memory,rng);end=offset+len(out['inner_events'])
                actual=[out['proposal'],out['coin'],out['current'],int(out['auxiliary_hold']),offset,end,
                        out['preparation_calls'],out['probe_calls'],out['outcome']]
                assert actual==expected
                assert out['inner_events']==channel['inner_events'][offset:end]
                assert out['target_coupled']==coupled and memory.target_coupled==coupled
                if out['current'] and memory.outer_position%12<previous%12:fair_downward+=1
                previous=memory.outer_position;counts[previous]+=1
                if out['tick'] is not None:clock.append(out['tick'])
                offset=end;reproduced+=1
            assert offset==len(channel['inner_events']) and counts==channel['occupation_counts'] and clock==channel['clock']
            y=math.sqrt(N/(4*W*STEPS))*(np.asarray(counts)-STEPS/N)
            for p,expected in zip(projectors,channel['band_readouts']):
                error=float(np.max(abs(y.reshape(-1,12)@p-np.asarray(expected))))
                band_error=max(band_error,error);assert error<1e-12
            native_comparison_calls+=memory.native_calls-native_before
    assert reproduced==262144 and native_comparison_calls==347533 and fair_downward>0

    # Changed endpoint-coin law, derived from actual native branch calls.
    biased=np.zeros((N,N),dtype=np.int64);biased_branch_calls=0
    for eid,(a,b,axis) in enumerate(ins.outer['edges']):
        for slot in (0,1):
            x=[0,0];x[slot]=1;y,rec=native(x,(0,1),0);biased_branch_calls+=1
            target=y.index(1)
            assert rec==(0,x[0]-y[0])
            if target!=slot:biased[(a,b)[slot],(a,b)[target]]+=1
    assert biased_branch_calls==12636
    directed=np.diag(biased.sum(axis=1))-biased
    icount=np.zeros((12,12),dtype=np.int64)
    for a,b in INTERNAL:icount[a,b]=biased[a,b]
    sinks=np.flatnonzero(icount.sum(axis=1)==0).tolist();assert sinks==[9,11]
    internal_eigenvalues=np.sort(2*icount.sum(axis=1)).tolist()
    cellbasis=np.kron(np.eye(81,dtype=np.int64),np.ones((12,1),dtype=np.int64))
    cellall=directed@cellbasis;cell=cellall[::12]
    assert np.array_equal(cellall,np.repeat(cell,12,axis=0))
    assert np.all(cell.sum(axis=0)==0) and np.all(cell.sum(axis=1)==0)
    controls=[]
    for momentum in [(1,1,0,0),(1,-1,0,0)]:
        p=2*np.pi*np.asarray(momentum)/3
        v=np.asarray([np.exp(1j*np.dot(p,z)) for z in ins.outer['cells']])/math.sqrt(81)
        pv=v-cell@v/W;z=complex(np.vdot(v,pv))
        assert np.linalg.norm(pv-z*v)<1e-12
        aa=float(np.sum(1-np.cos(p)));bb=float(np.sum(np.sin(p)))
        cb=(2*(1/(1-z)).real-1)/(4*W)
        target=aa/(2*(aa*aa+bb*bb))-1/(4*W)
        assert abs(cb-target)<1e-10
        cf=1/(2*aa)-1/(4*W)
        # Same local F_12 reader: derivative norm q1²+q2²=6.
        derivative=4*(math.sin(p[0]/2)**2+math.sin(p[1]/2)**2)
        controls.append(dict(momentum=momentum,biased_cell_eigenvalue=[z.real,z.imag],
            fair_covariance=cf,biased_covariance=cb,biased_formula=target,
            fair_F12_covariance=derivative*cf,biased_F12_covariance=derivative*cb,
            contact_adjusted_ratio=(cb+1/(4*W))/(cf+1/(4*W))))
    assert abs(controls[0]['contact_adjusted_ratio']-.75)<1e-10
    assert abs(controls[1]['contact_adjusted_ratio']-1)<1e-10
    assert controls[1]['biased_F12_covariance']-controls[0]['biased_F12_covariance']>.2

    rng=random.Random(BIAS_SEED);channels=[];biased_native_calls=0;biased_downward=0;biased_internal_moves=0
    for _ in range(4):
        start=rng.randrange(N);memory=ins.prepare(start,False);trace=[];counts=[0]*N;clock=[];previous=start
        for step in range(1,STEPS+1):
            out=ins.step(memory,OuterCoinZero(rng));assert out['coin']==0
            if out['current'] and memory.outer_position%12!=previous%12:
                biased_internal_moves+=1
                if memory.outer_position%12<previous%12:biased_downward+=1
            previous=memory.outer_position;counts[previous]+=1
            if out['tick'] is not None:clock.append(out['tick'])
            trace.append(out)
        reader=start;inner_reader=12*ins.inner['ci'][(0,0,0)];recount=[0]*N;reclock=[];ticks=0;origin=ins.outer['cells'][start//12][3]
        for step,out in enumerate(trace,1):
            for eid,coin,jj in out['inner_events']:
                ia,ib,_=ins.inner['edges'][eid]
                if jj:assert inner_reader==(ia if jj==1 else ib);inner_reader=ib if jj==1 else ia
            a,b,_=ins.outer['edges'][out['proposal']];j=out['current']
            if j:assert reader==(a if j==1 else b);reader=b if j==1 else a
            recount[reader]+=1
            if ins.outer['cells'][reader//12][3]==(origin+ticks+1)%3:
                ticks+=1;reclock.append([step,ticks,reader])
        assert recount==counts and reclock==clock and inner_reader==memory.inner_position
        assert memory.native_calls==sum(1+len(o['inner_events']) for o in trace)
        biased_native_calls+=memory.native_calls
        channels.append(dict(start=start,target_coupled=False,trace=trace,occupation_counts=counts,clock=clock,
            final_open_interval=STEPS-(clock[-1][0] if clock else 0),native_calls=memory.native_calls,
            censors=[[i+1,o['outcome']] for i,o in enumerate(trace) if 'timeout' in o['outcome']]))
    assert biased_downward==0 and biased_internal_moves>0
    ys=np.array([math.sqrt(N/(4*W*STEPS))*(np.asarray(c['occupation_counts'])-STEPS/N) for c in channels])
    af=(ys.reshape(4,81,12)@np.ones(12)/math.sqrt(12)).reshape(4,3,3,3,3);dirs=(3,0,1,2)
    curl=[(np.roll(af[nu],-1,dirs[mu])-af[nu]-np.roll(af[mu],-1,dirs[nu])+af[mu]).tolist()
          for mu in range(4) for nu in range(mu+1,4)]
    raw=enc(dict(kind='different outer-coin0 law',seed=BIAS_SEED,
        source_preparation='uncoupled c=0; fair inner coins; changed outer coin law',channels=channels,
        field_strength_readouts=curl))
    (HERE/'CONTROL_RECORDS.json.gz').write_bytes(gzip.compress(raw,mtime=0))
    result=dict(candidate='MFR-unified-native-field-v1',source_comparison_native_calls=native_comparison_calls,
        source_macro_histories_exactly_reproduced=reproduced,source_preparation_bits=[True,False],
        source_graphs_and_all_events_identical=True,corrected_projector_polynomial_error=band_error,
        initial_decoder_failure_retained_in_ancestor=True,
        biased_census_native_calls=biased_branch_calls,biased_tape_native_calls=biased_native_calls,
        total_native_calls=native_comparison_calls+biased_branch_calls+biased_native_calls,
        biased_absorbing_internal_ports=sinks,biased_internal_generator_eigenvalues=internal_eigenvalues,
        fair_observed_downward_port_moves=fair_downward,biased_observed_downward_port_moves=biased_downward,
        biased_observed_internal_moves=biased_internal_moves,
        changed_law_same_F12_population=controls,
        changed_population_is_actual_precision_estimate=False,
        biased_clock_ticks=[len(c['clock']) for c in channels],
        biased_final_open_intervals=[c['final_open_interval'] for c in channels],
        biased_censors=[c['censors'] for c in channels],
        all_registered_checks='pass',official='untested',
        moving_target_covariance_established=False,minimal_electric_coupling_established=False,
        control_raw_uncompressed_sha256=hashlib.sha256(raw).hexdigest())
    (HERE/'RECEIPT.json').write_bytes(enc(result))
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=='__main__':main()
