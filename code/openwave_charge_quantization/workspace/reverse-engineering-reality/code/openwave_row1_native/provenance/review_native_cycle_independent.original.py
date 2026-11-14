"""Independent absorbing-stage probability and source-record cycle-meter audit."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,importlib.util,json,sys
import numpy as np
import networkx as nx
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[3]
sys.path[:0]=[str(ROOT/'reverse-engineering-reality'),str(ROOT/'openwave-tuning/candidate')]
from model.observer.native import native
from oph_openwave.geometry import CARRIER_SEAMS
SOURCE=HERE.parent/'routes/auditor_bridge/native_cycle_meter.py'
spec=importlib.util.spec_from_file_location('meter_under_test',SOURCE);meter=importlib.util.module_from_spec(spec);spec.loader.exec_module(meter)
EDGES=tuple(CARRIER_SEAMS);GRAPH=nx.Graph();GRAPH.add_edges_from(EDGES)


def absorbing_reset(initial,retries):
    """Propagate actual fair-coin unresolved mass; no producer reset helper."""
    unresolved=Q(1);discoveries={};attempt_mass=Q()
    loads=tuple(int(v==initial) for v in range(12))
    for edge in EDGES:
        for attempt in range(retries):
            attempt_mass+=unresolved
            nextmass=Q()
            for coin in (0,1):
                state,(_,current)=native(loads,edge,coin)
                if current:
                    destination=state.index(1)
                    discoveries[destination]=discoveries.get(destination,Q())+unresolved/2
                else:
                    if state!=loads:raise ValueError('silent move without record')
                    nextmass+=unresolved/2
            unresolved=nextmass
    # Geometric bounded gambler for each requested one-load hop.
    # p0=0, p(r+1)=1/2 + p(r)/2; every failed step is retained.
    hop=Q()
    for _ in range(retries):hop=(1+hop)/2
    success=Q();routefailure=Q()
    for destination,weight in discoveries.items():
        d=nx.shortest_path_length(GRAPH,destination,0)
        success+=weight*hop**d
        routefailure+=weight*(1-hop**d)
    if success+routefailure+unresolved!=1:raise ValueError('reset not complete')
    return dict(success=success,discovery_timeout=unresolved,routing_timeout=routefailure,hop_success=hop,expected_discovery_attempts=attempt_mass)


def independently_validate_trace(initial,result):
    state=tuple(int(v==initial) for v in range(12));known=None;steps=0
    for row in result['records']:
        state,(_,j)=native(state,tuple(row['seam']),row['coin'])
        if j!=row['current']:raise ValueError('record differs from actual native law')
        if row['stage']!='discover' and known is None:raise ValueError('route before observed position')
        if j:known=row['seam'][1] if j>0 else row['seam'][0]
        if row['stage']=='measure' and j:
            steps+=row['requested_cycle_direction']
    if list(state)!=result['final_loads']:raise ValueError('hidden final state mismatch')
    if Q(result['fine_history_probability'])!=Q(1,2**len(result['records'])):raise ValueError('history weight mismatch')
    if result['success']:
        if known!=meter.CYCLE[0] or steps%6 or result['winding']!=steps//6:raise ValueError('open/inconsistent pointer loop')
    elif result['winding'] is not None:raise ValueError('timeout assigned a charge')


def run():
    # Source/current counter domains independently exhaust all G,n and regions.
    regions=np.array([[int(mask>>v&1) for v in range(12)] for mask in range(4096)])
    decode={0:0,1:1,5:-1};counterchecks=0;sourceevents=0
    for origin in range(12):
        G=np.eye(12,dtype=int)[origin]
        for occupied in range(12):
            n=np.eye(12,dtype=int)[occupied];q=G-n;c=q%6
            actual=regions@q;classes=(regions@c)%6
            measured=np.array([decode[int(x)] for x in classes])
            if not np.array_equal(actual,measured):raise ValueError('bounded regional decoder alias')
            counterchecks+=4096
            for e,(a,b) in enumerate(EDGES):
                for coin in (0,1):
                    updated,(_,j)=native(tuple(map(int,n)),(a,b),coin)
                    nextc=c.copy();nextc[a]=(nextc[a]+j)%6;nextc[b]=(nextc[b]-j)%6
                    if not np.array_equal(nextc,(G-np.array(updated))%6):raise ValueError('counter/current mismatch')
                    sourceevents+=1
    # Exact unconditional law, all pointer starting ports, several retry caps.
    rows=[];max_reset_error=Q();max_timeout_error=Q()
    for retries in (1,2,4,16):
        for initial in range(12):
            distribution=absorbing_reset(initial,retries)
            reported=meter.reset_success(initial,retries)
            max_reset_error=max(max_reset_error,abs(distribution['success']-reported))
            for signedsteps in (0,6,12,18):
                success=distribution['success']*distribution['hop_success']**signedsteps
                measurementfailure=distribution['success']*(1-distribution['hop_success']**signedsteps)
                failure=distribution['discovery_timeout']+distribution['routing_timeout']+measurementfailure
                if success+failure!=1:raise ValueError('unconditional probability omitted')
                max_timeout_error=max(max_timeout_error,abs(failure-(1-reported*distribution['hop_success']**signedsteps)))
                if retries==16 and signedsteps==6:
                    rows.append(dict(initial_port=initial,success=str(success),discovery_timeout=str(distribution['discovery_timeout']),routing_timeout=str(distribution['routing_timeout']),measurement_timeout=str(measurementfailure),total_timeout=str(failure)))
    rng=np.random.default_rng(1819);replays=0;failures=0;successes=0
    for retries in (1,2,4):
        for initial in range(12):
            for turns in ((0,),(1,),(-1,),(1,-1),(2,)):
                # Supplied finite coin tape is sampled only for trace checking;
                # exact probabilities above are not inferred from frequencies.
                coins=iter(map(int,rng.integers(2,size=4096)))
                result=meter.execute(initial,turns,coins,retries)
                independently_validate_trace(initial,result)
                if result['success'] and result['winding']!=sum(turns):raise ValueError('wrong nonconditioned success readout')
                replays+=1;successes+=int(result['success']);failures+=int(not result['success'])
    # Explicit favorable traces for every unknown port and each bounded class.
    for initial in range(12):
        for cls,winding in ((0,0),(1,1),(5,-1)):
            result=meter.execute_class(initial,cls,meter.favorable_coins(initial,(winding,),2),2)
            independently_validate_trace(initial,result)
            if not result['success'] or result['winding']!=winding:raise ValueError('bounded meter replay')
    # Positive complement labels are refused, not silently mapped to halfcharge.
    refused=[meter.execute_class(0,c,iter(()),2) for c in (2,3,4)]
    alias=decode[7%6]!=7
    checks=dict(independent_absorbing_reset_matches=max_reset_error==0,unconditional_stage_failures_complete=max_timeout_error==0,all_unknown_pointer_ports_covered=len(rows)==12,actual_native_trace_replay=replays==180 and failures>0 and successes>0,bounded_source_all_regions=counterchecks==589824,bounded_counter_all_native_events=sourceevents==8640,complement_labels_retained=all(not r['success'] and r['winding'] is None for r in refused),larger_charge_alias_rejected=alias,no_coherent_pulses_in_meter='compiler' not in SOURCE.read_text() and 'native_pulse' not in SOURCE.read_text())
    files=[Path(__file__),SOURCE,ROOT/'reverse-engineering-reality/model/observer/native.py',ROOT/'openwave-tuning/candidate/oph_openwave/geometry.py']
    receipt=SOURCE.with_name('NATIVE_CYCLE_METER.json')
    if receipt.exists():
        files.append(receipt)
        recorded=json.loads(receipt.read_text())
        actual_worst=max(Q(r['total_timeout']) for r in rows)
        if Q(recorded['compressed_primary_meter_worst_timeout'])!=actual_worst:
            raise ValueError('archived unconditional bounded timeout differs')
        for initial,r in enumerate(rows):
            hop=Q(1)-Q(1,2**16)
            if Q(recorded['exact_reset_success_by_unknown_initial_port'][initial])!=Q(r['success'])/hop**6:
                raise ValueError('archived exact reset differs')
        for name,expected in recorded['source_sha256'].items():
            if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=expected:
                raise ValueError('producer receipt source pin mismatch: '+name)
        files.append(HERE.parent/'routes/spatial/NATIVE_DUAL_WINDING.json')
        checks['archived_probabilities_and_source_custody_match']=True
    return dict(schema='oph.row1.native_cycle_independent_review.v1',reviewer='independent Codex worker reviewing external auditor source; not a Fable opposing-team review',checks=checks,source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},exact_unconditional_bounded_primary_law=rows,max_reset_probability_error=str(max_reset_error),max_complete_timeout_error=str(max_timeout_error),regional_comparisons=counterchecks,native_counter_events=sourceevents,trace_replays=replays,retained_timeout_replays=failures,successful_replays=successes,worst_bounded_primary_timeout=float(max(Q(r['total_timeout']) for r in rows)),assessment='Native-only apparatus removes proper-code/analog-phase/ideal-pointer premises. It reads an actual source-record regional class with an additional literal one-load native carrier and finite record-adaptive controller. The physical pointer path is apparatus; the spatial homology/current theorem supplies what its winding represents. Source controller access, retained G and bounded one-load decoder are explicit. No success postselection or guessed timeout charge.',official_rows_changed=[])

if __name__=='__main__':
    r=run()
    if not all(r['checks'].values()):raise ValueError(r['checks'])
    (HERE/'NATIVE_CYCLE_INDEPENDENT_REVIEW.json').write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps({k:r[k] for k in ('checks','worst_bounded_primary_timeout','trace_replays','retained_timeout_replays')},indent=2))
