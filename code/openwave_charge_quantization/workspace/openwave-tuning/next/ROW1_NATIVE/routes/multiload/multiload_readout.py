"""Two-load literal-C adapter. Frozen one-load pointer is reused unchanged."""
import hashlib, importlib.util, json, inspect
from pathlib import Path
from fractions import Fraction
import numpy as np
HERE=Path(__file__).resolve().parent
METER=HERE.parent/'auditor_bridge/native_cycle_meter.py'
spec=importlib.util.spec_from_file_location('frozen_cycle_meter', METER)
meter=importlib.util.module_from_spec(spec);spec.loader.exec_module(meter)
native=meter.native; EDGES=meter.EDGES

def require(ok,msg):
    if not ok: raise ValueError(msg)

def decode(residue):
    require(type(residue) is int and 0<=residue<6,'six-state class required')
    require(residue!=3,'class3 outside two-load interval; no arbitrary integer lift')
    return residue if residue<=2 else residue-6

def preparation(G):
    require(len(G)==12 and all(type(x) is int and x>=0 for x in G) and sum(G)==2,'exact nonnegative total-two preparation required')
    return tuple(G)

def execute(G,schedule,source_coins,region,pointer_initial,pointer_coins,retries=16):
    G=preparation(G); n=G; counters=[0]*12; records=[]
    require(len(schedule)==len(source_coins),'complete scheduled source coins required')
    region=tuple(region)
    require(len(set(region))==len(region) and all(type(v) is int and 0<=v<12 for v in region),'valid region required')
    for e,coin in zip(schedule,source_coins):
        require(type(e) is int and 0<=e<len(EDGES) and type(coin) is int and coin in (0,1),'native seam/coin alphabet required')
        a,b=EDGES[e];n,(_,j)=native(n,(a,b),coin)
        counters[a]=(counters[a]+j)%6;counters[b]=(counters[b]-j)%6
        records.append({'seam':e,'coin':coin,'current':j})
    # Controller uses only current-derived finite memory, never terminal loads.
    residue=sum(counters[v] for v in region)%6;w=decode(residue)
    result=meter.execute(pointer_initial,(w,),pointer_coins,retries)
    if result['success']:
        require(type(result['winding']) is int and result['winding']==w,'pointer output disagrees with commanded winding')
    p=meter.reset_success(pointer_initial,retries)*(1-Fraction(1,2**retries))**(6*abs(w))
    return {'source_records':records,'counter_class':residue,'meter':result,
            'decoded_spatial_degree':result['winding'] if result['success'] else None,
            'decoded_relative_charge':-result['winding'] if result['success'] else None,
            'conditional_meter_success':str(p),'conditional_meter_timeout':str(1-p),
            'combined_attempts':len(records)+len(result['records']),
            'fine_joint_probability':str(Fraction(1,2**(len(records)+len(result['records'])))),
            'audit_only':{'initial_G':G,'final_n':n,'expected_w':sum(G[v]-n[v] for v in region)}}

def run():
    states=[]
    for a in range(12):
        for b in range(a,12):
            n=[0]*12;n[a]+=1;n[b]+=1;states.append(tuple(n))
    cuts=np.array([[(m>>v)&1 for v in range(12)] for m in range(4096)],dtype=np.int64)
    occup=np.asarray(states,dtype=np.int64)@cuts.T
    pair_checks=0
    for i in range(78):
        differences=occup[i]-occup
        residues=differences%6
        decoded=np.where(residues<=2,residues,residues-6)
        require(np.array_equal(decoded,differences),'bounded decoder alias')
        require(not np.any(residues==3),'forbidden class reached')
        pair_checks+=differences.size
    events=0;moves=0;max_current=0
    for G in states:
        for e,(a,b) in enumerate(EDGES):
            for coin in (0,1):
                n,(_,j)=native(G,(a,b),coin);require(sum(n)==2 and min(n)>=0,'closure')
                delta=np.zeros(12,dtype=np.int64);delta[a]=j;delta[b]=-j
                require(np.array_equal(np.array(G)-np.array(n),delta),'actual native incidence')
                require(np.array_equal(cuts@delta,cuts@(np.array(G)-np.array(n))),'all-cut current')
                events+=1;moves+=int(j!=0);max_current=max(max_current,abs(j))
    examples=[]
    for w in range(-2,3):
        for x in range(12):
            trace=meter.execute(x,(w,),meter.favorable_coins(x,(w,),16),16)
            require(trace['success'] and trace['winding']==w,'signed pointer')
            examples.append({'w':w,'pointer_initial':x,'attempts':len(trace['records']),'fine_probability':trace['fine_history_probability']})
    # Concrete current-derived double-charge readout: (2,0)->(1,1).
    G=[0]*12;G[0]=2;e=next(i for i,ab in enumerate(EDGES) if 0 in ab);a,b=EDGES[e]
    n,(_,j)=native(tuple(G),(a,b),0);region=(b,);w=G[b]-n[b]
    joined=execute(G,[e],[0],region,7,meter.favorable_coins(7,(w,),16),16)
    require(joined['decoded_spatial_degree']==joined['audit_only']['expected_w'],'joined readout')
    # Complete two-attempt source law, including holds; actual degree two is reached.
    source_edges=[i for i,ab in enumerate(EDGES) if 0 in ab][:2]
    complete=[]
    for c0 in (0,1):
        for c1 in (0,1):
            n=tuple(G)
            for ei,c in zip(source_edges,(c0,c1)):n=native(n,EDGES[ei],c)[0]
            w=G[0]-n[0]
            r=execute(G,source_edges,[c0,c1],(0,),4,meter.favorable_coins(4,(w,),16),16)
            require(r['decoded_spatial_degree']==w,'actual degree-two sequence')
            r['source_probability']='1/4';complete.append(r)
    require(any(r['decoded_spatial_degree']==2 for r in complete),'degree two absent')
    require(sum(Fraction(r['source_probability']) for r in complete)==1,'source law incomplete')
    timeout=execute(G,[],[],(),7,iter([0]*2000),16)
    require(timeout['decoded_relative_charge'] is None,'timeout not zero')
    worst=max(1-meter.reset_success(x,16)*(1-Fraction(1,2**16))**12 for x in range(12))
    return {'schema':'oph.native.two_load_readout.v1','checks':{'all_two_load_preparations':len(states)==78,'all_native_branches':events==4680,'all_regions_and_preparation_pairs_injective':pair_checks==24920064,'all_signed_pointer_replays':len(examples)==60,'joined_actual_current_readout':True,'timeout_not_zero':True,'alias_plus3_minus3':3%6==(-3)%6},'state_count':78,'native_events':events,'moving_events':moves,'max_current':max_current,'regional_pair_checks':pair_checks,'pointer_replays':examples,'joined_example':joined,'complete_two_event_source_branches':complete,'timeout_example':timeout,'complete_conditional_law':'p_success(x,w)=p_reset(x)*(1-2^-R)^(6 abs(w)); timeout=1-p_success; every fine source+pointer history has 2^-combined_attempts','worst_timeout_R16':str(worst),'worst_timeout_float':float(worst),'meter_attempt_cap':(30+3+12)*16,'source_contract':'Exactly two nonnegative classical loads, known G, zero initial counters; arbitrary finite native source schedule. Exhaustive G,n pairs are decoder overapproximation, not asserted reachable histories.','row8_boundary':'Exact total load conservation and relative excess/deficit only; no annihilation-to-radiation or signed-species process established.','force_boundary':'Native C does not read F; no force feedback inferred.','official_rows_changed':[],'native_function_sha256':hashlib.sha256(Path(inspect.getfile(native)).read_bytes()).hexdigest(),'geometry_sha256':hashlib.sha256(Path(inspect.getfile(__import__('oph_openwave.geometry',fromlist=['CARRIER_SEAMS']))).read_bytes()).hexdigest(),'source_sha256':{str(p.relative_to(HERE.parent.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),METER)}}

if __name__=='__main__':
    out=run();(HERE/'MULTILOAD_READOUT.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out['checks']))
