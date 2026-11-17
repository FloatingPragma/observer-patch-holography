"""Read saved trajectories; no producer imports and no source regeneration."""
import argparse
from collections import Counter
from fractions import Fraction as F
import gzip
import hashlib
import json
from pathlib import Path
import random
import time

HERE=Path(__file__).resolve().parent
EDGES=((0,1),(1,2),(0,2))


def require(ok,message):
    if not ok:raise ValueError(message)


def canon(x):return json.dumps(x,sort_keys=True,separators=(',',':'))
def sha(x):return hashlib.sha256(x).hexdigest()
def read(path):return json.loads(gzip.decompress(path.read_bytes()))


def calibration_replay(observed, raw, target, sign, cap=64):
    """Rebuild every controller decision from the raw, ordered record tape."""
    index=0;position=1;phases=0;columns=[];decoded=[]
    def consume(edge,stage,group):
        nonlocal index,position
        require(index<len(raw),'calibration tape truncated')
        got_stage,record=raw[index];index+=1
        require(got_stage==stage and record[:2]==[target[i] for i in edge],'record-only requested action')
        a,b,c,j=record;before=position
        require(type(c) is int and c in (0,1) and type(j) is int,'record alphabet')
        q=[sign*int(i==position) for i in range(3)];q[edge[0]]-=j;q[edge[1]]+=j
        require(q.count(sign)==1 and q.count(0)==2,'record-only unit contradiction')
        position=q.index(sign)
        decoded.append(dict(stage=stage,group=list(group),before_position=before,
                            record=record,after_position=position))
        return c,j
    def route(destination,group):
        nonlocal phases
        phases+=1
        if position==destination:return True
        edge=tuple(sorted((position,destination)))
        for _ in range(cap):
            consume(edge,'calibration_route',group)
            if position==destination:return True
        return False
    outcome='calibrated';unfinished=None;stopped=False
    for p in range(3):
        if stopped:break
        for e,edge in enumerate(EDGES):
            seen={};group=[p,e]
            for _ in range(cap):
                if not route(p,group):outcome='route_timeout';stopped=True;break
                c,j=consume(edge,'calibration_probe',group)
                row=dict(input_position=p,edge=list(edge),coin=c,current=j,decoded_output=position)
                if c in seen:require(seen[c]==row,'inconsistent repeated column')
                else:seen[c]=row;columns.append(row)
                if len(seen)==2:break
            else:outcome='coupon_timeout';stopped=True
            if stopped:unfinished=group;break
    if not stopped and not route(1,['return']):outcome='return_timeout';unfinished=['return']
    if outcome=='calibrated':columns.sort(key=lambda r:(r['input_position'],EDGES.index(tuple(r['edge'])),r['coin']))
    expected=dict(outcome=outcome,target=target,sign=sign,entry_position=1,return_position=position,
                  cap=cap,columns=columns,events=decoded,routing_phases=phases,
                  endpoint_calls=index,unfinished_group=unfinished)
    require(observed==expected,'complete observer calibration differs from record-only reconstruction')
    return index,position,columns


def inspect_trial(trial,old,g,law,seed,forced=None):
    forced=forced or {};audit=trial['audit'];observer=trial['observer'];entry=observer['admission']
    require(entry==old['expected_controller'],'admission observation ancestry')
    prefix=old['expected_audit']['events'];n=len(prefix);events=audit['events']
    require(events[:n]==prefix,'continuing admission prefix')
    q=list(old['initial_audit_only']);allowed={tuple(sorted(e)) for e in g['edges']}
    ar=random.Random(old['admission_seed_audit_only']);cr=random.Random(seed)
    counts={'native':0,'always_swap':0};position=1
    background=list(entry.get('observed_final',q))
    if entry['outcome']=='signed_target_unit':background[g['target'][1]]-=entry['sign']
    for i,(stage,r) in enumerate(events):
        require(len(r)==4 and all(type(v) is int for v in r),'integer native event')
        a,b,c,j=r;require(tuple(sorted((a,b))) in allowed and c in (0,1),'source edge/coin')
        rng=cr if stage.startswith('calibration_') else ar
        expected_coin=forced[stage] if stage in forced else rng.randrange(2)
        require(c==expected_coin,'independent public coin stream')
        active_law=law if stage=='calibration_probe' else 'native'
        d=q[a]-q[b];expected_j=d if active_law=='always_swap' else (d-(2*c-1)*(abs(d)%2))//2
        require(j==expected_j,'literal native current')
        q[a]-=j;q[b]+=j;counts[active_law]+=1
        if i==n-1:require(q==entry.get('observed_final',old['expected_audit']['final']),'actual admission state')
        if i>=n:
            sig=entry['sign'];delta=[q[v]-background[v] for v in g['target']]
            require(delta.count(sig)==1 and delta.count(0)==2,'actual continuing unit')
            position=delta.index(sig);temporary=list(q);temporary[g['target'][position]]-=sig
            require(temporary==background,'full spectator background changed')
    require(audit['final']==q and audit['literal_endpoint_calls']==counts,'full source final/counts')
    if entry['outcome']=='no_unit':
        require(observer=={'admission':entry,'outcome':'not_admitted'} and len(events)==n,'zero terminal')
        return counts,n,0,None
    used,p,columns=calibration_replay(observer['calibration'],events[n:],g['target'],entry['sign'])
    outcome=observer['calibration']['outcome']
    if outcome=='calibrated':
        rest=events[n+used:];coins=set();records=[];blank_out='blank_timeout'
        for i,(stage,r) in enumerate(rest):
            require(stage=='calibration_blank' and r[:2]==g['blank'] and i<64,'blank action')
            records.append(r)
            if r[3]!=0:blank_out='reference_rejected';break
            coins.add(r[2])
            if len(coins)==2:blank_out='blank_checked';break
        require(len(records)==len(rest),'events after blank stopping condition')
        if blank_out=='blank_timeout':require(len(records)==64,'premature blank timeout')
        require(observer['blank']==dict(outcome=blank_out,records=records,endpoint_calls=len(records)),'blank observer')
        outcome='ready' if blank_out=='blank_checked' else blank_out
    else:require(n+used==len(events),'events after calibration timeout')
    require(observer['outcome']==outcome,'combined outcome')
    if outcome=='ready':
        require(p==1 and q==entry['observed_final'] and len(events)-n<=37568,'return/cost')
        require(trial['measured_table_sha256']==sha(canon(columns).encode()),'measured table digest')
    return counts,n,len(events)-n,columns


def check_menu(columns,menu,tickets,weights):
    require(len(columns)==18,'18 columns')
    indexed={(r['input_position'],tuple(r['edge']),r['coin']):(r['decoded_output'],r['current']) for r in columns}
    require(set(indexed)=={(p,e,c) for p in range(3) for e in EDGES for c in (0,1)},'exact unique column coverage')
    inactive=1-sum(weights);P=[[F(inactive*(i==j)) for j in range(3)] for i in range(3)]
    G=[[F(inactive*(i==j)) for j in range(3)] for i in range(3)];primitive=[]
    for t,e in enumerate(tickets):
        for c in (0,1):
            rows=[indexed[p,EDGES[e],c] for p in range(3)]
            for p,(y,j) in enumerate(rows):P[y][p]+=weights[t]/2
            for p in range(3):
                for q in range(3):G[p][q]+=weights[t]/2*(rows[p]==rows[q])
            primitive.append(dict(ticket=t,edge=list(EDGES[e]),coin=c,outputs=[r[0] for r in rows],currents=[r[1] for r in rows]))
    require(G==[[F(i==j) for j in range(3)] for i in range(3)],'fine current isometry')
    width=2*len(tickets)+1;R=[]
    for p in range(3):
        for h in range(width):
            if h==width-1:R.append(p*width+h)
            else:
                y,j=indexed[p,EDGES[tickets[h//2]],h%2]
                R.append(y*width+2*(h//2)+((h%2)^int(j!=0)))
    require(sorted(R)==list(range(3*width)) and all(R[R[i]]==i for i in range(len(R))),'source reversal involution')
    require(menu['tickets']==tickets and list(map(F,menu['ticket_probabilities']))==weights and F(menu['inactive_probability'])==inactive,'declared weights')
    require(menu['P']==[[str(v) for v in r] for r in P] and menu['fine_isometry_gram']==[[str(v) for v in r] for r in G],'measured mean/Gram')
    require(menu['R_folded_permutation']==R and menu['R_dimension']==len(R) and menu['primitive_columns']==primitive,'measured reversal/primitive')


def verify(evidence_path):
    start=time.monotonic();saved=read(evidence_path);inputs=read(HERE/'ADMISSION_INPUTS.json.gz')['rows']
    public=json.loads((HERE/'PUBLIC_GRAPHS.json').read_text());require(len(inputs)==1314,'registered input multiplicity')
    expected=Counter((i,law,seed) for i in range(1314) for law in ('native','always_swap') for seed in (72371001,72371002))
    actual=Counter((t['input_index_audit_only'],t['probe_law_audit_only'],t['calibration_seed_audit_only']) for t in saved['trials'])
    require(actual==expected,'exact 5256 trial coverage')
    total=Counter();prefix=calibration=ready=empty=max_calls=0;representatives={}
    for index,t in enumerate(saved['trials']):
        old=inputs[t['input_index_audit_only']];g=public[f"{old['B']}:{old['M']}"]
        for k in ('B','M','initial_audit_only','total_audit_only','admission_seed_audit_only'):require(t[k]==old[k],'audit input binding')
        counts,n,c,columns=inspect_trial(t,old,g,t['probe_law_audit_only'],t['calibration_seed_audit_only'])
        total.update(counts);prefix+=n;calibration+=c;max_calls=max(max_calls,c)
        if columns is None:empty+=1;continue
        require(t['observer']['outcome']=='ready','registered fair run failed');ready+=1
        key=t['probe_law_audit_only']+':'+str(t['observer']['admission']['sign'])
        if key not in representatives:representatives[key]=dict(table_sha256=t['measured_table_sha256'],columns=columns,source_trial_index=index,actual_blank_records=t['observer']['blank']['records'],observed_sign=t['observer']['admission']['sign'])
        require(columns==representatives[key]['columns'],'background/stream table consistency')
    require(representatives==saved['representative_measured_tables'],'representative evidence binding')
    require(set(saved['compiled_catalogs'])=={r['table_sha256'] for r in representatives.values()},'catalog coverage')
    menus=0
    for r in representatives.values():
        catalogs=saved['compiled_catalogs'][r['table_sha256']];require(set(catalogs)=={'cooling5','weighted_cooling_diagnostic','stars'},'catalog fields')
        check_menu(r['columns'],catalogs['cooling5'],[0,0,2,2,1],[F(1,5)]*5)
        check_menu(r['columns'],catalogs['weighted_cooling_diagnostic'],[0,0,2,2,1],[F(3,20)]*5);menus+=2
        require(len(catalogs['stars'])==132,'star menu multiplicity')
        for item,(m,slot) in zip(catalogs['stars'],((m,s) for m in (4,8,32) for s in range(3*m))):
            require((item['M'],item['slot'])==(m,slot),'star selector coverage')
            v,r0=divmod(slot,m);w=(v+1)%3;tickets=[e for e,edge in enumerate(EDGES) if v in edge]*(m-r0)+[e for e,edge in enumerate(EDGES) if w in edge]*r0
            check_menu(r['columns'],item,tickets,[F(1,2*m)]*(2*m));menus+=1
    require((len(saved['trials']),ready,empty)==(5256,5184,72),'outcome counts')
    require(saved['passed'] is True and (saved['trial_count'],saved['ready_trials'],saved['no_unit_trials'])==(5256,ready,empty),'primary scientific verdict')
    require(dict(total)==saved['literal_trial_calls'] and prefix==saved['literal_admission_prefix_calls'] and calibration==saved['literal_calibration_calls'] and max_calls==saved['maximum_observed_calibration_calls'],'literal counts')
    old=next(r for r in inputs if r['total_audit_only']>0);g=public[f"{old['B']}:{old['M']}"];controls=saved['timeout_controls'];cc=Counter()
    require([t['kind'] for t in controls]==['route','coupon','blank'],'timeout control coverage')
    for t in controls:
        require(t['forced_coin_control']=={'calibration_'+{'coupon':'probe'}.get(t['kind'],t['kind']):0},'forced control definition')
        counts,_,_,_=inspect_trial(t,old,g,'native',72371001,t['forced_coin_control']);cc.update(counts)
        require(t['observer']['outcome']=={'route':'route_timeout','coupon':'coupon_timeout','blank':'blank_timeout'}[t['kind']],'retained timeout')
    require(dict(cc)==saved['literal_timeout_control_calls'],'timeout literal counts')
    require(saved['calibration_bounds']==dict(cap=64,core_calibration_failure_upper=str(F(595,2**64)),with_blank_failure_upper=str(F(597,2**64)),success_probability_lower=str(1-F(597,2**64)),max_calibration_endpoint_calls=37568,expected_endpoint_calls_upper=86,conditional_expected_calls_upper=str(F(86)/(1-F(597,2**64)))),'bound arithmetic')
    return dict(passed=True,trials=5256,ready=ready,no_unit=empty,native_records=sum(total.values()),timeout_records=sum(cc.values()),compiled_menus=menus,producer_imports=False,source_regeneration=False,seconds=time.monotonic()-start)


def main():
    p=argparse.ArgumentParser();p.add_argument('--evidence',type=Path,default=HERE/'records/EVIDENCE.json.gz');p.add_argument('--record-to',type=Path);args=p.parse_args()
    if args.record_to is not None and args.record_to.exists():raise FileExistsError(args.record_to)
    prereg=json.loads((HERE/'PREREG.json').read_text())
    for name,digest in prereg['files'].items():require(sha((HERE/name).read_bytes())==digest,'frozen source/input '+name)
    result=verify(args.evidence)
    if args.record_to is not None:
        with args.record_to.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    elif args.evidence==HERE/'records/EVIDENCE.json.gz':
        old=json.loads((HERE/'INDEPENDENT_CHECK.json').read_text());a=dict(result);b=dict(old);a.pop('seconds');b.pop('seconds');require(a==b,'saved independent result')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
