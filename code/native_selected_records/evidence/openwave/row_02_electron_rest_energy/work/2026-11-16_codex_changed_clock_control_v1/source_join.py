"""One continuing native unit supplies both preparation and changed-clock tables."""
from fractions import Fraction as F
import gzip
import json
from pathlib import Path
import random
import sys

HERE=Path(__file__).resolve().parent
CALIBRATION=HERE.parent/'2026-11-16_codex_record_calibration_v1'
sys.path.insert(0,str(CALIBRATION))
from integration import graph
from admission import admit
from controller import calibrate,check_blank

TRI=((0,1),(1,2),(0,2))
TICKETS=((0,1),(0,1),(0,2),(0,2),(1,2))
STREAMS=(((72365001,72385002,72385003),(72365002,72385005,72385006)),
         ((72365011,72385102,72385103),(72365012,72385105,72385106)))
PROGRESS=[]


def require(ok,message):
    if not ok:raise ValueError(message)


def event(q,edge,coin,law):
    a,b=edge;out=list(q);total=q[a]+q[b]
    out[a]=q[b] if law=='swap' else (total+coin)//2;out[b]=total-out[a]
    return out,[a,b,coin,q[a]-out[a]]


def continue_source(initial,g,seeds,clock_law):
    q=list(initial);events=[];counts={'native':0,'swap':0};allowed={tuple(sorted(e)) for e in g['edges']}
    rng={key:random.Random(seed) for key,seed in zip(('admission','prep','clock'),seeds)}
    def request(edge,stage):
        nonlocal q
        require(tuple(sorted(edge)) in allowed,'physical edge')
        group='prep' if stage.startswith('prep_') else 'clock' if stage.startswith('clock_') else 'admission'
        law=clock_law if stage in ('clock_probe','clock_blank') else 'native'
        q,r=event(q,edge,rng[group].randrange(2),law);counts[law]+=1;events.append([stage,r]);return r
    observer={}
    try:
        entry=admit(request,g);observer['admission']=entry
        require(entry['outcome']=='signed_target_unit','common signed admission')
        for label in ('prep','clock'):
            def wrapped(edge,stage):return request(edge,stage.replace('calibration_',label+'_'))
            observed=calibrate(wrapped,g['target'],entry['sign'],1)
            observer[label]=observed
            if observed['outcome']!='calibrated':return dict(outcome=observed['outcome'],observer=observer,events=events,final=q,counts=counts)
            blank=check_blank(wrapped,g['blank']);observer[label+'_blank']=blank
            if blank['outcome']!='blank_checked':return dict(outcome=blank['outcome'],observer=observer,events=events,final=q,counts=counts)
            require(q==entry['observed_final'],'continuing target/background not restored')
        return dict(outcome='ready',observer=observer,events=events,final=q,counts=counts)
    except Exception as error:return dict(outcome='failed',error=repr(error),observer=observer,events=events,final=q,counts=counts)


def table(rows):
    require(len(rows)==18,'actual18 columns')
    result={(r['input_position'],tuple(r['edge']),r['coin']):r for r in rows}
    require(set(result)=={(p,e,c) for p in range(3) for e in TRI for c in (0,1)},'conditional column coverage')
    return result


def spectral_binding(rows):
    data=table(rows);P=[[F(0) for _ in range(3)] for _ in range(3)];R=[]
    for p in range(3):
        for h in range(11):
            if h==10:R.append(11*p+h);continue
            row=data[p,TICKETS[h//2],h%2];P[row['decoded_output']][p]+=F(1,10)
            R.append(11*row['decoded_output']+(h^int(row['current']!=0)))
    require(sorted(R)==list(range(33)) and all(R[R[i]]==i for i in range(33)),'measured clock reversal')
    means=[]
    for v in ((0,1,-1),(-2,1,1)):
        mv=[sum(P[i][j]*v[j] for j in range(3)) for i in range(3)];j=next(i for i,x in enumerate(v) if x);mu=mv[j]/v[j]
        require(mv==[mu*x for x in v],'actual optical invariant mode');means.append(mu)
    return {'P':[[str(x) for x in row] for row in P],'R_permutation':R,'mode_means':[str(x) for x in means]}


def joint(g,backgrounds,signs,columns,law,ledger):
    n=g['N'];base=backgrounds[0]+backgrounds[1];tables=[table(c) for c in columns]
    for x in range(33):
        for y in range(33):
            initial=list(base)
            for side,column in enumerate((x,y)):
                initial[side*n+g['target'][1]]-=signs[side];initial[side*n+g['target'][column//11]]+=signs[side]
            finals=[]
            for order in ((0,1),(1,0)):
                q=list(initial);outputs=[x,y];records=[]
                for side in order:
                    pos,choice=divmod(outputs[side],11)
                    if choice==10:
                        edge=tuple(side*n+p for p in g['blank']);q,r=event(q,edge,0,law)
                        require(r[3]==0,'actual blank identity')
                    else:
                        edge=TICKETS[choice//2];coin=choice%2;row=tables[side][pos,edge,coin]
                        q,r=event(q,tuple(side*n+g['target'][p] for p in edge),coin,law)
                        require(r[3]==row['current'],'global column differs from actual measured current')
                        outputs[side]=11*row['decoded_output']+(choice^int(r[3]!=0))
                    records.append(r)
                spectator=list(q)
                for side,out in enumerate(outputs):
                    spectator[side*n+g['target'][out//11]]-=signs[side];spectator[side*n+g['target'][1]]+=signs[side]
                require(spectator==base,'joint full background changed');finals.append((q,outputs))
                ledger.append({'input':[x,y],'order':list(order),'output':outputs,'records':records})
            require(finals[0]==finals[1],'joint clock address order')
    require(len(ledger)==2178,'joint coverage')
    return ledger


def run():
    PROGRESS.clear();g=graph(2,32)
    old=json.loads(gzip.decompress((HERE.parent/'2026-11-16_codex_amplified_common_pair_v1/SOURCE_JOIN.json.gz').read_bytes()))
    representatives={};calls={'native':0,'swap':0}
    for clock_law in ('native','swap'):
        for ancestor in old['rows']:
            row={'clock_law':clock_law,'pattern':ancestor['pattern'],'public_stream_pair':ancestor['public_stream_pair'],
                 'signs':ancestor['signs'],'sources':[]};PROGRESS.append(row)
            for side,initial in enumerate(ancestor['initial_states_audit_only']):
                seeds=STREAMS[ancestor['public_stream_pair']][side];result=continue_source(initial,g,seeds,clock_law)
                row['sources'].append({'initial_audit_only':initial,'public_seeds':seeds,**result})
                require(result['outcome']=='ready','continuing double calibration failed')
                sign=result['observer']['admission']['sign'];require(sign==ancestor['signs'][side],'actual sign')
                require(result['final']==ancestor['source_audits'][side]['final'],'same actual source background')
                for key,value in result['counts'].items():calls[key]+=value
                entry={'prep_columns':result['observer']['prep']['columns'],'clock_columns':result['observer']['clock']['columns'],
                       'prep_blank':result['observer']['prep_blank']['records'],'clock_blank':result['observer']['clock_blank']['records'],
                       'spectrum':spectral_binding(result['observer']['clock']['columns'])}
                key=f'{clock_law}:{sign}'
                if key in representatives:
                    require(entry['prep_columns']==representatives[key]['prep_columns'] and entry['clock_columns']==representatives[key]['clock_columns'],'background changes measured maps')
                else:representatives[key]=entry
            row['joint_columns']=[]
            joint(g,[r['final'] for r in row['sources']],row['signs'],[r['observer']['clock']['columns'] for r in row['sources']],clock_law,row['joint_columns'])
    require(len(PROGRESS)==32 and len(representatives)==4,'full changed-clock source menu')
    return {'passed':True,'graph':g,'pairs':PROGRESS,'representatives':representatives,'actual_continuing_sources':64,
            'literal_continuing_calls':calls,'literal_global_column_evaluations':32*2178*2,
            'per_source_timeout_upper':str(F(1816,2**64)),'independent_pair_entry_lower':str((1-F(1816,2**64))**2),
            'per_source_entry_attempt_max':38360+2*37568,
            'joint_graph_vertices':2*g['N'],'idle_connecting_edge':[g['blank'][1],g['N']+g['blank'][0]],
            'idle_connecting_edge_events':0,'quantum_scope':'Exact full-code joint operator binding, not sampled quantum histories.'}
