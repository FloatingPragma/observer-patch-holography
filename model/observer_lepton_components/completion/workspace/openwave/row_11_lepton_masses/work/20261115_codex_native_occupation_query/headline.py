#!/usr/bin/env python3
"""Exact native-current coherent occupation query; no physical mass claim."""
from pathlib import Path
from fractions import Fraction as F
from dataclasses import asdict,replace
from itertools import product
from collections import Counter
import sys,json,hashlib,argparse,time
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
sys.path.insert(0,str(ROOT/'reverse-engineering-reality'))
from model.observer import native
from model.native_electric_gauss import source as g
from model.defect_field_bridge import full_domain as car
SEAM=(0,10); TRI=((0,5),(0,11),(5,11)); LIMIT=1
PIN='dd7c87bcda73f55990690effd0633020ff05b9131205d2314793329f8e0c811d'
CALLS=Counter(); CAR_CALLS=0

def canonical(v):return json.dumps(v,sort_keys=True,indent=2)+'\n'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def loads(s):return tuple(int(i in s.pair) for i in range(12))
def enc(s,fl):
    other=next(i for i in s.pair if i!=s.tagged_port)
    return (1<<(6*s.tagged_port+fl[0]))|(1<<(6*other+fl[1]))
def statekey(s):return (s.tagged_port,s.memory,s.flux,s.pair,s.charge_number)
def eps(s):return 1 if s.tagged_port<next(i for i in s.pair if i!=s.tagged_port) else -1

def native_column(s,fl,coin,kind):
    global CAR_CALLS
    out,(_,j)=native.native(loads(s),SEAM,coin);CALLS[kind]+=1
    hop=g.apply_native_current(s,SEAM,abs(j),LIMIT)
    rows={}
    for law in ('C','E'):
        raw,mask,sign,jc=car.branch_action(enc(s,fl),law,coin,sites=12,flavors=6,seam=SEAM)
        CAR_CALLS+=1
        assert j==jc and mask.bit_count()==2
        rows[law]={'raw':raw,'proposal_mask':str(mask),'proposal_sign':sign,'current':jc}
        if not hop.overflow:
            assert out==loads(hop.state) and mask==enc(hop.state,fl)
            assert eps(s)*sign==eps(hop.state)
    assert g.gauss(hop.state)==g.gauss(s)
    assert rows['C']['proposal_sign']==rows['E']['proposal_sign']
    return hop,rows,{'coin':coin,'native_output':list(out),'current':j,
                     'proposed_signed_current':hop.proposed_signed_current,
                     'realized_signed_current':hop.realized_signed_current,
                     'overflow':hop.overflow,'records':rows}

def input_states():
    states=[]
    for pos,m,abc in product((0,5,11,7),range(6),product((-1,0,1),repeat=3)):
        flux=[0]*30;flux[g.EDGE_INDEX[(0,1)]]=1
        if pos==7:flux[g.EDGE_INDEX[(0,7)]]=-1
        for e,v in zip(TRI,abc):flux[g.EDGE_INDEX[e]]=v
        s=g.ElectricState(tuple(sorted((pos,1))),pos,1,m,tuple(flux))
        if not any(g.gauss(s)):states.append(s)
    return sorted(states,key=statekey)

def transfer(q,f,p):
    if q:
        for k in range(3):
            if (f,p)==(k,3):return q,0,k
            if (f,p)==(0,k):return q,k,3
    return q,f,p

def run():
    global CAR_CALLS
    CALLS.clear();CAR_CALLS=0;checks={}
    def check(k,v):
        checks[k]=bool(v)
        if not v:raise AssertionError(k)
    pins=json.loads((HERE/'SOURCES.json').read_text())
    check('all_frozen_source_bindings',all(sha(ROOT/p)==h for p,h in pins.items()))
    check('literal_native_pin',sha(Path(native.__file__))==PIN)
    A=input_states();check('full60_gauss_input_bases',len(A)==60)
    check('public_auxiliary_empty_and_zero_link',all(10 not in s.pair and s.flux[g.EDGE_INDEX[SEAM]]==0 for s in A))
    forward=[];inverse=[];images=set();temporary=[]
    for i,s in enumerate(A):
        for fl in product(range(6),repeat=2):
            hop,raw,rec=native_column(s,fl,0,'forward')
            check_n=int(s.tagged_port==0)
            assert rec['current']==check_n and not hop.overflow
            assert all((r['raw']%4 in (1,2))==bool(check_n) for r in raw.values())
            temporary.append((i,s,fl,hop,raw,rec))
            if check_n:
                images.add(hop.state)
                back,braw,brec=native_column(hop.state,fl,1,'reverse_image')
                assert back.state==s and not back.overflow
                assert all(braw[l]['proposal_sign']*raw[l]['proposal_sign']==1 for l in ('C','E'))
                inverse.append({'input':i,'flavors':fl,'moved_state':asdict(hop.state),**brec})
    states=A+sorted(images-set(A),key=statekey);index={s:i for i,s in enumerate(states)}
    check('full78_working_base_chart',len(states)==78)
    tables={law:{} for law in ('C','E')};transpositions={law:[] for law in ('C','E')}
    for i,s,fl,hop,raw,rec in temporary:
        j=index[hop.state]
        forward.append({'input':i,'output':j,'flavors':fl,**rec})
        if rec['current']:
            for law in ('C','E'):
                r=raw[law]['raw'];sign=raw[law]['proposal_sign']
                left=(i,fl,0);right=(j,fl,r)
                assert left not in tables[law] and right not in tables[law]
                tables[law][left]=(j,r,sign);tables[law][right]=(i,0,sign)
                transpositions[law].append({'input':i,'image':j,'flavors':fl,'raw':r,'sign':sign})
    check('all2160_native_forward_rows',len(forward)==2160)
    check('all648_native_reverse_rows',len(inverse)==648)
    def Q(law,key,fl,hold=False):
        i,r,b=key
        if hold:return key,1
        j,t,sign=tables[law].get((i,fl,r),(i,r,1))
        return (j,t,b),sign
    def copy_current(key,enabled=True):
        i,r,b=key;flag=int(r%4 in (1,2)) if enabled else 0
        return i,r,b^flag
    full_counts={};query_rows={}
    for law in ('C','E'):
        full=0;queries=[]
        for fl in product(range(6),repeat=2):
            outputs=set()
            for i,r,b in product(range(len(states)),range(16),(0,1)):
                k,sg=Q(law,(i,r,b),fl);back,sb=Q(law,k,fl)
                assert back==(i,r,b) and sg*sb==1 and k not in outputs
                outputs.add(k);full+=1
            assert len(outputs)==len(states)*32
            for i,b in product(range(len(A)),(0,1)):
                k,s1=Q(law,(i,0,b),fl);copied=copy_current(k);final,s2=Q(law,copied,fl)
                want=(i,0,b^int(A[i].tagged_port==0))
                assert final==want and s1*s2==1
                queries.append({'source':i,'flavors':fl,'input_query':b,'after_first_Q':k,
                                'after_record_copy':copied,'output':final,'phase':s1*s2})
        full_counts[law]=full;query_rows[law]=queries
    check('complete_signed_involutions',all(n==89856 for n in full_counts.values()))
    check('all8640_coherent_query_columns',sum(len(v) for v in query_rows.values())==8640)
    # Exact sparse coherent tests. Integer amplitudes have global norm squared2.
    S=g.pair_string((0,1),1,0);R=g.pair_string((7,0,1),1,0)
    si,ri=index[S],index[R];fl=(0,0);dress=eps(S)*eps(R)
    initial={(si,0,0):dress,(ri,0,0):1}
    def move(vec,operation):
        ans={}
        for k,a in vec.items():
            k2,s=operation(k);ans[k2]=ans.get(k2,0)+a*s
        return {k:a for k,a in ans.items() if a}
    def query(vec,law,hold=False,copy=True):
        vec=move(vec,lambda k:Q(law,k,fl,hold))
        vec=move(vec,lambda k:(copy_current(k,copy),1))
        return move(vec,lambda k:Q(law,k,fl,hold))
    def phase(vec):return move(vec,lambda k:(k,(-1)**k[2]))
    def pointer(vec):
        mass=F(sum(a*a for a in vec.values()),2);contrast=F(0)
        for (i,r,b),a in vec.items():
            if i in (si,ri):
                partner=(ri if i==si else si,r,b)
                contrast+=F(dress*a*vec.get(partner,0),2)
        return {'weight':str(mass),'contrast':str(contrast),
                'plus':str((mass+contrast)/2),'minus':str((mass-contrast)/2),
                'source_query_state':[{'key':k,'amplitude_integer':a} for k,a in sorted(vec.items())]}
    interference={}
    for law in ('C','E'):
        cases={}
        for name,hold,copy in [('native',False,True),('HOLD',True,True),('missed_copy',False,False)]:
            vec=query(initial,law,hold,copy);vec=phase(vec);vec=query(vec,law,hold,copy)
            cases[name]=[pointer(vec)]
        first=move(initial,lambda k:Q(law,k,fl))
        branches={r:{k:a for k,a in first.items() if k[1]==r} for r in sorted({k[1] for k in first})}
        committed=[]
        for r,vec in branches.items():
            vec=move(vec,lambda k:(copy_current(k),1));vec=move(vec,lambda k:Q(law,k,fl))
            vec=phase(vec);vec=query(vec,law)
            committed.append({'committed_raw':r,**pointer(vec)})
        cases['committed_raw']=committed
        minus={(si,0,0):-dress,(ri,0,0):1}
        cases['minus_plant']=[pointer(query(phase(query(minus,law)),law))]
        sums={k:sum(F(x['plus']) for x in rows) for k,rows in cases.items()}
        assert sums=={'native':F(0),'HOLD':F(1),'missed_copy':F(1),'committed_raw':F(1,2),'minus_plant':F(1)}
        assert all(sum(F(x['plus'])+F(x['minus']) for x in rows)==1 for rows in cases.values())
        interference[law]={'cases':cases,'plus_probabilities':{k:str(v) for k,v in sums.items()}}
    check('native_phase_unquery_and_frozen_controls',True)
    # Admission belongs to the outer instrument, never to intermediate Q steps.
    def admitted(s):return s in index and index[s]<len(A)
    witnesses=[]
    outside=[('auxiliary_positive',g.pair_string((10,0,1),1,0)),
             ('both_query_endpoints_occupied',g.pair_string((0,10),1,0))]
    flux=list(S.flux)
    for e,v in [((0,10),-1),((10,11),-1),((0,11),1)]:flux[g.EDGE_INDEX[e]]+=v
    outside.append(('outward_auxiliary_boundary',replace(S,flux=tuple(flux))))
    for name,s in outside:
        assert not admitted(s) and not any(g.gauss(s))
        hop,raw,rec=native_column(s,(0,0),0,'outside_witness')
        witnesses.append({'name':name,'input_state':asdict(s),'hypothetical_native':rec,
                          'actual_query_outcome':'invalid_domain_identity','actual_output':asdict(s),
                          'query_bit_unchanged':True})
        if name=='outward_auxiliary_boundary':
            assert hop.overflow and hop.state==s and rec['current']==1 and hop.realized_signed_current==0
        if name=='both_query_endpoints_occupied':assert rec['current']==0 and 0 in s.pair
    collision=(1<<0)|(1<<1)
    for law in ('C','E'):
        raw,out,sg,j=car.branch_action(collision,law,0,sites=12,flavors=6,seam=SEAM);CAR_CALLS+=1
        assert raw==car.COLLISION_RAW and out==collision and sg==1 and j==0
        witnesses.append({'name':'nonbinary_full_Fock_'+law,'mask':str(collision),'native_CAR_record':raw,
                          'actual_query_outcome':'invalid_domain_identity','actual_output_mask':str(collision)})
    check('auxiliary_boundary_and_Fock_complements_retained',len(witnesses)==5)
    datamap=[]
    for q,f,p in product((0,1),range(3),range(4)):
        out=transfer(q,f,p);assert transfer(*out)==(q,f,p)
        if p==3:assert out==((1,0,f) if q else(q,f,3))
        datamap.append({'input':[q,f,p],'output':out,'sign':1})
    check('neutral_family_transfer_full24_permutation',len({tuple(x['output']) for x in datamap})==24)
    check('spin_and_Gauss_factorization',all(r['records']['C']['proposal_sign']==r['records']['E']['proposal_sign'] for r in forward))
    return {'schema':'oph.row11.native_occupation_query.v1','candidate':'A-native-occupation-query','class':'N+',
            'checks':checks,'all_pass':all(checks.values()),'actual_native_calls':dict(CALLS),
            'total_actual_native_calls':sum(CALLS.values()),'coherent_C_E_columns':CAR_CALLS,
            'source_bindings':pins,'source_input_count':len(A),'source_working_count':len(states),
            'source_states':[asdict(s) for s in states],'source_valid_inputs':list(range(len(A))),
            'forward_records':forward,'reverse_image_records':inverse,
            'completion_transpositions':transpositions,'completion_identity_extension':
            'Each listed signed transposition acts in both directions; every unlisted source/flavor/raw/query address is identity. Query bit is spectator for Q.',
            'full_completion_columns_checked':full_counts,'query_columns':query_rows,
            'interference':interference,'complement_witnesses':witnesses,'neutral_data_transfer':datamap,
            'all_physical_Gauss_and_spin_preserved':True,'new_M6_pulse_execution':False,
            'combined_vacuum_filter_executed':False,'physical_row_established':False,
            'compression':'Complete Q rows are losslessly represented by native transpositions plus identity; conjugated-query rows are explicitly retained. Ramsey tests reuse exact executed source-column tables, with no new endpoint calls counted.',
            'coherent_source_scope':'On admitted blank-current inputs, Q equals the fixed-coin native C/E isometry with its actual CAR/Gauss image. The signed involutive extension and inverse are explicit coherent services; classical bytes alone do not select them. Audited column records are not simultaneous classical copies of the uncomputed quantum archive.',
            'scope':'Finite native-current coherent occupation query on the public empty/zero-link auxiliary domain; conditional neutral-data embedding only. No physical electric unit, clock, energy or lepton minimum.'}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',type=Path);ap.add_argument('--check',type=Path);args=ap.parse_args()
    started=time.perf_counter();r=run();blob=canonical(r)
    if args.write:
        with args.write.open('x') as f:f.write(blob)
    if args.check and args.check.read_text()!=blob:raise SystemExit('Byte replay mismatch')
    print(canonical({'all_pass':r['all_pass'],'checks':len(r['checks']),'native_calls':r['total_actual_native_calls'],
      'coherent_C_E_columns':r['coherent_C_E_columns'],'completion_columns':r['full_completion_columns_checked'],
      'source_inputs':r['source_input_count'],'working_source_states':r['source_working_count'],
      'plus_probabilities':{k:v['plus_probabilities'] for k,v in r['interference'].items()},
      'seconds':time.perf_counter()-started,'physical_row_established':False}))
