#!/usr/bin/env python3
"""Exact retained N=1 source-family query, without measured exports or M6 claims."""
from pathlib import Path
from fractions import Fraction as F
from dataclasses import dataclass, asdict, replace
from itertools import product
from collections import Counter
import sys, json, hashlib, argparse, time
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
sys.path.insert(0,str(ROOT/'reverse-engineering-reality'))
import sympy as sy
from model.observer import native
from model.native_electric_gauss import source as g
from model.defect_field_bridge import full_domain as car
from model.finite_family_clock import primitives as fp, family
from model.lepton_register_round import carrier, walk as wk

PIN='dd7c87bcda73f55990690effd0633020ff05b9131205d2314793329f8e0c811d'
SEAM=(7,8); TRI=((0,5),(5,11),(0,11)); LIMIT=1
CALLS=Counter(); CAR_CALLS=0

@dataclass(frozen=True,slots=True)
class K:
    """a+b sqrt(2)+i c+i d sqrt(2), all four coefficients exact rational."""
    v: tuple
    def __add__(self,o):
        if not isinstance(o,K):o=k(o)
        return K(tuple(a+b for a,b in zip(self.v,o.v)))
    __radd__=__add__
    def __neg__(self):return K(tuple(-a for a in self.v))
    def __sub__(self,o):return self+-k(o)
    def __rsub__(self,o):return k(o)+-self
    def __mul__(self,o):
        o=k(o);r=[F()]*4
        for i,a in enumerate(self.v):
            if a:
                for j,b in enumerate(o.v):
                    if b:r[i^j]+=a*b*(2 if i&j&1 else 1)*(-1 if i&j&2 else 1)
        return K(tuple(r))
    __rmul__=__mul__
    def __truediv__(self,o):return self*k(F(1,o))
    def __bool__(self):return any(self.v)
    def conj(self):return K((self.v[0],self.v[1],-self.v[2],-self.v[3]))
    def data(self):return [str(x) for x in self.v]
def k(a=0,b=0,c=0,d=0):
    if isinstance(a,K):return a
    return K(tuple(map(F,(a,b,c,d))))
ZERO=k();ONE=k(1);SQ2=k(0,1);Z=k(F(3,5),0,F(4,5));ZINV=Z.conj()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def canonical(x):return json.dumps(x,indent=2,sort_keys=True)+'\n'
def add(v,i,x):
    v[i]=v.get(i,ZERO)+x
    if not v[i]:del v[i]
def dot(a,b):return sum((x.conj()*b.get(i,ZERO) for i,x in a.items()),ZERO)
def scaled(v,a):return {i:x*a for i,x in v.items() if x*a}
def plus(a,b,scale=ONE):
    c=dict(a)
    for i,x in b.items():add(c,i,x*scale)
    return c
def vecdata(v):return [[i,x.data()] for i,x in sorted(v.items())]
def matdata(m):return [[x.data() for x in row] for row in m]
def index(ic,ie,p):return (12*ic+ie)*4+p
def unpack(i):
    obj,p=divmod(i,4);ic,ie=divmod(obj,12);return ic,ie,p
def step(ic):return (1 if ic%2==0 else -1) if (ic//2)%2 else 0

def native_col(s,fl,coin,kind):
    global CAR_CALLS
    loads=tuple(int(i in s.pair) for i in range(12))
    out,(_,j)=native.native(loads,SEAM,coin);CALLS[kind]+=1
    h=g.apply_native_current(s,SEAM,abs(j),LIMIT)
    other=next(i for i in s.pair if i!=s.tagged_port)
    enc=lambda t:(1<<(6*t.tagged_port+fl[0]))|(1<<(6*next(i for i in t.pair if i!=t.tagged_port)+fl[1]))
    rows={}
    for law in ('C','E'):
        raw,mask,sign,jc=car.branch_action(enc(s),law,coin,sites=12,flavors=6,seam=SEAM);CAR_CALLS+=1
        assert jc==j
        if not h.overflow:
            assert mask==enc(h.state) and out==tuple(int(i in h.state.pair) for i in range(12))
            assert (1 if s.tagged_port<other else -1)*sign==(1 if h.state.tagged_port<other else -1)
        rows[law]={'raw':raw,'CAR_sign':sign,'mask':str(mask)}
    assert g.gauss(s)==g.gauss(h.state)
    return h,rows,{'coin':coin,'native_current':j,'overflow':h.overflow,
        'realized_signed_current':h.realized_signed_current,'C_E':rows}

def arm_query(check):
    states=[]
    for pos,m,abc in product((0,5,11,7),range(6),product((-1,0,1),repeat=3)):
        flux=[0]*30;flux[g.EDGE_INDEX[(0,1)]]=1
        if pos==7:flux[g.EDGE_INDEX[(0,7)]]=-1
        for e,v in zip(TRI,abc):flux[g.EDGE_INDEX[e]]=v
        s=g.ElectricState(tuple(sorted((pos,1))),pos,1,m,tuple(flux))
        if not any(g.gauss(s)):states.append(s)
    states.sort(key=lambda s:(s.tagged_port,s.memory,s.flux))
    check('all60_source_bases',len(states)==60)
    check('native_auxiliary_edge_and_empty_zero_link',SEAM in g.EDGE_INDEX and all(8 not in s.pair and all(s.flux[j]==0 for e,j in g.EDGE_INDEX.items() if 8 in e) for s in states))
    forward=[];backward=[];moving=[];images=set()
    for i,s in enumerate(states):
        for fl in product(range(6),repeat=2):
            h,rows,record=native_col(s,fl,0,'arm_forward')
            assert record['native_current']==int(s.tagged_port==7) and not h.overflow
            assert all((r['raw']%4 in (1,2))==bool(s.tagged_port==7) for r in rows.values())
            forward.append({'input':i,'flavors':fl,'image':asdict(h.state),**record})
            if s.tagged_port==7:
                back,rr,br=native_col(h.state,fl,1,'arm_reverse')
                assert back.state==s and all(rr[l]['CAR_sign']*rows[l]['CAR_sign']==1 for l in rows)
                backward.append({'input':i,'flavors':fl,**br})
                moving.append((i,fl,h.state,rows));images.add(h.state)
    work=states+sorted(images-set(states),key=lambda s:(s.tagged_port,s.memory,s.flux));idx={s:i for i,s in enumerate(work)}
    tables={l:{} for l in ('C','E')};pairs={l:[] for l in tables}
    for i,fl,t,rows in moving:
        for l,r in rows.items():
            a=(i,fl,0);b=(idx[t],fl,r['raw']);sg=r['CAR_sign']
            assert a not in tables[l] and b not in tables[l]
            tables[l][a]=(idx[t],r['raw'],sg);tables[l][b]=(i,0,sg)
            pairs[l].append([i,list(fl),idx[t],r['raw'],sg])
    columns={};query_rows={};controls={}
    for law in tables:
        def Q(key,fl):
            i,r,b=key;j,t,sg=tables[law].get((i,fl,r),(i,r,1));return(j,t,b),sg
        n=0;qr=[];cr=[]
        for fl in product(range(6),repeat=2):
            image=set()
            for i,r,b in product(range(len(work)),range(16),(0,1)):
                mid,sg=Q((i,r,b),fl);end,sg2=Q(mid,fl)
                assert end==(i,r,b) and sg*sg2==1 and mid not in image
                image.add(mid);n+=1
            for i,b in product(range(len(states)),(0,1)):
                mid,sg=Q((i,0,b),fl)
                copied=(mid[0],mid[1],mid[2]^int(mid[1]%4 in (1,2)))
                out,sg2=Q(copied,fl);out=(out[0],out[1],out[2]^1)
                assert out==(i,0,b^int(states[i].tagged_port!=7)) and sg*sg2==1
                qr.append([i,list(fl),b,list(mid),list(copied),list(out),sg*sg2])
                # Replace this charged source completion by whole-event HOLD.
                # Its raw record remains0, the frozen copy gives0, and the
                # final X makes chi1 on BOTH source and reference arms.
                held=(i,0,b);hcopied=(held[0],held[1],held[2]^int(held[1]%4 in (1,2)))
                hout=(hcopied[0],hcopied[1],hcopied[2]^1)
                assert hout==(i,0,b^1)
                cr.append([i,list(fl),b,list(hout),1])
        columns[law]=n;query_rows[law]=qr;controls[law]=cr
    check('all179712_signed_completion_columns',sum(columns.values())==179712)
    check('all8640_whole_arm_query_columns',sum(map(len,query_rows.values()))==8640)
    # Complete identity-invalid instrument, with real source boundary witnesses.
    ref=g.pair_string((7,0,1),1,0)
    flux=list(ref.flux)
    for e,v in (((7,8),-1),((6,8),1),((6,7),-1)):flux[g.EDGE_INDEX[e]]+=v
    outside=[('positive_at_auxiliary',g.pair_string((8,7,0,1),1,0)),
        ('both_query_endpoints_occupied',g.pair_string((7,8),1,0)),
        ('outward_auxiliary_boundary',replace(ref,flux=tuple(flux)))]
    witnesses=[]
    for name,s in outside:
        assert s not in states and not any(g.gauss(s))
        h,rows,r=native_col(s,(0,0),0,'arm_outside_witness')
        if name=='outward_auxiliary_boundary':assert h.overflow and h.state==s and r['native_current']==1 and h.realized_signed_current==0
        if name=='both_query_endpoints_occupied':assert r['native_current']==0
        witnesses.append({'name':name,'source':asdict(s),'hypothetical_native':r,'actual_query':'invalid_identity'})
    for law in ('C','E'):
        mask=(1<<(7*6))|(1<<(7*6+1));raw,out,sg,j=car.branch_action(mask,law,0,sites=12,flavors=6,seam=SEAM)
        assert raw==car.COLLISION_RAW and out==mask and sg==1 and j==0
        witnesses.append({'name':'nonbinary_Fock_'+law,'raw':raw,'actual_query':'invalid_identity'})
    check('all_query_complements_retained',len(witnesses)==5)
    transfer=[]
    def W(q,f,p):
        if q:
            for j in range(3):
                if(f,p)==(j,3):return q,0,j
                if(f,p)==(0,j):return q,j,3
        return q,f,p
    for q,f,p in product((0,1),range(3),range(4)):
        out=W(q,f,p);assert W(*out)==(q,f,p)
        transfer.append({'input':[q,f,p],'output':out})
    check('full24_equal_charge_transfer_columns',len({tuple(x['output']) for x in transfer})==24)
    return {'input_states':[asdict(s) for s in states],'working_states':[asdict(s) for s in work],
        'forward':forward,'reverse':backward,'signed_pairs':pairs,'full_completion_columns':columns,
        'query_rows':query_rows,'HOLD_query_rows':controls,'invalid_witnesses':witnesses,'neutral_transfer':transfer}

def count_source(check):
    w=carrier.carrier_walk();cfg=wk.config(2,2);records=[];calls=[]
    for cell in ('X','Y'):
        rows=[]
        for j,(a,b) in enumerate(w.addresses):
            im=wk.images(cfg,0,a,b);assert len(im)==1 and im[0][2]==1
            t,r,amp=im[0];assert wk.images(t,r,a,b)==((cfg,0,F(1)),)
            rows.append(r)
            for coin in (0,1):
                out,(_,current)=native.native(cfg,(a,b),coin);CALLS['count_native']+=1
                assert out==t and current==(0 if r==0 else 1 if r==2 else -1)
                calls.append({'cell':cell,'address':j,'coin':coin,'seam':[a,b],'raw':r,'current':current,'output':out})
        records.append(rows)
    flags=([int(r==2) for r in records[0]],[int(r in (2,3)) for r in records[1]])
    R=[[F(1,30)-int(i==j) for j in range(60)] for i in range(60)]
    projector_rows=[]
    for side,flag in enumerate(flags):
        mat=[]
        for col in range(60):
            a=w.step({(cfg,0,col):F(1)})
            a=wk.restrict(a,lambda key:(key[1]==2 if side==0 else key[1] in (2,3)))
            out=w.inverse_step(a)
            assert all(c==cfg and rec==0 for c,rec,j in out)
            v=[out.get((cfg,0,j),F()) for j in range(60)]
            want=[sum((R[i][l]*flag[l]*R[l][col] for l in range(60)),F()) for i in range(60)]
            assert v==want;mat.append(v)
        projector_rows.append([[str(x) for x in row] for row in mat])
    check('all120_actual_reducing_source_columns',True)
    check('native240_count_calls',len(calls)==240)
    # Full local24 completion and all60 address columns prove the full source
    # complement, not only prepared fibres. Matrix products are exact fractions.
    C=wk.completion()
    assert all(sum(C[l][i]*C[l][j] for l in range(24))==int(i==j) for i in range(24) for j in range(24))
    assert all(sum(R[l][i]*R[l][j] for l in range(60))==int(i==j) for i in range(60) for j in range(60))
    check('complete_local24_and_R60_unitaries',True)
    rows=[];counts=Counter();vectors={}
    for kind in ('xor','and'):
        q1=[]
        for ix,iy in product(range(60),repeat=2):
            fx,fy=flags[0][ix],flags[1][iy]
            x=F(sum(flags[0]),30)-fx;y=F(sum(flags[1]),30)-fy
            h=x*y if kind=='and' else x+y-2*x*y
            q1.append(h)
            if kind=='xor':
                rows.append([ix,iy,records[0][ix],records[1][iy],fx,fy,fx^fy,fx*fy]);counts[fx^fy]+=1
        p=F(1,72) if kind=='and' else F(2,9)
        assert sum(q1,F())/3600==p and sum((x*x for x in q1),F())/3600==p
        assert sum((x*(1-x) for x in q1),F())==0
        vectors[kind]={'q1_numerator':[str(x) for x in q1],
            'q0_numerator':[str(1-x) for x in q1],'common_denominator':60,
            'q_weights':[str(1-p),str(p)]}
    check('all3600_labels_and_orthogonal_retained_count_vectors',counts[1]==800 and len(rows)==3600)
    return {'addresses':w.addresses,'original_native_records':calls,'joint_labels':rows,
        'count_flags':flags,'one_cell_projector_columns':projector_rows,
        'local24_completion':[[str(x) for x in row] for row in C],
        'R60':[[str(x) for x in row] for row in R], 'retained_count_vectors':vectors,
        'full_count_source_dimension_per_cell':78*4*60,
        'restricted_fibre_dimension_per_cell':60,
        'full_address_projector_definition':'E_kl=(R60 tensor R60)|kl><kl|(R60 tensor R60); sum all3600 E_kl=I',
        'full_source_complement_definition':'Use complete Walk B=Qrepair R60 on all78 two-load configurations,4 records,60 addresses; Qcount=Bpair† Pi_XOR Bpair. It reduces the blank-double fibre; its orthogonal complement is retained by this same formula.'}

def family_query(check):
    coherent=family.native_event_check()
    check('actual360_coherent_native_columns',coherent['actual_Observer_repair_columns']==360 and coherent['exact_mismatches']==0)
    B=fp.oriented_record()*fp.repair()*fp.address_reflection()
    assert B.T*B==sy.eye(12)
    columns=lambda m:[[(i,k(F(m[i,j]))) for i in range(12) if m[i,j]] for j in range(12)]
    bcols=columns(B);bicols=columns(B.T)
    phi={}
    for obj,v in enumerate(fp.analyzer()):
        if v:
            x=sy.simplify(v*sy.sqrt(24));assert x in (1,sy.sqrt(2));phi[obj]=ONE if x==1 else SQ2
    po=fp.state(sy.sqrt(sy.Rational(1,3)),0,sy.sqrt(sy.Rational(2,3)))
    prep={}
    for obj,v in enumerate(sy.kronecker_product(po,po)):
        if v:
            x=sy.simplify(v*6);assert x in (1,sy.sqrt(2),2);prep[obj]=SQ2 if x==sy.sqrt(2) else k(int(x))
    assert sum((x*x for x in phi.values()),ZERO)==k(24)
    assert sum((x*x for x in prep.values()),ZERO)==k(36)
    def factor(vec,which,inverse=False):
        out={};tab=bicols if inverse else bcols
        for i,a in vec.items():
            ic,ie,p=unpack(i)
            for t,co in tab[ic if which==0 else ie]:add(out,index(t if which==0 else ic,t if which else ie,p),a*co)
        return out
    def shift(vec,inverse=False):
        out={}
        for i,a in vec.items():
            ic,ie,p=unpack(i);p2=(p+(-1 if inverse else 1)*step(ic))%3 if p<3 else p
            add(out,index(ic,ie,p2),a)
        return out
    def phase(vec,q,inverse=False):
        if q==0:return dict(vec)
        out={}
        for i,a in vec.items():
            st=step(unpack(i)[0])*(-1 if inverse else 1)
            out[i]=a*(Z if st==1 else ZINV if st==-1 else ONE)
        return out
    def tick(vec,q,inverse=False):
        if inverse:
            return factor(factor(shift(phase(vec,q,True),True),1,True),0,True)
        return phase(shift(factor(factor(vec,0),1)),q)
    def pi(vec):
        coefficients=[ZERO]*4
        for i,a in vec.items():
            ic,ie,p=unpack(i);coefficients[p]+=phi.get(12*ic+ie,ZERO)*a
        out={}
        for obj,a in phi.items():
            for p,c in enumerate(coefficients):
                if c:out[4*obj+p]=a*c/24
        return out
    def project(vec,q):return tick(pi(tick(vec,q)),q,True)
    # Inverse is the same factors reversed, not an adjoint of a reduced map.
    inverse_columns=0
    for q,i in product((0,1),range(576)):
        v={i:ONE};assert tick(tick(v,q),q,True)==v;inverse_columns+=1
    check('all1152_full_family_source_inverse_columns',inverse_columns==1152)
    preps=[{4*obj+p:a for obj,a in prep.items()} for p in range(4)]
    phis=[{4*obj+p:a for obj,a in phi.items()} for p in range(4)]
    factors={};Ys={};As={};nulls={}
    for q in (0,1):
        vs=[tick(v,q,True) for v in phis]
        assert all(dot(vs[i],vs[j])==k(24*int(i==j)) for i in range(4) for j in range(4))
        # This checks the source-word projector on all576 columns against its
        # retained rank4 factor, rather than applying a compressed effect.
        for i in range(576):
            basis={i:ONE};direct=project(basis,q);rank={}
            for v in vs:
                co=dot(v,basis)/24
                if co:rank=plus(rank,scaled(v,co))
            assert direct==rank
        ys=[project(v,q) for v in preps]
        amat=[[dot(preps[i],ys[j])/36 for j in range(4)] for i in range(4)]
        assert all(dot(ys[i],ys[j])/36==amat[i][j] for i in range(4) for j in range(4))
        assert all(project(v,q)==v for v in vs)
        zero_obj=next(obj for obj in range(144) if obj not in phi)
        null=tick({4*zero_obj:ONE},q,True);assert project(null,q)=={}
        factors[str(q)]=[vecdata(v) for v in vs];Ys[q]=ys;As[q]=amat;nulls[str(q)]=vecdata(null)
    check('all1152_full_projector_columns_and_both_query_complements',True)
    check('four_source_rank_vectors_orthogonal_per_count_block',True)
    check('analyzer_range_and_null_plants_detected',True)
    def expected_A(q):
        z=Z if q else ONE;zb=z.conj()
        cc=SQ2*z/27+zb*zb/54
        ans=[[ZERO for _ in range(3)] for _ in range(3)]
        for j in range(3):ans[j][j]=k(F(2,27));ans[(j+1)%3][j]=cc;ans[(j-1)%3][j]=cc.conj()
        return ans
    for q in (0,1):assert [row[:3] for row in As[q][:3]]==expected_A(q)
    mixed={}
    for name,p in (('native_xor',F(2,9)),('changed_AND',F(1,72)),('count_HOLD',F(0))):
        a=[[(1-p)*As[0][i][j]+p*As[1][i][j] for j in range(3)] for i in range(3)]
        assert sum((a[i][i] for i in range(3)),ZERO)==k(F(2,9))
        # Frozen normalized family ray (1,1,1)/sqrt3.
        prob=sum((x for row in a for x in row),ZERO)/3
        mixed[name]={'compression':matdata(a),'query0_invariant_family_ray':prob.data(),
            'query1_invariant_family_ray':(ONE-prob).data()}
    check('full_source_contraction_equals_A1',True)
    check('frozen_AND_and_count_HOLD_controls_detected',mixed['native_xor']['query0_invariant_family_ray']!=mixed['changed_AND']['query0_invariant_family_ray'] and mixed['native_xor']['query0_invariant_family_ray']!=mixed['count_HOLD']['query0_invariant_family_ray'])
    # Complete original query effects (q0=P, q1=I-P), all internal outputs retained.
    for q in (0,1):
        for i,j in product(range(4),repeat=2):
            qi=plus(preps[i],Ys[q][i],-ONE);qj=plus(preps[j],Ys[q][j],-ONE)
            assert(dot(Ys[q][i],Ys[q][j])+dot(qi,qj))/36==k(int(i==j))
            assert dot(Ys[q][i],qj)==ZERO
    check('all_reference_safe_query_effects_sum_identity',True)
    tracial=sum((F(7,9)*As[0][i][i]+F(2,9)*As[1][i][i] for i in range(3)),ZERO)/3
    assert tracial==k(F(2,27))
    return {'actual_coherent_native_columns':coherent,'B12':[[str(x) for x in row] for row in B.tolist()],
        'full_family_dimension':576,'inverse_columns':inverse_columns,
        'phi_numerator':[[i,a.data()] for i,a in sorted(phi.items())],'phi_denominator_squared':24,
        'prep_numerator':[[i,a.data()] for i,a in sorted(prep.items())],'prep_denominator':6,
        'projector_rank4_numerators':factors,'rank_denominator_squared':24,
        'query0_family_factors':{str(q):[vecdata(v) for v in Ys[q]] for q in (0,1)},
        'null_plants':nulls,'effects_and_controls':mixed,
        'whole_arm_readout':{'native_reference_query0':'0','HOLD_arm_query_reference_tracial_query0':'2/27',
            'native_active_tracial_query0':'2/27','native_active_tracial_query1':'25/27',
            'retained_count_inverse_return':'1','committed_count_Q_inverse_return':'53/81','committed_count_Q_inverse_nonreturn':'28/81'},
        'computed_tracial_query0_coefficients':tracial.data(),
        'source_unitary_definition':'U_A=exp(i beta J Qcount)(U0 tensor Icount), z=(3+4i)/5; each full inverse executes the displayed factors in reverse order.',
        'query_definition':'S_A=U_A†[Z_query chi Pi+X_query(I-chi Pi)]U_A; exact current query and W are used to compute chi and transfer family before, then decoded/uncomputed after.',
        'complete_output_factorization':'On active arm and input query0, q0 coefficient at count(k,l),family row r,input f is [h0(k,l)Y0_f(r)+h1(k,l)Y1_f(r)]/360; q1 is [prep_f(r)-same]/360. For input query1 exchange outputs and negate the P output. All3600 addresses and576 rows retained. Physical family decoding W† puts processor back in3 on prepared inputs.'}

def calculate():
    global CAR_CALLS
    CALLS.clear();CAR_CALLS=0;checks={}
    def check(name,value):
        checks[name]=bool(value)
        if not value:raise AssertionError(name)
    bindings=json.loads((HERE/'SOURCES.json').read_text())
    check('all_frozen_sources',all(sha(ROOT/p)==v for p,v in bindings.items()))
    check('literal_native_pin',sha(Path(native.__file__))==PIN)
    arm=arm_query(check);count=count_source(check);fam=family_query(check)
    count_p=F(sum(row[6] for row in count['joint_labels']),len(count['joint_labels']))
    assert (1-count_p)**2+count_p**2==F(53,81)
    # The actual reference query has chi0; the source-HOLD changed law's
    # executed query table has chi1. Compose those columns with the SAME
    # measured full-projector tracial probability, with no normalized success.
    tracial=K(tuple(map(F,fam['computed_tracial_query0_coefficients'])))
    for law in ('C','E'):
        for actual,held in zip(arm['query_rows'][law],arm['HOLD_query_rows'][law],strict=True):
            source_index,fl,b=actual[:3]
            if b==0 and arm['input_states'][source_index]['tagged_port']==7:
                assert actual[5][2]==0 and held[3][2]==1
                assert actual[5][2]*tracial==ZERO and held[3][2]*tracial==k(F(2,27))
    check('retained_count_and_committed_Q_control',True)
    check('charged_reference_HOLD_control',fam['whole_arm_readout']['native_reference_query0']!=fam['whole_arm_readout']['HOLD_arm_query_reference_tracial_query0'])
    return {'schema':'oph.row11.retained_family_oracle.v1','candidate':'A-retained-family-oracle-N1','class':'F with N+ native occupation/count components',
        'checks':checks,'all_pass':all(checks.values()),'source_bindings':bindings,
        'literal_native_calls':dict(CALLS),'total_literal_native_calls':sum(CALLS.values()),'CAR_C_E_column_evaluations':CAR_CALLS,
        'whole_arm_native_query':arm,'retained_count_source':count,'retained_family_query':fam,
        'new_M6_pulse_execution':False,'old_CP_instrument_lift_inferred':False,'physical_row_established':False,
        'scope':'Exact logical retained source factors and full inverse on the same tagged charged domain. Coherent source, auxiliary chart/preparation/analyzer, reflection, equal-charge W and witness beta are declared. No physical clock, minimum, beta selection or mass ratio.'}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--write',type=Path);ap.add_argument('--check',type=Path);args=ap.parse_args()
    start=time.perf_counter();r=calculate();blob=canonical(r)
    if args.write:
        with args.write.open('x') as f:f.write(blob)
    if args.check:assert args.check.read_text()==blob,'byte replay mismatch'
    print(canonical({'all_pass':r['all_pass'],'checks':len(r['checks']),'native_calls':r['total_literal_native_calls'],
        'coherent_native_columns':r['retained_family_query']['actual_coherent_native_columns']['actual_Observer_repair_columns'],
        'full_family_inverse_columns':r['retained_family_query']['inverse_columns'],
        'full_projector_columns':1152,'original_count_labels':3600,
        'reference_query0':r['retained_family_query']['whole_arm_readout'],
        'seconds':time.perf_counter()-start,'physical_row_established':False}))

if __name__=='__main__':main()
