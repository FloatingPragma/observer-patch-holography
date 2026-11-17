"""Independent row07 review: no imports from the packet under review."""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal as D, localcontext
import hashlib,json,math
import numpy as np

HERE=Path(__file__).resolve().parent
P=HERE.parents[1]/'work/2026-11-15_codex_pr723_rotation'
EXPECTED='05af1f898eb98a68413084b3b01586e600c79aa231d2f0318f0f90b887158f32'
assert hashlib.sha256((P/'MANIFEST.sha256').read_bytes()).hexdigest()==EXPECTED
for directory in (P,P/'berry'):
    for line in (directory/'MANIFEST.sha256').read_text().splitlines():
        h,name=line.split('  ',1)
        assert hashlib.sha256((directory/name).read_bytes()).hexdigest()==h


def endpoint(q,a,b,c,law):
    assert type(c)is int and c in (0,1)
    y=list(q)
    if law=='native':
        lo=(q[a]+q[b])//2;hi=q[a]+q[b]-lo
        y[a],y[b]=(lo,hi) if c==0 else (hi,lo)
    elif law=='move_always':y[a],y[b]=q[b],q[a]
    elif law!='hold':raise AssertionError(law)
    return tuple(y),q[a]-y[a]


raw=json.loads((P/'RECORDS.json').read_text())
states=[(0,)*12]+[tuple(int(i==p)-int(i==r) for i in range(12)) for p in range(12) for r in range(12) if p!=r]
index={q:i for i,q in enumerate(states)}
edges=json.loads((P/'seams.json').read_text())
old_counts={}
for law,table in raw.items():
    seen=set();counts=np.zeros((133,133),dtype=int)
    for z,e,c,j,out in table:
        assert (z,e,c)not in seen;seen.add((z,e,c))
        y,current=endpoint(states[z],*edges[e],c,law)
        assert j==current and out==index[y]
        counts[z,out]+=1
    assert len(seen)==133*30*2
    old_counts[law]=counts
native=old_counts['native']
assert np.all(np.diag(native)>np.sum(native,axis=1)-np.diag(native))
assert np.all(np.diag(native)[1:]==50) and native[0,0]==60

B=((1,-1,0),(1,1,-2))
TEDGES=((0,1),(1,2),(2,0))
cache={}
total_body=0


def triangle(weights,table,law='native'):
    global total_body
    tickets=[e for e,n in enumerate(weights) for _ in range(n)]
    assert len(tickets)==64
    seen=set();counts=np.zeros((3,3),dtype=int)
    fine=np.zeros((384,3,3),dtype=int)
    for z,t,c,j,out in table:
        assert (z,t,c)not in seen;seen.add((z,t,c))
        q=tuple(int(i==z) for i in range(3))
        y,current=endpoint(q,*TEDGES[tickets[t]],c,law)
        assert current==j and y.index(1)==out and -1<=j<=1
        counts[out,z]+=1
        fine[(2*t+c)*3+j+1,out,z]=1
    assert len(seen)==384
    assert np.array_equal(np.einsum('rai,raj->ij',fine,fine),128*np.eye(3,dtype=int))
    M=[[F(sum(B[i][a]*int(counts[a,b])*B[k][b] for a in range(3) for b in range(3)),128*(2 if i==0 else 6)) for k in range(2)] for i in range(2)]
    w=np.ones(384)/math.sqrt(384);u=w.copy();u[0]-=1;u/=np.linalg.norm(u)
    K=fine/math.sqrt(128)
    E=K-2*u[:,None,None]*np.einsum('r,rij->ij',u,K)
    assert np.max(abs(E[0]-counts/128/math.sqrt(3)))<1e-13
    assert np.max(abs(np.einsum('rai,raj->ij',E,E)-np.eye(3)))<1e-13
    total_body+=len(table)
    if law=='native':cache[tuple(weights)]=(M,E)
    return M


berry=P/'berry'
rr=json.loads((berry/'RECEIPT.json').read_text())
rd=json.loads((berry/'RECORDS.json').read_text())
for law in ('native','move_always','hold'):
    for w,table in zip(rr['weight_loop'],rd['source'][law]['loop']):triangle(w,table,law)
    triangle((32,16,16),rd['source'][law]['reference'],law)
assert len(rd['native_scheduler'])==384
for ticket,bit,c,j in rd['native_scheduler']:
    assert c==(ticket>>bit)&1 and j==0


def multiply(v,M):return tuple(sum(M[i][j]*v[j] for j in range(2)) for i in range(2))
def norm2(v):return v[0]*v[0]/3+v[1]*v[1]
def scalar(v,w):return v[0]*w[0]/3+v[1]*w[1]
def log_fraction(x):return math.log(x.numerator)-math.log(x.denominator)
def products(sequence,reference=None):
    v=(F(0),F(1));w=v
    if reference is None:reference=[(32,16,16)]*len(sequence)
    for a,b in zip(sequence,reference):v=multiply(v,cache[tuple(a)][0]);w=multiply(w,cache[tuple(b)][0])
    vis=2*scalar(v,w)/(norm2(v)+norm2(w))
    success=(norm2(v)+norm2(w))/(2*3**len(sequence))
    return {'visibility':float(vis),'log10_success':log_fraction(success)/math.log(10),
            'overlap':float(v[1])/math.sqrt(float(norm2(v))),'coefficient_sign':(v[1]>0)-(v[1]<0)}


rotor=json.loads((berry/'ROTOR_RECEIPT.json').read_text())
rot_records=json.loads((berry/'ROTOR_RECORDS.json').read_text())
sequences={}
for coupling in ('polar','apolar'):
    rows=[r for r in rotor['geometry_records'] if r['coupling']==coupling]
    assert [r['step'] for r in rows]==list(range(1,97))
    sequences[coupling]=[tuple(r['native_ticket_weights']) for r in rows]
    for row,table in zip(rows,rot_records[coupling]):triangle(row['native_ticket_weights'],table)
signs={}
for coupling,seq in sequences.items():
    signs[coupling]=products(seq)
    want=rotor['rotor_products'][coupling]['1']
    assert abs(signs[coupling]['visibility']-want['relative_interference_visibility'])<2e-10
    assert abs(signs[coupling]['log10_success']-want['log10_complete_success_probability'])<2e-8
two=products(sequences['polar']*2)
two_saved=json.loads((berry/'ROTOR_AUDIT.json').read_text())['two_turns']
assert abs(two['visibility']-two_saved['relative_visibility'])<2e-10
assert abs(two['log10_success']-two_saved['log10_original_success_probability'])<2e-8
assert signs['polar']['coefficient_sign']==-1 and signs['apolar']['coefficient_sign']==1 and two['coefficient_sign']==1

# Recompute complete first-rejection probabilities using exact rational optical
# state at every stage, not repeated floating multiplication or eigenvectors.
def outcomes(sequence,reference,stored=None):
    a=(F(0),F(1));b=a
    max_prob_error=0.;max_reach_error=0.;failure=0.
    for i,(x,y) in enumerate(zip(sequence,reference)):
        na,nb=norm2(a),norm2(b)
        reach=(na+nb)/(2*3**i)
        wa=float(na/(na+nb));wb=1-wa
        def vector(v):
            denominator=math.sqrt(6*float(norm2(v)))
            return np.array([float(v[0])*B[0][k]+float(v[1])*B[1][k] for k in range(3)])/denominator
        va,vb=vector(a),vector(b)
        ax=np.einsum('rij,j->ri',cache[tuple(x)][1],va)
        bx=np.einsum('rij,j->ri',cache[tuple(y)][1],vb)
        probs=wa*np.sum(ax*ax,axis=1)+wb*np.sum(bx*bx,axis=1)
        assert abs(sum(probs)-1)<1e-12 and min(probs)>-1e-15
        logreach=log_fraction(reach)/math.log(10)
        if stored:
            stage=stored['stages'][i]
            assert stage['stage']==i+1
            max_prob_error=max(max_prob_error,max(abs(probs-np.array(stage['conditional_probabilities_by_eraser_flag0_through383']))))
            max_reach_error=max(max_reach_error,abs(logreach-stage['log10_probability_of_reaching_stage']))
        failure+=float(reach)*sum(probs[1:])
        a=multiply(a,cache[tuple(x)][0]);b=multiply(b,cache[tuple(y)][0])
    success=(norm2(a)+norm2(b))/(2*3**len(sequence))
    assert abs(failure+float(success)-1)<1e-12
    if stored:
        assert len(stored['stages'])==len(sequence)
        assert max_prob_error<2e-10 and max_reach_error<2e-8
        assert abs(stored['log10_final_success_probability']-log_fraction(success)/math.log(10))<2e-8
    return {'stages':len(sequence),'flag_probabilities_checked':384*len(sequence),
            'max_probability_error':max_prob_error,'max_reach_log_error':max_reach_error,
            'log10_success':log_fraction(success)/math.log(10),'all_first_rejection_mass':failure}


outcome_saved=json.loads((berry/'ROTOR_OUTCOMES.json').read_text())
checks={}
for name,seq in [('polar96',sequences['polar']),('apolar96',sequences['apolar']),
                 ('polar192',[w for w in sequences['polar'] for _ in range(2)])]:
    checks[name]=outcomes(seq,[(32,16,16)]*len(seq),outcome_saved[name])
checks['repeated_polar96_twice']=outcomes(sequences['polar']*2,[(32,16,16)]*192)
base_loop=[tuple(w) for w in rr['weight_loop']]
all_outcomes=json.loads((berry/'ALL_OUTCOMES.json').read_text())
for name,seq in [('original48',base_loop),('original192',[w for w in base_loop for _ in range(4)])]:
    checks[name]=outcomes(seq,[(32,16,16)]*len(seq),all_outcomes[name])
vertices=[(32,16,16),(16,32,16),(32,16,16),(16,16,32),(16,32,16),(16,16,32),(32,16,16)]
zero=[tuple(a[k]+step*(b[k]-a[k])//16 for k in range(3))
      for a,b in zip(vertices,vertices[1:]) for step in range(1,17)]
one=[w for w in base_loop for _ in range(2)]
assert sorted(zero)==sorted(one)==sorted(base_loop*2)
checks['balanced96_zero_winding']=outcomes(one,zero,all_outcomes['balanced96_zero_winding'])
balanced=json.loads((berry/'BALANCED_RECEIPT.json').read_text())
for name,reference in [('zero_winding',zero),('two_windings',base_loop*2)]:
    measured=products(one,reference)
    assert abs(measured['visibility']-balanced['summaries'][name]['relative_visibility'])<2e-10
    assert (measured['visibility']<-.9)==balanced['balanced_fringe_results'][name]
amendment=P.parent/'2026-11-15_codex_pr723_two_turn_rejections'
amendment_hash='1861b420260787509cc944f291919af901a87a7f5dc7ca8c86c91c7ab4cfe512'
assert hashlib.sha256((amendment/'MANIFEST.sha256').read_bytes()).hexdigest()==amendment_hash
for line in (amendment/'MANIFEST.sha256').read_text().splitlines():
    h,name=line.split('  ',1)
    assert hashlib.sha256((amendment/name).read_bytes()).hexdigest()==h
checks['repeated_polar96_twice_amendment']=outcomes(sequences['polar']*2,[(32,16,16)]*192,
    json.loads((amendment/'TWO_TURN_OUTCOMES.json').read_text()))

# A constant component is not rejected by the source transfer. Test its effect
# algebraically on this same declared complete instrument, without attributing
# an alternative preparation to the frozen executed run.
va=(F(0),F(1));vb=va
for w in sequences['polar']:
    va=multiply(va,cache[w][0]);vb=multiply(vb,cache[(32,16,16)][0])
epsilon=F(1,10**18)
leak_visibility=2*(scalar(va,vb)+epsilon**2)/(norm2(va)+norm2(vb)+2*epsilon**2)
assert float(leak_visibility)>.99

import copy
mutation_checks={}
for name,edit in [('invalid_coin',lambda t:t[0].__setitem__(2,2)),
                  ('wrong_current_flag',lambda t:t[0].__setitem__(3,t[0][3]+1))]:
    bad=copy.deepcopy(rd['source']['native']['reference']);edit(bad)
    try:triangle((32,16,16),bad)
    except AssertionError:mutation_checks[name]=True
    else:mutation_checks[name]=False
bad=copy.deepcopy(json.loads((amendment/'TWO_TURN_OUTCOMES.json').read_text()))
bad['stages'][0]['conditional_probabilities_by_eraser_flag0_through383'][1]+=.01
try:outcomes(sequences['polar']*2,[(32,16,16)]*192,bad)
except AssertionError:mutation_checks['wrong_rejection_probability']=True
else:mutation_checks['wrong_rejection_probability']=False
assert all(mutation_checks.values())

PI=D('3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679')
def cos(x):
    x=x%(2*PI)
    if x>PI:x-=2*PI
    s=term=D(1)
    for n in range(1,400):
        term*=-x*x/D((2*n-1)*(2*n));s+=term
        if abs(term)<D('1e-85'):return s
    raise AssertionError('series')
ties=[]
with localcontext()as ctx:
    ctx.prec=90
    for row in rotor['geometry_records']:
        winding=1 if row['coupling']=='polar' else 2
        theta=D(row['step'])*2*PI/D(96)
        real=[D(64)*(1+cos(winding*theta-D(e)*2*PI/3)/2)/3 for e in range(3)]
        for i,x in enumerate(real):
            if abs(x-x.to_integral_value())<D('1e-70'):real[i]=x.to_integral_value()
        floors=[int(x//1) for x in real]
        frac=[x-f for x,f in zip(real,floors)]
        added=[w-f for w,f in zip(row['native_ticket_weights'],floors)]
        assert sum(row['native_ticket_weights'])==64 and all(v in (0,1) for v in added)
        for a in range(3):
            for b in range(3):
                if added[a] and not added[b]:
                    assert frac[a]>=frac[b]-D('1e-65')
                    if abs(frac[a]-frac[b])<D('1e-65'):ties.append([row['coupling'],row['step'],a,b])

result={'bound_manifest':EXPECTED,'source_tables':{'signed_native':7980,'signed_changed':15960,
        'triangle_and_rotor_body_records':total_body,'triangle_native_body':49*384,'triangle_changed_body':2*49*384,
        'native_ticket_records':384,'rotor_native':2*96*384},
        'native_calls_total':7980+49*384+384+2*96*384,
        'exact_signs':signs,'exact_two_turn':two,'all_outcome_audits':checks,
        'bound_amendment_manifest':amendment_hash,'mutation_checks':mutation_checks,
        'constant_input_amplitude_sensitivity':{'amplitude':'1e-18','visibility':float(leak_visibility),
          'scope':'counterfactual preparation sensitivity of the declared instrument; no new native calls'},
        'decimal_ties':ties,'isometry_and_complete_eraser':True,'source_selection_derived':False,
        'physical_rotation_derived':False,'same_group':True,'G6_pass':False,'G9_pass':False}
(HERE/'AUDIT.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
print(json.dumps({'all_checks':True,'native_calls_total':result['native_calls_total'],'polar':signs['polar'],'two_turn':two,'tie_count':len(ties)},sort_keys=True))
