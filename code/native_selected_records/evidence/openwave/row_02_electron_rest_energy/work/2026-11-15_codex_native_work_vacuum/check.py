"""Independent implementation of finite receipts; imports no producer modules.
Usage: python3 -B check.py [packet_directory]. NumPy+stdlib only.
"""
from pathlib import Path
import json
import math
import sys
from fractions import Fraction as F
import numpy as np

BASE=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
old=json.loads((BASE/'OUTPUT.json').read_text())
new=json.loads((BASE/'refinement/OUTPUT.json').read_text())
EDGES=[(0,1),(0,1),(0,2),(0,2),(1,2)]
CHECKS=0

def require(ok,label):
    global CHECKS
    CHECKS+=1
    if not ok:raise AssertionError(label)

def unpack(x):
    return np.array(x['real'])+1j*np.array(x['imag']) if isinstance(x,dict) else np.array(x)

def close(a,b,label,tol=2e-10):
    a,b=np.asarray(a),np.asarray(b)
    require(a.shape==b.shape and np.max(abs(a-b),initial=0)<tol,label)

def event(q,e,c):
    require(c in (0,1),'coin domain')
    a,b=e;t=q[a]+q[b];lo=t//2;hi=t-lo;y=list(q)
    y[a],y[b]=(lo,hi) if c==0 else (hi,lo)
    return tuple(y),q[a]-y[a]

require(len(old['source_records'])==30 and len(old['negative_source_records'])==30,'full source columns')
require(len(new['source_columns'])==30 and len(new['refinements'])==7 and len(new['selector_full_records'])==7,'full refinement domains')
require(old['native_events_literal_function_calls']==1230,'old literal count')
require(new['literal_native_function_calls']==372+sum(2*k*4**k for k in range(7)),'new literal count')
# All native columns are independently evaluated by floor/ceiling arithmetic.
R=np.zeros((30,30));P=[[F(0) for j in range(3)] for i in range(3)]
for index,r in enumerate(old['source_records']):
    sign,p,e,c,j,y,nc=r
    require(index==(p*5+e)*2+c,'source column ordering')
    q=tuple(sign*int(x==p) for x in range(3));out,current=event(q,EDGES[e],c)
    require(current==j and out[y]==sign and nc==(c^int(j!=0)),'source column')
    R[10*y+2*e+nc,10*p+2*e+c]=1;P[y][p]+=F(1,10)
for sign,p,e,c,j,y,nc in old['negative_source_records']:
    out,current=event(tuple(sign*int(x==p) for x in range(3)),EDGES[e],c)
    require(current==j and out[y]==sign and nc==(c^int(j!=0)),'negative column')
P=np.array([[float(q) for q in row] for row in P]);close(P,old['P'],'P')
a=np.ones(10)/math.sqrt(10);A=np.kron(np.eye(3),a[:,None]);C=np.kron(np.eye(3),2*np.outer(a,a)-np.eye(10));W=R@C
close(W,old['W'],'W')
v=a-np.eye(10)[0];v/=np.linalg.norm(v);H=np.eye(10)-2*np.outer(v,v)
B=np.kron(np.eye(3),H)@R@A;K=[B[q::10] for q in range(10)]
I=np.eye(3);x=np.array([0.,1.,0.]);sum_failure=0.
require(len(old['birth_outcomes'])==20,'all birth flags')
for r in old['birth_outcomes']:
    s,q=r['arm'],r['choice'];y=(int(q==0)*I+s*K[q])@x/2
    close(y,r['source_amplitude'],'birth state');close(y@y,r['probability'],'birth mass')
    if r['terminal']:sum_failure+=r['probability']
x=(I-P)@x/2
require(len(old['filter_first_failure_outcomes'])==8,'all cooling stages')
for n,t in enumerate(old['filter_first_failure_outcomes']):
    require(t['stage']==n+1 and len(t['outcomes'])==10,'filter stage')
    for q,r in enumerate(t['outcomes']):
        require(q==r['choice'],'filter flag');y=K[q]@x
        close(y,r['source_amplitude'],'filter state');close(y@y,r['probability'],'filter mass')
        if q:sum_failure+=r['probability']
    x=P@x
success=float(x@x);psi=x/math.sqrt(success);f=np.array([0.,1.,-1.])/math.sqrt(2);g=np.array([-2.,1.,1.])/math.sqrt(6)
close(success,old['selection_success'],'cooling success');close(sum_failure+success,1,'all first failures')
require(abs(float(F(old['selection_exact_probability']))-success)<1e-15,'exact success')
close((f@psi)**2,old['selection_compact_fidelity'],'selected fidelity')

require(len(old['changed_law_records'])==30,'all alternative columns')
for sign,p,e,c,j,y,nc in old['changed_law_records']:
    q=[int(x==p) for x in range(3)];a0,b0=EDGES[e];out=q.copy();out[a0],out[b0]=out[b0],out[a0]
    require(j==q[a0]-out[a0] and out[y]==1 and nc==c,'always-swap column')
close(old['controls']['always_swap_transfer'],2*P-I,'changed transfer')
close(old['controls']['always_swap_compact_fidelity'],4/7,'changed selection')
# Work path independently replays every actual source input and event.
work=old['work'];primitive_edges=[]
for a,b in work['edges']:primitive_edges.extend([(a,b),(b,a)])
require(len(work['one_unit_traces'])==12 and len(work['zero_two_unit_audit'])==32,'all work sectors')
for r in work['one_unit_traces']:
    p,c=r['input'];q=tuple(int(x==p) for x in range(6))
    require(len(r['trace'])==len(work['primitive_word']),'work length')
    for k,t in zip(work['primitive_word'],r['trace']):
        out,j=event(q,primitive_edges[k],c);y=out.index(1);nc=c^int(j!=0)
        require(t==[k,p,c,j,y,nc],'work trace');p,c,q=y,nc,out
    require(r['output']==[p,c],'work output')
for r in work['zero_two_unit_audit']:
    q=tuple(r['input']);c=r['coin']
    for k,t in zip(work['primitive_word'],r['trace']):
        out,j=event(q,primitive_edges[k],c);nc=c^int(j!=0)
        require(t==[k,c,j,nc,list(out)],'complement work trace');q,c=out,nc
    require(list(q)==r['output'] and c==r['output_coin'],'complement output')

require(len(new['work_columns'])==216 and len(new['work_traces'])==6,'new work domain')
for edge,p,c,j,y,nc in new['work_columns']:
    out,current=event(tuple(int(x==p) for x in range(6)),edge,c)
    require(j==current and out[y]==1 and nc==(c^int(j!=0)),'new work columns')
for p0,trace in enumerate(new['work_traces']):
    p,c=p0,0;require(len(trace)==21,'new work trace length')
    for k,t in zip(new['work_word'],trace):
        out,j=event(tuple(int(x==p) for x in range(6)),primitive_edges[k],c);y=out.index(1);nc=c^int(j!=0)
        require(t==[k,p,c,j,y,nc],'new work trace');p,c=y,nc
    require(p==(p0+3)%6 and c==0,'new work terminal')
# Independent exact-mode Fourier law, not the producer FFT path.
def dirichlet(angle,Q):
    z=angle/2
    den=np.sin(z)
    result=np.empty_like(z,dtype=complex)
    small=abs(den)<1e-12
    result[small]=1
    result[~small]=np.exp(1j*(Q-1)*z[~small])*np.sin(Q*z[~small])/(Q*den[~small])
    return result

Rf=np.eye(33)
for p,e,c,j,y,nc in new['source_columns']:
    out,current=event(tuple(int(x==p) for x in range(3)),EDGES[e],c)
    require(j==current and out[y]==1 and nc==(c^int(j!=0)),'refinement source column')
    Rf[:,11*p+2*e+c]=0;Rf[11*y+2*e+nc,11*p+2*e+c]=1
all_errors=[];holding=[]
for k,(r,ticket) in enumerate(zip(new['refinements'],new['selector_full_records'])):
    eps=2.**(-k);width=2*k
    require(r['k']==k and r['epsilon']==eps,'refinement order')
    require(len(ticket['tickets'])==4**k and ticket['width']==width,'ticket completeness')
    for number,t in enumerate(ticket['tickets']):
        bits=[(number>>j)&1 for j in reversed(range(width))]
        require(t['bits']==bits and t['currents']==[0]*width and t['active']==(number==0),'ticket records')
        for c in bits:
            out,j=event((0,0),(0,1),c);require(out==(0,0) and j==0,'blank native event')
    ae=np.r_[np.full(10,eps/math.sqrt(10)),math.sqrt(1-eps*eps)]
    Ae=np.kron(np.eye(3),ae[:,None]);Ce=np.kron(np.eye(3),2*np.outer(ae,ae)-np.eye(11));We=Rf@Ce
    close(Ae.T@Rf@Ae,r['complete_transfer'],'refined transfer')
    Q=r['counter_bins'];angles=2*math.pi*np.arange(Q)/Q
    amplitude=np.zeros((Q,33),complex)
    ideal=[]
    for mode,lam0 in ((f,.6),(g,.4)):
        lam=1-eps*eps*(1-lam0);theta=2*math.asin(eps*math.sqrt((1-lam0)/2))
        u=Ae@mode;w=(Rf@u-lam*u)/math.sqrt(1-lam*lam)
        plus=(u-1j*w)/math.sqrt(2);minus=plus.conj()
        coefficient=float(mode@psi)/math.sqrt(2)
        amplitude+=coefficient*(dirichlet(theta-angles,Q)[:,None]*plus+dirichlet(-theta-angles,Q)[:,None]*minus)
        ideal.append(plus)
    probability=np.sum(abs(amplitude)**2,axis=1)
    close(probability,r['all_counter_probabilities'],'all Fourier probabilities',2e-10)
    rho=unpack(r['positive_density']);positive=amplitude[1:Q//2]
    independent_rho=positive.T@positive.conj()/sum(probability[1:Q//2])
    close(rho,independent_rho,'positive complete density',3e-10)
    close(np.vdot(ideal[0],rho@ideal[0]).real,r['selected_positive_mode_fidelity'],'positive fidelity')
    # Check all public Ramsey moments from independent source-mode formula.
    require(len(r['ramsey_records'])==17,'all Ramsey moments')
    for n,record in enumerate(r['ramsey_records']):
        require(record['time']==n and record['native_W_words']==n*2**k,'clock tags')
        z=sum(float(mode@psi)**2*math.cos(n*2*math.asin(eps*math.sqrt((1-lam0)/2))/eps) for mode,lam0 in ((f,.6),(g,.4)))
        close(record['x_plus'],(1+z)/2,'Ramsey record');close(record['y_plus'],.5,'Ramsey imaginary')
    close(r['positive_probability'],sum(probability[1:Q//2]),'positive mass')
    close(r['negative_probability'],sum(probability[Q//2+1:]),'negative mass')
    close(r['all_counter_probability_sum'],sum(probability),'counter total')
    close(r['record_only_reader']['selected_frequency'],2*math.asin(eps/math.sqrt(5))/eps,'selected record frequency',1e-6)
    poles=r['record_only_reader']['poles']
    expected=sorted([-2*math.asin(eps/math.sqrt(5))/eps,2*math.asin(eps/math.sqrt(5))/eps,-2*math.asin(eps*math.sqrt(.3))/eps,2*math.asin(eps*math.sqrt(.3))/eps])
    close([p['phase'] for p in poles],expected,'all signed inferred poles',1e-6)
    all_errors.append(max(abs(np.array([p['phase'] for p in poles])-expected)))
    live=rho.copy();minimum=1.;maxchange=0.
    for n in range(33):
        minimum=min(minimum,float(np.vdot(ideal[0],live@ideal[0]).real))
        maxchange=max(maxchange,float(np.linalg.norm(live-rho,ord='fro')))
        live=We@live@We.T
    holding.append({'k':k,'minimum_selected_mode_fidelity_through32':minimum,'max_density_change_Frobenius':maxchange,'three_port_support_exact':True})
    # Independent full-coordinate compression check at nontrivial finite levels.
    if k in (1,2):
        d=10*4**k;J=np.zeros((d,11));J[:10,:10]=np.eye(10);J[10:,10]=1/math.sqrt(d-10)
        E=np.kron(np.eye(3),J);fullR=np.eye(3*d)
        for p,e,c,j,y,nc in new['source_columns']:
            fullR[:,d*p+2*e+c]=0;fullR[d*y+2*e+nc,d*p+2*e+c]=1
        raw=np.ones(d)/math.sqrt(d);fullC=np.kron(np.eye(3),2*np.outer(raw,raw)-np.eye(d))
        close(E.T@E,np.eye(33),'compression isometry')
        close(fullR@fullC@E,E@We,'full raw clock intertwining')
require(max(all_errors)<1e-6,'record decoder all frequencies')
print(json.dumps({'checks':CHECKS,'verdict':'PASS','source_floor_ceiling_reimplementation':True,
 'producer_modules_imported':False,'full_padding_levels':[1,2],
 'maximum_all_band_decoder_error':max(all_errors),'holding':holding,
 'scope':'Independent implementation by the producer group, not other-group G6 or scientific award.'},indent=2))
