"""Native blank-ticket clock refinement; complete finite record populations.
Self-contained NumPy+stdlib, no file/network/process access or packet imports.
"""
import itertools
import json
import math
import numpy as np

CALLS=0
EDGES=[(0,1),(0,1),(0,2),(0,2),(1,2)]

def native(q,e,c):
    global CALLS
    CALLS+=1
    a,b=e;d=q[a]-q[b];j=(d-(2*c-1)*(abs(d)%2))//2
    out=list(q);out[a]-=j;out[b]+=j
    return tuple(out),(c,j)

def encode(x):
    x=np.asarray(x)
    return {'real':x.real.tolist(),'imag':x.imag.tolist()} if np.iscomplexobj(x) else x.tolist()

# Folded source arcs: ten actual active microchoices and one inactive aggregate.
R=np.eye(33);source=[]
for p in range(3):
    for e,edge in enumerate(EDGES):
        for c in (0,1):
            q=tuple(int(x==p) for x in range(3));out,(coin,j)=native(q,edge,c)
            y=out.index(1);nc=c^int(j!=0);i=11*p+2*e+c;o=11*y+2*e+nc
            R[:,i]=0;R[o,i]=1;source.append([p,e,c,j,y,nc])
assert np.max(abs(R@R-np.eye(33)))==0
avec=np.r_[np.ones(10)/math.sqrt(10),0.]
A=np.kron(np.eye(3),avec[:,None]);P=A.T@R@A
initial=np.array([0.,1.,0.]);s=(np.eye(3)-P)@initial/2
for n in range(8):s=P@s
cooling_probability=float(s@s);s/=math.sqrt(cooling_probability)
f=np.array([0.,1.,-1.])/math.sqrt(2)

# Work compiler on six physical positions with an independent native work coin.
work_edges=[(0,1),(0,2),(1,2),(3,4),(3,5),(4,5),(0,3),(1,4),(2,5)]
primitives=[];work_columns=[]
for e in work_edges:
    for edge in (e,tuple(reversed(e))):
        perm=[]
        for p in range(6):
            for c in (0,1):
                q=tuple(int(x==p) for x in range(6));out,(coin,j)=native(q,edge,c)
                y=out.index(1);nc=c^int(j!=0);perm.append(2*y+nc)
                work_columns.append([list(edge),p,c,j,y,nc])
        moved=[i for i in range(12) if perm[i]!=i]
        primitives.append((edge,perm,moved))
adj={i:[] for i in range(12)}
for k,(e,perm,moved) in enumerate(primitives):
    a,b=moved;adj[a].append((b,k));adj[b].append((a,k))
word=[]
for x,y in ((0,3),(1,4),(2,5)):
    start,goal=2*x,2*y;todo=[start];prev={start:None}
    for p in todo:
        for q,k in adj[p]:
            if q not in prev:prev[q]=(p,k);todo.append(q)
    path=[];p=goal
    while p!=start:p,k=prev[p];path.append(k)
    path.reverse();word+=path+path[-2::-1]
permutation=[];work_traces=[]
for p0 in range(6):
    p,c=p0,0;trace=[]
    for k in word:
        q=tuple(int(x==p) for x in range(6));out,(coin,j)=native(q,primitives[k][0],c)
        y=out.index(1);nc=c^int(j!=0);trace.append([k,p,c,j,y,nc]);p,c=y,nc
    assert c==0 and p==(p0+3)%6
    permutation.append(p);work_traces.append(trace)
transport=np.zeros((66,66))
for p,y in enumerate(permutation):
    for h in range(11):transport[11*y+h,11*p+h]=1

# This reader uses only measured binary Ramsey probabilities, not W/eigenvectors.
def read_poles(records):
    z=np.array([(2*r['x_plus']-1)+1j*(2*r['y_plus']-1) for r in records])
    mat=np.array([[z[n+j] for j in range(4)] for n in range(4)])
    coefficients=np.linalg.solve(mat,-z[4:8])
    poles=np.roots(np.r_[1.,coefficients[::-1]])
    order=np.argsort(np.angle(poles));poles=poles[order]
    weights=np.linalg.solve(np.array([[p**n for p in poles] for n in range(4)]),z[:4])
    predicted=np.array([sum(w*p**n for p,w in zip(poles,weights)) for n in range(len(z))])
    rows=[{'phase':float(np.angle(p)),'radius':float(abs(p)),
           'weight':[float(w.real),float(w.imag)]} for p,w in zip(poles,weights)]
    positive=[r for r in rows if r['phase']>0 and r['weight'][0]>.01]
    selected=min(positive,key=lambda r:r['phase'])
    return {'poles':rows,'selected_frequency':selected['phase'],
            'fit_condition_number':float(np.linalg.cond(mat)),
            'heldout_max_error':float(max(abs(predicted[8:]-z[8:])))}

rows=[];all_tickets=[]
for k in range(7):
    epsilon=2.**(-k);width=2*k;tickets=[]
    for bits in itertools.product((0,1),repeat=width):
        currents=[]
        for c in bits:
            out,(coin,j)=native((0,0),(0,1),c)
            assert out==(0,0) and coin==c and j==0
            currents.append(j)
        tickets.append({'bits':list(bits),'currents':currents,'active':not any(bits)})
    all_tickets.append({'k':k,'width':width,'equal_probability_denominator':4**k,'tickets':tickets})
    ae=np.r_[np.full(10,epsilon/math.sqrt(10)),math.sqrt(1-epsilon**2)]
    Ae=np.kron(np.eye(3),ae[:,None]);Pe=Ae.T@R@Ae
    Ce=np.kron(np.eye(3),2*np.outer(ae,ae)-np.eye(11));W=R@Ce
    assert np.max(abs(Pe-(np.eye(3)-epsilon**2*(np.eye(3)-P))))<1e-12
    u=Ae@s;stride=2**k;clock=np.eye(33,dtype=complex)
    for j in range(stride):clock=W@clock
    raw=[];live=u.astype(complex)
    for n in range(17):
        z=np.vdot(u,live)
        raw.append({'time':n,'native_W_words':n*stride,'x_plus':float((1+z.real)/2),'y_plus':float((1+z.imag)/2)})
        live=clock@live
    decoded=read_poles(raw)
    # Complete generic phase counter; no eigenphase is supplied to its circuit.
    Q=256*stride;waves=[];live=u.astype(complex)
    for n in range(Q):waves.append(live.copy());live=W@live
    amps=np.fft.fft(np.asarray(waves),axis=0)/Q
    probabilities=np.sum(abs(amps)**2,axis=1)
    pos=amps[1:Q//2];positive_probability=float(sum(probabilities[1:Q//2]))
    rho=pos.T@pos.conj()/positive_probability
    lam=1-.4*epsilon**2;theta=math.acos(lam)
    uf=Ae@f;wf=(R@uf-lam*uf)/math.sqrt(1-lam*lam);plus=(uf-1j*wf)/math.sqrt(2)
    fidelity=float(np.vdot(plus,rho@plus).real)
    # Apply the actual native copy permutation to the entire selected density.
    common=np.kron(np.eye(2),W);before=np.zeros((66,66),complex);before[:33,:33]=rho
    after=transport@before@transport.T
    phase_before=np.trace(before@common);phase_after=np.trace(after@common)
    joint=np.kron(plus,plus).reshape(33,33);joint_after=W@joint@W.T
    pair_error=float(np.linalg.norm(joint_after-np.exp(2j*theta)*joint))
    expected=2*math.asin(epsilon/math.sqrt(5))/epsilon
    error=decoded['selected_frequency']-2/math.sqrt(5)
    rows.append({'k':k,'epsilon':epsilon,'active_probability':epsilon**2,
      'raw_full_choice_count':10*4**k,'folded_choice_count':11,
      'orthogonal_inactive_padding_count':max(0,10*(4**k-1)-1),
      'source_selector_native_events_per_preparation':width,
      'complete_transfer':encode(Pe),'ramsey_records':raw,'record_only_reader':decoded,
      'predicted_frequency':expected,'limit_frequency':2/math.sqrt(5),
      'analytic_frequency_excess_upper_bound':epsilon**2/24,
      'decoded_frequency_error_from_prediction':decoded['selected_frequency']-expected,
      'counter_bins':Q,'all_counter_probabilities':probabilities.tolist(),
      'positive_probability':positive_probability,'negative_probability':float(sum(probabilities[Q//2+1:])),
      'positive_density':encode(rho),'selected_positive_mode_fidelity':fidelity,
      'single_word_phase':theta,'pair_word_phase':2*theta,'pair_direct_residual':pair_error,
      'work_phase_before':[float(phase_before.real),float(phase_before.imag)],
      'work_phase_after':[float(phase_after.real),float(phase_after.imag)],
      'work_commutator':float(np.max(abs(transport@common-common@transport))),
      'all_counter_probability_sum':float(sum(probabilities)),
      'integer_shift_alias_frequency':expected+2*math.pi/epsilon,
      'finite_frequency_bound':1.2,
      'gates':{'reader_prediction':abs(decoded['selected_frequency']-expected)<1e-5,
       'heldout_reader':decoded['heldout_max_error']<1e-8,
       'source_bound':-1e-10<=expected-2/math.sqrt(5)<=epsilon**2/24+1e-10,
       'positive_mode':fidelity>.99,'probability_complete':abs(sum(probabilities)-1)<1e-10,
       'pair_direct':pair_error<1e-10,'work_phase':abs(phase_before-phase_after)<1e-12,
       'alias_excluded_by_uniform_bound':expected+2*math.pi/epsilon>1.2}})
output={'candidate':'A-native-triangle-refined-work-v1','source_columns':source,
 'selector_full_records':all_tickets,'refinements':rows,'work_columns':work_columns,
 'work_word':word,'work_traces':work_traces,'work_events':len(word),
 'source_cooling_probability':cooling_probability,'literal_native_function_calls':CALLS,
 'separated_pair_anchored_keV':1022.,'second_band_limit_predicted_keV':511*math.sqrt(1.5),
 'coherence':'Unread current reversal, source-choice reflection and Fourier/Ramsey memory are declared services.',
 'clock':'epsilon per completed W word;2k native selector preparation attempts are separately charged.',
 'scope':'Source-generated finite-frequency refinement and selected localized excitation. No hardware duration, inertia, ambient mass gap, or external G6 claim.'}
print(json.dumps(output,indent=2,sort_keys=True,default=lambda x:x.item()))
assert all(all(r['gates'].values()) for r in rows)
