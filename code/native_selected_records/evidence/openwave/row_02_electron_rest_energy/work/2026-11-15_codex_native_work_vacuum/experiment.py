"""Finite native triangle selection/clock/exchange. Offline; stdout is the receipt.
Only NumPy and the standard library. No files, network, or process operations.
Native columns are literal calls; composed instruments are population execution.
"""
import json
import math
from collections import deque
from fractions import Fraction as F
import numpy as np

SEAMS = [(0,1),(0,1),(0,2),(0,2),(1,2)]
D, STEPS, BINS = 10, 8, 256
CALLS = 0

def native(loads, seam, coin):
    # Literal source function body at source SHA dd7c87...d0c811d.
    global CALLS
    CALLS += 1
    a,b = seam
    d = loads[a] - loads[b]
    current = (d - (2*coin-1)*(abs(d)%2)) // 2
    out = list(loads)
    out[a] -= current
    out[b] += current
    return tuple(out), (coin,current)

def encode(z):
    z=np.asarray(z)
    if np.iscomplexobj(z):
        return {'real':z.real.tolist(),'imag':z.imag.tolist()}
    return z.tolist()

def columns(sign=1, changed=False):
    R=np.zeros((30,30)); records=[]
    for x in range(3):
        for e,seam in enumerate(SEAMS):
            for coin in (0,1):
                q=tuple(sign*int(p==x) for p in range(3))
                out,(c,j)=native(q,seam,coin)
                if changed:
                    # Competing reversible always-swap law, same seams/coin alphabet.
                    out=list(q);a,b=seam;out[a],out[b]=out[b],out[a]
                    out=tuple(out);j=q[a]-out[a]
                    newcoin=coin
                else:
                    newcoin=coin ^ int(j!=0)
                y=next(p for p in range(3) if out[p])
                R[y*D+2*e+newcoin,x*D+2*e+coin]=1
                records.append([sign,x,e,coin,j,y,newcoin])
    return R,records

R,records=columns()
Rn,negative_records=columns(-1)
Rc,control_records=columns(changed=True)
a=np.ones(D)/math.sqrt(D)
A=np.kron(np.eye(3),a[:,None])
v=a-np.eye(D)[0];v/=np.linalg.norm(v)
H=np.eye(D)-2*np.outer(v,v)
C=np.kron(np.eye(3),2*np.outer(a,a)-np.eye(D))
P=A.T@R@A
Pcontrol=A.T@Rc@A
B=np.kron(np.eye(3),H)@R@A
K=[B[q::D] for q in range(D)]
W=R@C
Wcontrol=Rc@C
I=np.eye(3)
birth=[(int(s),q,(int(q==0)*I+s*K[q])/2) for s in (1,-1) for q in range(D)]
source_checks={
 'R_involution':float(np.max(abs(R@R-np.eye(30)))),
 'W_unitarity':float(np.max(abs(W.T@W-np.eye(30)))),
 'negative_same_transfer':float(np.max(abs(A.T@Rn@A-P))),
 'P_expected':float(np.max(abs(P-np.array([[.6,.2,.2],[.2,.7,.1],[.2,.1,.7]])))),
 'filter_completeness':float(np.max(abs(sum(k.T@k for k in K)-I))),
 'birth_completeness':float(np.max(abs(sum(k.T@k for s,q,k in birth)-I)))
}
f=np.array([0,1,-1])/math.sqrt(2)
g=np.array([-2,1,1])/math.sqrt(6)
omega=np.ones(3)/math.sqrt(3)
initial=np.array([0.,1.,0.])
first=[]
for s,q,k in birth:
    y=k@initial
    first.append({'arm':s,'choice':q,'probability':float(y@y),'terminal':not(s==-1 and q==0),'source_amplitude':y.tolist()})
state=(I-P)@initial/2
failures=[]
for n in range(STEPS):
    outcomes=[]
    for q,k in enumerate(K):
        y=k@state
        outcomes.append({'choice':q,'probability':float(y@y),'source_amplitude':y.tolist()})
    failures.append({'stage':n+1,'outcomes':outcomes})
    state=P@state
survival=float(state@state)
psi=state/math.sqrt(survival)
selected_fidelity=float((f@psi)**2)
selection_total=sum(r['probability'] for r in first if r['terminal'])+sum(r['probability'] for t in failures for r in t['outcomes'] if r['choice']!=0)+survival
u=A@psi
waves=[];live=u.astype(complex)
for n in range(BINS):
    waves.append(live.copy());live=W@live
phase_amplitudes=np.fft.fft(np.asarray(waves),axis=0)/BINS
phase_probabilities=np.sum(abs(phase_amplitudes)**2,axis=1)
positive=phase_amplitudes[1:BINS//2]
positive_mass=float(sum(phase_probabilities[1:BINS//2]))
rho=positive.T@positive.conj()/positive_mass
lam=.6;theta=math.acos(lam)
uf=A@f;wf=(R@uf-lam*uf)/math.sqrt(1-lam*lam)
plus=(uf-1j*wf)/math.sqrt(2)
positive_fidelity=float(np.vdot(plus,rho@plus).real)
ramsey=[];power=np.eye(30,dtype=complex)
for n in range(33):
    amplitude=np.trace(rho@power)
    ramsey.append({'words':n,'x_plus':float((1+amplitude.real)/2),'y_plus':float((1+amplitude.imag)/2),'amplitude':[float(amplitude.real),float(amplitude.imag)]})
    power=W@power
amplitude=complex(*ramsey[1]['amplitude'])
measured_phase=math.atan2(amplitude.imag,amplitude.real)

# Every eigenphase bin is retained. Positive orientation is a measured outcome.
phase_receipt={'bins':BINS,'probabilities':phase_probabilities.tolist(),
 'amplitudes':encode(phase_amplitudes),'positive_bins':[1,BINS//2-1],
 'positive_probability_given_cooling':positive_mass,
 'complete_positive_source_density':encode(rho),
 'overall_positive_probability':survival*positive_mass,
 'fidelity_to_source_positive_mode':positive_fidelity}

# Actual compound W word, using a factorized density whose equivalence is proved.
pairs=[]
for r in ramsey:
    z=complex(*r['amplitude']);ind=z*z;exclusive=z
    pairs.append({'words':r['words'],
      'independent_x_plus':float((1+ind.real)/2),'independent_y_plus':float((1+ind.imag)/2),
      'exclusive_x_plus':float((1+exclusive.real)/2),'exclusive_y_plus':float((1+exclusive.imag)/2)})
# Direct two-copy vector check, independently of the scalar factorization.
probe=np.kron(plus,plus).reshape(30,30)
joint_after=W@probe@W.T
pair_vector_residual=float(np.linalg.norm(joint_after-np.exp(2j*theta)*probe))

# Coherent source transport. R_e is reconstructed from literal native columns.
work_edges=[(0,1),(0,2),(1,2),(3,4),(3,5),(4,5),(0,3),(1,4),(2,5)]
primitives=[];work_columns=[]
for e,(a0,b0) in enumerate(work_edges):
    for reverse in (False,True):
        edge=(b0,a0) if reverse else (a0,b0)
        perm=[]
        for x in range(6):
            for coin in (0,1):
                q=tuple(int(i==x) for i in range(6));out,(c,j)=native(q,edge,coin)
                y=out.index(1);nc=coin^int(j!=0);perm.append(2*y+nc)
                work_columns.append([e,int(reverse),x,coin,j,y,nc])
        assert all(perm[perm[i]]==i for i in range(12))
        moved=[i for i in range(12) if perm[i]!=i]
        assert len(moved)==2
        primitives.append((edge,perm,moved))
adj={i:[] for i in range(12)}
for i,(edge,perm,moved) in enumerate(primitives):
    a0,b0=moved;adj[a0].append((b0,i));adj[b0].append((a0,i))
word=[];paths=[]
for x,y in ((0,3),(1,4),(2,5)):
    start,goal=2*x,2*y;todo=deque([start]);previous={start:None}
    while goal not in previous:
        p=todo.popleft()
        for q,k in adj[p]:
            if q not in previous:previous[q]=(p,k);todo.append(q)
    path=[];p=goal
    while p!=start:p,k=previous[p];path.append(k)
    path.reverse();palindrome=path+list(reversed(path[:-1]))
    paths.append({'positions':[x,y],'path':path,'word':palindrome});word+=palindrome
work_traces=[];work_permutation=[]
for x in range(6):
    for coin in (0,1):
        state_index=2*x+coin;trace=[]
        for k in word:
            edge=primitives[k][0];p,c=divmod(state_index,2)
            q=tuple(int(i==p) for i in range(6));out,(oldc,j)=native(q,edge,c)
            y=out.index(1);nc=c^int(j!=0)
            trace.append([k,p,c,j,y,nc]);state_index=2*y+nc
        work_permutation.append(state_index)
        work_traces.append({'input':[x,coin],'trace':trace,'output':list(divmod(state_index,2))})
expected=[2*((x+3)%6)+c if c==0 else 2*x+c for x in range(6) for c in (0,1)]
assert work_permutation==expected

# Full binary zero/two-unit audit: no one-unit transfer claim is borrowed here.
other_sector=[]
inputs=[(0,)*6]
for x in range(6):
    for y in range(x+1,6):inputs.append(tuple(int(i in (x,y)) for i in range(6)))
for q0 in inputs:
    for coin in (0,1):
        q=q0;c=coin;trace=[]
        for k in word:
            out,(oldc,j)=native(q,primitives[k][0],c);nc=c^int(j!=0)
            trace.append([k,c,j,nc,list(out)]);q,c=out,nc
        ideal=q0[3:]+q0[:3]
        other_sector.append({'input':list(q0),'coin':coin,'output':list(q),'output_coin':c,'trace':trace,'is_copy_swap':q==ideal and c==coin,'work_authorized':sum(q0)<=1})

# A unitary input-output test transports the entire unknown internal clock state.
transport=np.zeros((60,60))
for cell in (0,1):
    for p in range(3):
        y=work_permutation[2*(3*cell+p)]//2
        for h in range(D):transport[y*D+h,(3*cell+p)*D+h]=1
clock=np.kron(np.eye(2),W)
loaded=np.zeros((60,60),complex);loaded[:30,:30]=rho
unloaded=transport@loaded@transport.T
before_phase=np.trace(loaded@clock);after_phase=np.trace(unloaded@clock)
work={'edges':work_edges,'native_columns':work_columns,'paths':paths,'primitive_word':word,
 'word_native_events':len(word),'one_unit_traces':work_traces,'zero_two_unit_audit':other_sector,
 'permutation':work_permutation,'operator_commutator':float(np.max(abs(transport@clock-clock@transport))),
 'battery_occupation_before':float(np.trace(loaded[:30,:30]).real),'system_occupation_after':float(np.trace(unloaded[30:,30:]).real),
 'phase_before':[float(before_phase.real),float(before_phase.imag)],'phase_after':[float(after_phase.real),float(after_phase.imag)],
 'phase_energy_transfer_in_511keV_units':1.0,'empty_vacuum_phase':0.0,
 'meaning':'Native transport of a clock-bearing charge between identical empty cells; not calorimetry or inertia.'}

# Same observer/filter under the geometry-preserving always-swap alternative.
control_state=(I-Pcontrol)@initial/2
control_norm=float(control_state@control_state)
for n in range(STEPS):control_state=Pcontrol@control_state
control_survival=float(control_state@control_state)
control_fidelity=float((f@control_state)**2/control_survival)
# Committing all arc coordinates destroys eigenphase interference.
dephased=np.diag(np.diag(rho));dephased_amplitude=np.trace(dephased@W)
# Planted double-tick and opposite sign are sensitivity tests, not emergence.
planted=np.trace(rho@W@W)
controls={'always_swap_transfer':encode(Pcontrol),'always_swap_cooling_success':control_survival,
 'always_swap_compact_fidelity':control_fidelity,'target_compact_fidelity':selected_fidelity,
 'dephased_phase_amplitude':[float(dephased_amplitude.real),float(dephased_amplitude.imag)],
 'double_tick_amplitude':[float(planted.real),float(planted.imag)],
 'negative_phase_probability':float(sum(phase_probabilities[BINS//2+1:])),
 'aliases':{'energy_theta':theta,'energy_theta_plus_2pi':theta+2*math.pi,'all_integer_word_records_identical':True}}

predicted_success=(F(2,25)*F(3,5)**(2*STEPS)+F(3,50)*F(2,5)**(2*STEPS))/4
predicted_fidelity=1/(1+F(3,4)*F(2,3)**(2*STEPS))
measured_pair_phase=math.atan2((amplitude*amplitude).imag,(amplitude*amplitude).real)
e6=-1260*math.log1p(-6/1260);e12=-1260*math.log1p(-12/1260)
gates={'source_operator_residuals':max(source_checks.values())<1e-12,
 'all_selection_outcomes_sum_one':abs(selection_total-1)<1e-12,
 'cooling_prediction':abs(survival-float(predicted_success))<1e-15,
 'compact_selection':selected_fidelity>.998,
 'positive_mode_fidelity':positive_fidelity>.99,
 'phase_after_selection':abs(measured_phase-theta)<.01,
 'pair_phase_additive':abs(measured_pair_phase-2*measured_phase)<1e-12,
 'direct_compound_word':pair_vector_residual<1e-12,
 'actual_copy_exchange':work['operator_commutator']<1e-12 and abs(after_phase-before_phase)<1e-12,
 'changed_law_fails_selection':control_fidelity<.8,
 'physical_rest_energy_claim':False}
output={'candidate':'A-native-triangle-vacuum-work-v1','source_records':records,'negative_source_records':negative_records,
 'changed_law_records':control_records,'P':encode(P),'W':encode(W),'source_checks':source_checks,
 'preparation':'One supplied unit at a leaf; no compact ray or spectral projector supplied.',
 'birth_outcomes':first,'filter_first_failure_outcomes':failures,'selection_success':survival,
 'selection_exact_probability':str(predicted_success),'selection_compact_fidelity':selected_fidelity,
 'selection_exact_fidelity':str(predicted_fidelity),'selected_position_state':psi.tolist(),
 'phase':phase_receipt,'ramsey':ramsey,'pairs':pairs,'pair_vector_residual':pair_vector_residual,
 'clock_theta':theta,'measured_one_word_phase':measured_phase,'measured_pair_phase':measured_pair_phase,
 'anchored_pair_keV':511*measured_pair_phase/measured_phase,
 'second_band_predicted_keV':511*math.acos(.4)/theta,
 'native_events_literal_function_calls':CALLS,'cooling_represented_attempts_per_success_path':STEPS+1,
 'phase_controlled_W_requests':BINS-1,'composed_record_laws_not_realized_samples':True,
 'work':work,'controls':controls,'old_511keV_pair_excess_keV':511*(e12/e6-2),'gates':gates,
 'physical_scope':'Source-selected localized operational quasienergy, same-source phase transfer and bounded pair composition. Integer-clock alias remains. No physical electron mass award or full G6.'}
print(json.dumps(output,indent=2,sort_keys=True))
assert all(v for k,v in gates.items() if k!='physical_rest_energy_claim')
