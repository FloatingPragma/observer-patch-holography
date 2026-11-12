"""All-Fock phase selection and a fourth-order coherent rest-mass compiler."""
from fractions import Fraction as F
from itertools import product
from math import factorial
import mpmath as mp
import numpy as np
from model.lepton_mass_compiler import analyzer, compiler as source


def _add(a,b):
    n=max(len(a),len(b));return tuple((a[i] if i<len(a) else F(0))+(b[i] if i<len(b) else F(0)) for i in range(n))


def _mul(a,b):
    c=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return tuple(c)


def _scale(a,x):return tuple(v*x for v in a)


def _remainder(a):
    a=list(a)
    # p=1/(4-cuberoot(4)), so 60p^3-48p^2+12p-1=0.
    for d in range(len(a)-1,2,-1):
        x=a[d];a[d]=F(0)
        a[d-1]+=x*F(4,5);a[d-2]-=x*F(1,5);a[d-3]+=x*F(1,60)
    return tuple(a[:3])


def formal_order_certificate():
    """Exact noncommutative word coefficients, without matrix spot checks."""
    answer={(): (F(1),)}
    palindrome=((0,F(1,8)),(1,F(1,8)),(2,F(1,8)),(3,F(1,4)),(2,F(1,8)),(1,F(1,8)),(0,F(1,8)))
    for duration in ((F(0),F(1)),(F(0),F(1)),(F(1),F(-4)),(F(0),F(1)),(F(0),F(1))):
        for letter,weight in palindrome:
            coefficient=_scale(duration,weight);powers=[(F(1),)]
            for degree in range(1,5):powers.append(_scale(_mul(powers[-1],coefficient),F(1,degree)))
            updated={}
            for word,a in answer.items():
                for degree in range(5-len(word)):
                    key=word+(letter,)*degree
                    updated[key]=_add(updated.get(key,(F(0),)),_mul(a,powers[degree]))
            answer=updated
    wrong=0
    for degree in range(5):
        target=F(1,4**degree*factorial(degree))
        for word in product(range(4),repeat=degree):
            delta=_add(answer.get(word,(F(0),)),(-target,))
            if any(_remainder(delta)):raise ArithmeticError('fourth-order identity failed')
            # p=2/5 violates the cubic cancellation; this must change some words.
            if sum(v*F(2,5)**j for j,v in enumerate(delta)):wrong+=1
    if wrong==0:raise ArithmeticError('altered composition coefficient was not detected')
    return {'exact_noncommuting_words':sum(4**d for d in range(5)),
            'orders_checked':list(range(5)),'p_polynomial':'60p^3-48p^2+12p-1',
            'wrong_p_two_fifths_nonzero_word_errors':wrong}


def trace_error_bound(rounds,time=F(1),sites=1):
    if type(rounds) is not int or type(sites) is not int or rounds<1 or sites<1:
        raise ValueError('positive integer rounds and sites required')
    if type(time) not in (int,F) or time<0:raise ValueError('nonnegative rational time required')
    t=F(time)
    if rounds<14*sites*t:raise ValueError('Taylor radius exceeds one')
    return min(F(2),F((14*sites)**5,10)*t**5/rounds**4)


def rounds_for_error(epsilon=F(1,10**12),time=F(1),sites=1):
    if type(epsilon) not in (int,F) or epsilon<=0:raise ValueError('positive exact error required')
    if type(sites) is not int or sites<1 or type(time) not in (int,F) or time<0:
        raise ValueError('positive integer sites and nonnegative rational time required')
    radius=14*sites*F(time)
    low=max(1,(radius.numerator+radius.denominator-1)//radius.denominator);high=low
    while trace_error_bound(high,time,sites)>epsilon:high*=2
    while low<high:
        middle=(low+high)//2
        if trace_error_bound(middle,time,sites)<=epsilon:high=middle
        else:low=middle+1
    return low


def active_matrices(effect):
    """Literal four-mode CAR reduction of the source's rank-one W block."""
    e=float(effect);w=np.array([[e,np.sqrt(e*(1-e))],[np.sqrt(e*(1-e)),1-e]])
    c=[source._annihilator(4,i) for i in range(4)]
    v=sum((w[i,j]*c[i].conj().T@c[2+j].conj().T for i in range(2) for j in range(2)),np.zeros((16,16),complex))
    v=v+v.conj().T
    count=np.array([((mask>>1)&1)+((mask>>3)&1) for mask in range(16)])
    rotations=[np.diag((1j)**(r*count)) for r in range(4)]
    hs=[r.conj().T@v@r for r in rotations]
    physical=[0,1,4,5];j=np.eye(16)[:,physical]
    avg=sum(hs)/4;target=j.conj().T@v@j
    if np.linalg.norm(avg@j-j@target)>1e-14:raise ArithmeticError('phase cycle did not select the exact compression')
    bad=(hs[0]+hs[2])/2
    wrong=np.linalg.norm(bad@j-j@target)
    if wrong<1e-3:raise ArithmeticError('two-phase mutation failed')
    if np.linalg.norm(v@v@v-v)>1e-14:raise ArithmeticError('exact one-pair exponent identity failed')
    return v,hs,j,target,wrong


def high_precision_cycle(effect,rounds,precision=65):
    """Execute every active Fock column; retained ancilla, no resets or postselection."""
    with mp.workdps(precision):
        e=mp.mpf(effect.numerator)/effect.denominator if isinstance(effect,F) else mp.mpf(str(effect));w=mp.matrix([[e,mp.sqrt(e*(1-e))],[mp.sqrt(e*(1-e)),1-e]])
        c=[mp.matrix([[mp.mpc(z.real,z.imag) for z in row] for row in source._annihilator(4,i)]) for i in range(4)]
        v=mp.zeros(16)
        for i in range(2):
            for j in range(2):v+=w[i,j]*(c[i].H*c[2+j].H)
        v+=v.H
        count=[((mask>>1)&1)+((mask>>3)&1) for mask in range(16)]
        rs=[mp.diag([mp.j**(r*x) for x in count]) for r in range(4)]
        hs=[r.H*v*r for r in rs]
        physical=[0,1,4,5];j=mp.matrix([[int(row==col) for col in physical] for row in range(16)])
        target=j.H*v*j
        p=1/(4-mp.root(4,3));dt=mp.mpf(1)/rounds
        palindrome=((0,mp.mpf(1)/8),(1,mp.mpf(1)/8),(2,mp.mpf(1)/8),(3,mp.mpf(1)/4),(2,mp.mpf(1)/8),(1,mp.mpf(1)/8),(0,mp.mpf(1)/8))
        cache={}
        def step(duration):
            key=mp.nstr(duration,precision)
            if key in cache:return cache[key]
            out=mp.eye(16)
            for r,weight in palindrome:
                x=duration*weight
                pulse=mp.eye(16)-mp.j*mp.sin(x)*hs[r]+(mp.cos(x)-1)*(hs[r]*hs[r])
                out=pulse*out
            cache[key]=out
            return out
        positive=step(p*dt);negative=step((1-4*p)*dt)
        one=positive*positive*negative*positive*positive
        actual=(one**rounds)*j;ideal=j*mp.expm(-mp.j*target)
        error=mp.sqrt(sum(abs(z)**2 for z in actual-ideal))
        gram_error=mp.sqrt(sum(abs(z)**2 for z in actual.H*actual-mp.eye(4)))
        # Maximally entangled input tests the entire physical/reference map.
        outcome=[]
        for mask in (0,2,8,10):
            outcome.append(sum(abs(actual[word|mask,col])**2 for col in range(4) for word in physical)/4)
        if abs(sum(outcome)-1)>mp.mpf('1e-50'):raise ArithmeticError('complete retained probabilities fail')
        bound=trace_error_bound(rounds)
        if 2*error>mp.mpf(bound.numerator)/bound.denominator:raise ArithmeticError('executed error exceeds proved all-Fock bound')
        ground=mp.matrix([1/mp.sqrt(2),0,0,-1/mp.sqrt(2)])
        particle=mp.matrix([0,1,0,0])
        ag=(j*ground).H*(actual*ground);ap=(j*particle).H*(actual*particle)
        contrast=mp.conj(ag[0])*ap[0]
        probabilities=[(1+contrast.real)/4,(1-contrast.real)/4,(1-contrast.imag)/4,(1+contrast.imag)/4]
        measured=mp.atan2(-contrast.imag,contrast.real)
        if abs(measured-e)>4*mp.mpf(bound.numerator)/bound.denominator:
            raise ArithmeticError('complete Ramsey mass bound failed')
        return {'rounds':rounds,'complete_source_dimension':16,'physical_input_columns':4,
                'untouched_reference_dimension':4,'native_Hamiltonian_segments':35*rounds,
                'physical_embedding_Frobenius_error':mp.nstr(error,20),
                'isometry_Gram_error':mp.nstr(gram_error,12),
                'complete_final_ancillary_probabilities':[mp.nstr(x,25) for x in outcome],
                'nonvacuum_probability':mp.nstr(sum(outcome[1:]),20),
                'complete_unhalved_trace_bound_exact':str(bound),
                'arbitrary_reference_bound':True,'vacuum_postselection':False,
                'floating_execution_is_interval_certified':False,
                'complete_Ramsey_Xplus_Xminus_Yplus_Yminus':[mp.nstr(x,25) for x in probabilities],
                'executed_Ramsey_energy':mp.nstr(measured,28),
                'Ramsey_energy_error_diagnostic':mp.nstr(abs(measured-e),15),
                'Ramsey_reference':'idle compensating charge/parity reference; fixed ground/positive-charge eigenspaces; four effects retain all complements'}


def run():
    algebra=formal_order_certificate()
    actual=analyzer.certificate()
    intervals=analyzer.rational_effect_intervals()['effect_intervals_k0_k1_k2']
    effects=actual['effect_values_k0_k1_k2'];n=rounds_for_error()
    rows=[]
    for k,e in enumerate(effects):
        v,hs,j,target,wrong=active_matrices(e)
        midpoint=sum(map(F,intervals[k]),F(0))/2
        row=high_precision_cycle(midpoint,n)
        row.update({'character':k,'native_effect':e,'native_effect_interval':intervals[k],
                    'two_phase_unwanted_ancillary_pair_norm':float(wrong)})
        rows.append(row)
    return {'schema':'oph.openwave.lepton_phase_cycle.v1',
            'exact_formal_product_certificate':algebra,
            'native_analyzer_all_columns':actual['all_inverse_columns_checked'],
            'source_phase_exact':actual['phase_exact'],'rows':rows,
            'all_Fock_bound_domain':'90 physical modes and 324000 enlarged modes per rest site, all occupations and untouched references; ||V||=6 per site',
            'old_reset_rounds_for_same_trace_budget':108*10**12,
            'coherent_cycle_rounds_for_same_trace_budget':n,
            'native_Hamiltonian_segments_per_cycle':35,
            'single_initial_ancillary_vacuum':True,'intermediate_measurements_or_resets':False,
            'finite_precision_source_masses_unchanged':True,
            'electron_anchored_source_values':source.electron_calibration(),
            'controlled_mass_intervals':mass_intervals(n),
            'negative_time_control':'Conjugate V by parity of its complete left lepton-doublet leg; all source evolution durations are positive.',
            'scope':'Declared coherent phase-cycle service for the native lepton rest pairing at H0=0. Exact source coefficients and retrospective phase ancestry are unchanged. Negative interaction segments use a charged-leg sign conjugation; total source-pulse duration is (8p-1)t. This is not the former reset protocol, autonomous evolution, a physical clock, selected Higgs vacuum or an interacting lepton identification.',
            'verdict':'PASS_SCOPED_COHERENT_REST_COMPILER'}


def mass_intervals(rounds):
    """Reuse the source Ramsey error lemma, with the coherent compiler bound."""
    import json
    from pathlib import Path
    root=Path(__file__).resolve().parents[2]
    data=json.loads((root/'model/lepton_cycle/experimental_inputs.json').read_text())
    anchor=F(data['electron_MeV']['value'])/1000
    boxes=[tuple(map(F,x)) for x in analyzer.rational_effect_intervals()['effect_intervals_k0_k1_k2']]
    epsilon=trace_error_bound(rounds);delta=4*epsilon
    assert delta<F(1,20000)
    rows=[]
    for name,k in [('muon',2),('tau',0)]:
        ratio=((boxes[k][0]-delta)/(boxes[1][1]+delta),(boxes[k][1]+delta)/(boxes[1][0]-delta))
        rows.append({'label':name,'ratio_exact':[str(x) for x in ratio],
                     'GeV_exact':[str(anchor*x) for x in ratio],
                     'GeV_decimal':[format(float(anchor*x),'.15g') for x in ratio]})
    return {'source_energy_error_exact':str(delta),'rows':rows,
            'electron_anchor_GeV_exact':str(anchor),'empirical_anchor_uncertainty_included':False,
            'finite_shot_or_physical_model_error_included':False,
            'meaning':'Complete ideal record-law bound, with source coefficients exact; floating active-span execution is a separate numerical check.'}
