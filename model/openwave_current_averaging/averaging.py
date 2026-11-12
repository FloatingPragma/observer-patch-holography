"""Exact CAR-word resonances and a uniform full-source averaging theorem."""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import comb
import numpy as np
from model.openwave_weak_decay.decay import interaction,mass_data,index,labels
from model.native_flavor_orbit.intervals import I,sqrt_rational

G=F(1,2**72)
GAP=F(1,50000)

def _root(a):return I(sqrt_rational(a.lo).lo,sqrt_rational(a.hi).hi)

def _mono(variable=None,negative=False,rest=False):
    exponents=[0]*6
    if variable is not None:exponents[variable]=1
    return (tuple(exponents),int(rest),F(-1 if negative else 1,2 if rest else 1))

def _multiply(coefficients,base):
    exponent=tuple(sum(v[0][k]for v in coefficients)for k in range(6))
    roots=sum(v[1]for v in coefficients);coefficient=F(base)*2**(roots//2)
    for _,_,c in coefficients:coefficient*=c
    return exponent,roots%2,coefficient

@lru_cache(maxsize=1)
def source():
    m,heavy,boxes=mass_data();q=(1-m)/2;mass=(m,1.,heavy)
    mb=I(*boxes[1])/I(*boxes[2]);hb=I(*boxes[0])/I(*boxes[2]);qb=(1-mb)/2
    energy_boxes=[mb,I(1),hb,*[_root(qb*qb+x*x)for x in(mb,I(1),hb)]]
    energies=[];vectors=[];mapping={};numeric={}
    for i in range(270):
        p,spin,internal,f=labels(i)
        if internal in(9,11):
            E=np.sqrt((q*p)**2+mass[f]**2);v=[0]*6;v[f+(3 if p else 0)]=2
        else:E=abs(p)*q;v=[-abs(p),abs(p),0,0,0,0]
        energies.append(E);vectors.append(tuple(v))
    for i in range(270):
        p,spin,internal,f=labels(i);sgn=1 if spin==0 else -1;epsilon=-sgn*q*p
        if internal==11:
            right=(-p+1)*90+(1-spin)*45+27+f
            u=np.sqrt((energies[i]+epsilon)/(2*energies[i]));v=sgn*np.sqrt((energies[i]-epsilon)/(2*energies[i]))
            numeric[i]=[(i,False,u),(right,True,-v)];numeric[right]=[(right,False,u),(i,True,v)]
            uc=_mono(rest=True)if p==0 else _mono(2*f+(0 if epsilon>0 else 1))
            vc=_mono(negative=sgn<0,rest=True)if p==0 else _mono(2*f+(1 if epsilon>0 else 0),negative=sgn<0)
            minus_v=(vc[0],vc[1],-vc[2])
            mapping[i]=[(i,False,uc),(right,True,minus_v)];mapping[right]=[(right,False,uc),(i,True,vc)]
        elif internal!=9:
            mapping[i]=[(i,epsilon<0,_mono())];numeric[i]=[(i,epsilon<0,1.)]
    assert len(mapping)==len(numeric)==270
    return dict(m=m,heavy=heavy,q=q,energy_boxes=energy_boxes,energies=energies,
                energy_vectors=vectors,mapping=mapping,numeric_mapping=numeric,
                mass_intervals=[mb.strings(),hb.strings()])

@lru_cache(maxsize=1)
def words():
    native,*_=interaction();d=source();all_words=[];resonant=[];frequencies=set()
    for(a,b,c,e),coefficient in native.items():
        choices=[]
        for k,dagger in((a,True),(b,True),(e,False),(c,False)):
            choices.append([(i,creation!=dagger,z)for i,creation,z in d['mapping'][k]])
        for choice in product(*choices):
            operations=tuple((i,creation)for i,creation,_ in choice)
            frequency=tuple(sum((1 if creation else -1)*d['energy_vectors'][i][k]
                                for i,creation in operations)for k in range(6))
            polynomial=_multiply([v[2]for v in choice],coefficient)
            item=(operations,polynomial,frequency);all_words.append(item);frequencies.add(frequency)
            if not any(frequency):resonant.append(item)
    return all_words,resonant,sorted(frequencies),sum(abs(v)for v in native.values())

def act(mask,operations):
    sign=1
    for i,create in reversed(operations):
        if bool((mask>>i)&1)==create:return None,0
        if (mask&((1<<i)-1)).bit_count()%2:sign=-sign
        mask^=1<<i
    return mask,sign

def reduce_polynomial(poly):
    """Canonical remainder modulo y_f^2=1-x_f^2, with sqrt(2) retained exactly."""
    result=defaultdict(F)
    for (powers,root),coefficient in poly.items():
        choices=[]
        for f in range(3):
            x,y=powers[2*f:2*f+2];n=y//2
            choices.append([(x+2*k,y%2,F((-1)**k*comb(n,k)))for k in range(n+1)])
        for term in product(*choices):
            exponent=tuple(v for x,y,_ in term for v in(x,y));factor=coefficient
            for _,_,z in term:factor*=z
            result[exponent,root]+=factor
    return {k:v for k,v in result.items()if v}

def resonant_column(mask):
    output=defaultdict(lambda:defaultdict(F))
    for operations,(powers,root,coefficient),_ in words()[1]:
        target,sign=act(mask,operations)
        if sign:output[target][powers,root]+=sign*coefficient
    reduced={k:reduce_polynomial(v)for k,v in output.items()}
    return {k:v for k,v in reduced.items()if v}

def parent():return 1<<index(0,1,1,1)
def receiver():return (1<<index(0,0,1,0))|(1<<index(1,1,0,1))|(1<<index(1,0,0,0))

@lru_cache(maxsize=1)
def exact_orbit():
    initial=parent();todo=[initial];seen={initial};columns={};constant=((0,)*6,0)
    while todo:
        mask=todo.pop();column=resonant_column(mask);columns[mask]=column
        for target in column:
            if target not in seen:seen.add(target);todo.append(target)
        assert len(seen)<20,'native resonant orbit changed'
    assert seen=={initial,receiver()}
    basis=[initial,receiver()];matrix=[]
    for a in basis:
        row=[]
        for b in basis:
            value=columns[b].get(a,{})
            assert set(value)<= {constant},'nonconstant Bogoliubov coefficient in resonant block'
            row.append(value.get(constant,F(0)))
        matrix.append(row)
    assert matrix==[[F(173,2),F(1)],[F(1),F(85)]]
    return {'basis_masks':list(map(str,basis)),'all_resonant_current_columns':[[str(x)for x in row]for row in matrix],
            'exact_invariant_orbit_dimension':len(basis),'coefficient_identity':'x_f^2+y_f^2=1; all nonconstant coefficients cancel exactly',
            'parent_free_energy':'1','receiver_free_energy':'m+2q=1',
            'parent_full_energy':'1+(173/2)G','receiver_full_energy':'1+85G',
            'resonant_detuning':'3G/2','resonant_off_diagonal':'+G in the canonical quasiparticle Fock phase',
            'original_receiver_probability':'(16/25)sin^2(5Gt/4)',
            'first_peak_scaled_time':'Gt=2pi/5','first_peak_receiver_probability':'16/25',
            'all_other_averaged_source_outcomes':'0 by the exact full-current invariant orbit'}

def frequency_certificate():
    d=source();all_words,resonant,freq,L=words();rows=[];minimum=None
    for v in freq:
        if not any(v):continue
        box=sum((c*energy/2 for c,energy in zip(v,d['energy_boxes'])),I(0))
        assert box.lo>0 or box.hi<0,'unresolved native Bohr-frequency relation'
        lower=min(abs(box.lo),abs(box.hi));minimum=lower if minimum is None else min(minimum,lower)
        rows.append({'doubled_energy_coefficients':list(v),'frequency_interval':box.strings()})
    assert minimum>GAP
    return {'native_quartic_words':len(interaction()[0]),'transformed_ordered_words':len(all_words),
            'zero_frequency_ordered_words':len(resonant),'frequency_labels':len(freq),
            'all_nonzero_frequency_enclosures':rows,'minimum_frequency_lower_diagnostic':float(minimum),
            'certified_uniform_nonzero_frequency_floor':str(GAP),
            'native_current_norm_upper':str(L),'expanded_word_l1_upper':str(4*L),
            'oscillatory_primitive_norm_upper':str(8*L/GAP)}

def full_source_bound():
    f=frequency_certificate();L=F(f['native_current_norm_upper']);primitive=F(f['oscillatory_primitive_norm_upper'])
    # pi<22/7 encloses every scaled time up to the first peak.
    tau=F(44,35);error=G*primitive*(1+2*L*tau)
    assert error<F(1,2500000)
    old=F(1,64)*primitive*(1+2*L*tau);assert old>1
    return {'declared_source_G':str(G),'scaled_time_domain':'0 <= Gt <= 2pi/5',
            'unitary_error_rule':'G*primitive*(1+2*G*L*t)',
            'uniform_full_source_unitary_error_upper':str(error),
            'original_receiver_peak_probability_interval':[str(F(16,25)-2*error),str(F(16,25)+2*error)],
            'original_receiver_peak_probability_decimal_interval':[float(F(16,25)-2*error),float(F(16,25)+2*error)],
            'parent_peak_probability_interval':[str(F(9,25)-2*error),str(F(9,25)+2*error)],
            'all_other_source_outcomes_probability_upper':str(error*error),
            'all_other_source_outcomes_probability_decimal_upper':float(error*error),
            'original_G_1_over_64_comparison_bound':str(old),
            'original_G_1_over_64_averaging_certified':False,
            'final_full_record_energy_work_absolute_upper':str(2*G*L),
            'long_full_fock_evolution_numerically_executed':False,
            'arbitrary_quantum_references':'operator-norm bound is stable under tensoring identity; complete parent/receiver/remainder instrument retained'}

def small_block_records():
    h=np.array([[86.5,1.],[1.,85.]]);eig,U=np.linalg.eigh(h);rows=[]
    for tau in(0.,.2,.5,1.,2*np.pi/5):
        state=U@(np.exp(-1j*tau*eig)*U[0]);prob=abs(state)**2
        analytic=16/25*np.sin(5*tau/4)**2
        assert abs(prob[1]-analytic)<3e-14 and abs(prob.sum()-1)<3e-14
        rows.append({'scaled_time_Gt':tau,'parent_probability':float(prob[0]),'receiver_probability':float(prob[1]),
                     'all_other_records_probability':0.,'analytic_receiver_probability':float(analytic),
                     'current_energy_expectation':float(np.vdot(state,h@state).real),
                     'complete_record_work_divided_by_G':float(np.dot(prob,np.diag(h))-86.5)})
    return rows

def certificate():
    orbit=exact_orbit();frequency=frequency_certificate();bound=full_source_bound();records=small_block_records()
    return {'schema':'oph.full_native_current_uniform_averaging.v1',
            'source':'the complete270-CAR free Hamiltonian and31464 native current-contact words, with the supplied weak coupling varied to2^-72',
            'reference':'same filled negative free bands, zero bands vacant; prepared rest parent and idle opposite-charge reference',
            'all_native_current_spectators_retained':True,'added_occupation_penalty_or_two_state_source':False,
            'mass_intervals':source()['mass_intervals'],'exact_resonant_orbit':orbit,
            'full_source_frequency_certificate':frequency,'uniform_complete_source_comparison':bound,
            'executed_exact_resonant_records':records,
            'zero_coupling_control_receiver_probability':'0',
            'zero_contact_diagonal_correction_control':'Deleting the actual3G/2 detuning changes the resonant maximum from16/25 to1.',
            'scope':'Complete unfiltered native-current averaging on the stated three-momentum box, with a derived invariant resonant orbit and uniform finite-time/reference error. No installed continuum rate, observed weak coupling, physical lifetime, irreversible decay, or neutral-flight source transfer.',
            'verdict':'PASS_UNFILTERED_NATIVE_CURRENT_RESONANT_RECORDS_WITH_UNIFORM_ERROR'}
