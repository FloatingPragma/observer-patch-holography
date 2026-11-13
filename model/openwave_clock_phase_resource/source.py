"""Complete native phase resources and energy-covariant branch controls."""
from fractions import Fraction as F
from itertools import product
import numpy as np
import sympy as s
from model.openwave_passive_clock_records.clock import pair_operators,native_attachment
from model.openwave_autonomous_clock_channel.clock import STEP,PUMP

PROB=STEP*PUMP

def exact_source():
    n,nr,a,_,_=pair_operators()
    N=s.Matrix([[s.Integer(round(z.real)) for z in row] for row in n])
    A=s.Matrix([[s.Integer(round(z.real)) for z in row] for row in a])
    p=s.Rational(PROB.numerator,PROB.denominator)
    K0=s.eye(4)+(s.sqrt(1-p)-1)*(A.H*A)
    Kp=s.sqrt(p)*A;Km=s.zeros(4)
    ks=(K0,Kp,Km)
    assert sum((k.H*k for k in ks),s.zeros(4))==s.eye(4)
    frequencies=(0,1,-2)
    for k,f in zip(ks,frequencies):assert N*k-k*N==f*k
    # Incoming battery energy is one. All three output energies are explicit.
    for k,e in zip(ks,(1,0,3)):assert N*k+e*k==k*(N+s.eye(4))
    return N,A,ks


def complete_word_records(steps=4):
    if type(steps) is not int or not 0<=steps<=7:raise ValueError('bounded complete word inventory required')
    rows=[]
    # This is the exact probability action of the actual four-state source.
    for initial in range(4):
        sub=[]
        for word in product(range(3),repeat=steps):
            state=initial;p=F(1);battery=0
            for event in word:
                battery+=(1,0,3)[event]
                if event==0:
                    if state==0:p*=1-PROB
                elif event==1:
                    if state!=0:p=F(0)
                    else:p*=PROB;state=3
                else:p=F(0)
            energies=(0,0,1,1)
            change=energies[state]-energies[initial]
            if p:assert change+battery-steps==0
            sub.append({'battery_word':list(word),'probability_exact':str(p),
                        'terminal_native_occupation':state,
                        'source_energy_change_over_native_mass_exact':str(change) if p else '0',
                        'total_battery_output_energy_over_native_mass_exact':str(battery),
                        'branch_energy_coherence_exact':'0'})
        assert sum(F(x['probability_exact']) for x in sub)==1
        rows.append({'initial_native_occupation':initial,'all_battery_words':sub,
                     'retained_word_count':3**steps})
    return rows


def battery_basis_control():
    N,A,ks=exact_source();rho=s.diag(1,0,0,0)
    kp,km=(ks[0]+ks[1])/s.sqrt(2),(ks[0]-ks[1])/s.sqrt(2)
    operators=(kp,km,ks[2]);labels=('battery_plus','battery_minus','battery_loss')
    assert sum((k.H*k for k in operators),s.zeros(4))==s.eye(4)
    rows=[]
    for label,k in zip(labels,operators):
        state=s.simplify(k*rho*k.H)
        coherence=s.simplify(s.trace((A+A.H)*state))
        probability=s.simplify(s.trace(state))
        rows.append({'outcome':label,'probability_exact':str(probability),
          'unnormalized_native_X_exact':str(coherence),
          'energy_invariant_branch':N*state-state*N==s.zeros(4)})
    p=s.Rational(PROB.numerator,PROB.denominator)
    assert s.simplify(s.sympify(rows[0]['unnormalized_native_X_exact'])**2-p*(1-p))==0
    assert not rows[0]['energy_invariant_branch']
    return {'records':rows,'measurement_mixes_battery_energies_over_mass':[1,0],
       'source_mean_energy_after_readout_over_mass_exact':str(PROB),
       'battery_mean_energy_before_readout_over_mass_exact':str(1-PROB),
       'battery_mean_energy_after_readout_over_mass_exact':'1/2',
       'readout_mean_work_over_mass_exact':str(PROB-F(1,2)),
       'readout_commutes_with_battery_energy':False,
       'readout_requires_phase_resource':True}


def relational_records():
    N,A,_=exact_source();I=s.eye(4)
    n1=s.kronecker_product(N,I);n2=s.kronecker_product(I,N)
    a1=s.kronecker_product(A,I);a2=s.kronecker_product(I,A)
    H=n1+n2;R=a1*a2.H+a1.H*a2;R2=R*R
    assert H*R-R*H==s.zeros(16)
    assert R2*R==R
    start=s.zeros(16,1);start[12]=1 # |11>_1 |00>_2
    rows=[]
    for label,c,si in [('0',1,0),('pi/4',s.sqrt(2)/2,s.sqrt(2)/2),
                       ('pi/2',0,1),('pi',-1,0)]:
        U=I16=s.eye(16)
        U=I16+(c-1)*R2-s.I*si*R
        assert s.simplify(U.H*U)==I16
        v=s.simplify(U*start);rho=s.simplify(v*v.H)
        assert H*rho-rho*H==s.zeros(16)
        assert s.trace(a1*rho)==s.trace(a2*rho)==0
        prob=[s.simplify(x*s.conjugate(x)) for x in v]
        assert sum(prob)==1
        rows.append({'interaction_phase':label,
          'all_sixteen_native_occupation_probabilities_exact':list(map(str,prob)),
          'relative_coherence_exact':str(s.simplify(s.trace(a1.H*a2*rho))),
          'mean_source_energy_over_native_mass_exact':'1',
          'source_energy_change_exact':'0','switch_on_and_off_mean_work_exact':'0',
          'free_total_energy_commutator_exact':'0'})
    return {'complete_dimension':16,'records':rows,
       'source':'two full four-occupation native pairs and coupling lambda*(A1 A2*+A1* A2)',
       'free_native_frequency_of_relative_coherence_exact':'0',
       'population_Rabi_angular_frequency':'2*lambda',
       'lambda_is_the_native_rest_mass':False,
       'coupling_is_declared_architecture':True}


def certificate():
    N,A,ks=exact_source()
    return {'schema':'oph.native_clock_phase_resource.v1',
       'verdict':'PASS_COMPLETE_ENERGY_COVARIANT_PHASE_AND_RELATIONAL_CONTROLS',
       'native_attachments':native_attachment(),
       'collision_probability_exact':str(PROB),
       'complete_operator_identities':{'four_state_Gram':'identity',
           'battery_total_energy_intertwiner':'(H+Hb)V=V(H+m)',
           'Kraus_source_frequency_over_mass':[0,1,-2],
           'inactive_native_occupations_retained':[1,2]},
       'complete_energy_basis_battery_records':complete_word_records(),
       'energy_incompatible_battery_readout':battery_basis_control(),
       'two_copy_relational_control':relational_records(),
       'general_theorem':{
           'premises':['initial joint state commutes with total energy',
                       'joint unitary conserves total energy',
                       'every recorded battery effect commutes with battery energy'],
           'conclusion':'Every unnormalized conditional source state commutes with its Hamiltonian; every fixed-observable free expectation is constant.',
           'includes_adaptive_energy_covariant_branches':True,
           'degenerate_relative_coherence_is_allowed':True,
           'phaseful_inputs_or_noncommuting_readouts_excluded':True,
           'no_go_for_all_OPH_architectures':False},
       'scope':{'autonomous_rest_mass_phase_onset_established':False,
                'energy_invariant_finite_collision_route_excluded':True,
                'thermal_or_infinite_spontaneous_symmetry_breaking_limit_proved':False}}
