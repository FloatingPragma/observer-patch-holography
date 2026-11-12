"""Source partial-isometry algebra, complete bins and autonomous beat bounds."""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json

import sympy as s

from model.openwave_energy_resolved_transport.response import (native_pole_residue,source_bounds,
    encoded_interval,upper,lower)
from model.coherent_defect_wave.chiral_joint.interacting_charged_response.model import operator_certificate
from model.coherent_defect_wave.chiral_joint.interacting_charged_response.transfer import data as charged_data
from model.certified_reference.interval import sine_cosine,pi,SCALE
from model.particle_common_contract.provenance import import_closure

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def partial_isometry_algebra():
    # Two copies of a native input subspace and one retained spectator.
    # These exact identities are dimension independent after tensoring by
    # the source P0 and its unchanged environment, since A*A=P and A^2=0.
    P=s.diag(1,1,0,0,0);Q=s.diag(0,0,1,1,0)
    A=s.zeros(5);A[2,0]=A[3,1]=1
    X=A+A.H;Y=s.I*(A-A.H);Z=P-Q;K=P+Q
    assert A.H*A==P and A*A.H==Q and A*A==s.zeros(5)
    assert X*X==Y*Y==Z*Z==K and X*Y-Y*X==2*s.I*Z
    preparation=s.eye(5)-K+(K+A-A.H)/s.sqrt(2)
    assert preparation.H*preparation==s.eye(5)
    assert preparation*P==(P+A)/s.sqrt(2)
    en,omega,other=s.symbols('en omega other',real=True)
    H=s.diag(en,en,en+omega,en+omega,other)
    assert H*A-A*H==omega*A
    assert s.I*(H*X-X*H)==omega*Y
    assert (s.I*(H*Y-Y*H)+omega*X).applyfunc(s.simplify)==s.zeros(5)
    assert (P*(preparation.H*H*preparation-H)*P-omega*P/2).applyfunc(s.simplify)==s.zeros(5)
    effects=[(s.eye(5)+X)/4,(s.eye(5)-X)/4,(s.eye(5)+Y)/4,(s.eye(5)-Y)/4]
    assert sum(effects,s.zeros(5))==s.eye(5)
    assert all(set(effect.eigenvals()).issubset({0,s.Rational(1,4),s.Rational(1,2)}) for effect in effects)
    number=s.diag(0,0,1,1,2)
    assert number*A-A*number==A
    assert s.I*(H*P-P*H)==s.zeros(5)
    return {'identities':'A*=adjoint; A*A=P, AA*=Q, A^2=0, X^2=Y^2=Z^2=P+Q',
            'actual_source_A':'S_a B_ec (P0 tensor reference_vacuum)/sqrt(w)',
            'source_isometry_premise':'P0 B_ec* S_a B_ec P0=w P0, from the actual native pole residue',
            'fixed_observables':'X=A+A*, Y=i(A-A*), zero on the retained orthogonal complement',
            'single_initial_unitary':'I-P-Q+(P+Q+A-A*)/sqrt(2)',
            'atomic_Heisenberg_relations':['i[H_atomic,X]=Omega Y','i[H_atomic,Y]=-Omega X'],
            'initial_preparation_mean_energy_increase':'Omega/2, including the reference contribution',
            'ongoing_phase_rotations_or_periodic_drive':False,
            'initial_coherent_service_and_terminal_measurement_are_declared':True,
            'passive_stationary_input_self_starts':False}


def interchannel_algebra():
    P=s.diag(1,1,0,0,0,0,0)
    AL=s.zeros(7);AL[2,0]=AL[3,1]=1
    AE=s.zeros(7);AE[4,0]=AE[5,1]=1
    QL=AL*AL.H;QE=AE*AE.H
    assert AL.H*AL==AE.H*AE==P and AE.H*AL==s.zeros(7)
    D=AE*AL.H
    assert D.H*D==QL and D*D.H==QE and D*D==s.zeros(7)
    mu,en,dl,de,other=s.symbols('mu en dl de other',real=True)
    N=s.diag(0,0,1,1,1,1,2)
    H=s.diag(en,en,en+mu+dl,en+mu+dl,en+mu+de,en+mu+de,other)
    assert H*D-D*H==(de-dl)*D and N*D==D*N
    initial_swap=s.eye(7)-P-QL+AL-AL.H
    mix=s.eye(7)-QL-QE+(QL+QE+D-D.H)/s.sqrt(2)
    assert initial_swap.H*initial_swap==s.eye(7) and mix.H*mix==s.eye(7)
    assert mix*initial_swap*P==(AL+AE)/s.sqrt(2)
    X=D+D.H;Y=s.I*(D-D.H)
    assert QL*X*QL==QE*X*QE==s.zeros(7)
    assert QL*Y*QL==QE*Y*QE==s.zeros(7)
    assert (P*((mix*initial_swap).H*H*(mix*initial_swap)-H)*P-(mu+(de+dl)/2)*P).applyfunc(s.simplify)==s.zeros(7)
    return {'actual_partial_isometries':'A_r=S_r B_r P0/sqrt(w_r), r=e^c,L, with full reference-family labels retained',
            'orthogonality':'A_ec* A_L=0 from distinct occupied native reference families',
            'transition':'D=A_ec A_L*, D*D=Q_L, DD*=Q_ec, D^2=0',
            'initial_isometry':'(A_L+A_ec)/sqrt(2)',
            'two_declared_initial_unitaries':'swap P0 into Q_L; then quarter-turn mix Q_L and Q_ec',
            'exact_reference_occupation_on_both_arms':'1',
            'exact_reference_cost_cancellation':True,
            'atomic_frequency':'Delta_ec-Delta_L',
            'one_time_preparation_energy_increase':'4+(Delta_ec+Delta_L)/2 relative to P0',
            'full_joint_Gauss_and_even_parity_intertwiners':True,
            'reference_family_dephasing_kills_both_quadratures_for_all_time':True,
            'L_channel_assigned_a_measured_neutrino_or_muon':False}


def trig_box(angle):
    """Integer-only phase reduction and certified sine/cosine intervals."""
    lo,hi=map(F,angle)
    if lo>hi:raise ValueError('ordered angle interval required')
    p=pi();plo,phi=F(p.lo,SCALE),F(p.hi,SCALE);pmid=(plo+phi)/2
    mid=(lo+hi)/2
    turns=(mid+pmid)//(2*pmid)
    offsets=(2*turns*plo,2*turns*phi)
    rlo,rhi=lo-max(offsets),hi-min(offsets)
    center=(rlo+rhi)/2;radius=(rhi-rlo)/2
    if abs(center)>4:raise ValueError('phase reduction left the certified Taylor domain')
    sine,cosine=sine_cosine(center)
    widen=lambda box:(max(F(-1),F(box.lo,SCALE)-radius),min(F(1),F(box.hi,SCALE)+radius))
    return widen(sine),widen(cosine)


@lru_cache(maxsize=1)
def parameters():
    native_pole_residue();values=source_bounds();data=values['data'];local=data['local']
    ground=local['e_interval'];atom=local['rows'][0]['energy_interval']
    delta=(atom[0]-ground[1],atom[1]-ground[0])
    omega=(delta[0]+4,delta[1]+4)
    neutral=(8*ground[0],8*ground[1]);center=sum(neutral)/2
    radius=s.Rational(1,1024)
    neutral_bin=(center-radius,center+radius)
    charged_bin=data['window']
    d0=radius-(neutral[1]-neutral[0])/2
    d1=data['distances'][0]
    eta=upper(data['v']/min(d0,d1))
    assert eta<s.Rational(1,1000)
    frequency=(charged_bin[0]-neutral_bin[1],charged_bin[1]-neutral_bin[0])
    return {'source_gap':delta,'omega':omega,'v':data['v'],'neutral_bin':neutral_bin,
            'charged_bin':charged_bin,'frequency_support':frequency,
            'herald_vector_error_upper':eta,'herald_success_lower':lower(1-eta**2),
            'conditional_herald_trace_error_upper':2*eta,
            'physical_preparation_success_lower':values['preparation_success_lower'],
            'P0_success_lower':values['P0_success_lower']}


@lru_cache(maxsize=1)
def interchannel_parameters():
    ec=parameters();lepton=charged_data('L')
    ground=lepton['local']['e_interval'];atom=lepton['local']['rows'][0]['energy_interval']
    dl=(atom[0]-ground[1],atom[1]-ground[0])
    omega=(ec['source_gap'][0]-dl[1],ec['source_gap'][1]-dl[0])
    assert omega[0]>0
    gap=min(lepton['distances'][0],charged_data('ec')['distances'][0])
    eta=upper(ec['v']/gap)
    return {'source_gap':omega,'omega':omega,'v':ec['v'],
            'neutral_bin':lepton['window'],'charged_bin':ec['charged_bin'],
            'frequency_support':(ec['charged_bin'][0]-lepton['window'][1],ec['charged_bin'][1]-lepton['window'][0]),
            'herald_vector_error_upper':eta,'herald_success_lower':lower(1-eta**2),
            'conditional_herald_trace_error_upper':2*eta,
            'physical_preparation_success_lower':ec['physical_preparation_success_lower'],
            'P0_success_lower':ec['P0_success_lower'],
            'L_pump_weight_interval':lepton['local']['rows'][0]['weight_interval']}


def observation_bounds(time,*,state_trace_error=0,generator_error=0,channel='ground_ec'):
    t,state_error,drive=map(s.Rational,(time,state_trace_error,generator_error))
    if t<0 or not 0<=state_error<=2 or drive<0:raise ValueError('nonnegative clock and error controls required')
    if channel not in ('ground_ec','ec_L'):raise ValueError('a proved native transition is required')
    p=parameters() if channel=='ground_ec' else interchannel_parameters()
    angle=(F(p['omega'][0]*t),F(p['omega'][1]*t))
    sine,cosine=trig_box(angle)
    error=min(s.Integer(2),p['conditional_herald_trace_error_upper']+state_error+2*(p['v']+drive)*t)
    boxes=[]
    # Standard Y gives <Y>=-sin(Omega t); the energy phase is atan2(-Y,X).
    for lo,hi in (cosine,(-cosine[1],-cosine[0]),(-sine[1],-sine[0]),sine):
        boxes.append((max(s.Integer(0),(1+s.Rational(lo)-error)/4),
                      min(s.Rational(1,2),(1+s.Rational(hi)+error)/4)))
    return {'source_time':str(t),'channel':channel,'ideal_sine_interval':list(map(str,sine)),
            'ideal_cosine_interval':list(map(str,cosine)),
            'complete_four_conditional_record_intervals':[encoded_interval(row,100) for row in boxes],
            'expectation_error_upper':str(error),'expectation_error_decimal':str(s.N(error,16)),
            'state_trace_error':str(state_error),'extra_generator_norm':str(drive),
            'state_error_domain':'normalized accepted clock state after the common herald',
            'changed_generator_preserves_original_frequency_window_asserted':False}


def complete_records():
    q,a,h=s.symbols('q a h');px=s.symbols('p0:4')
    total=(1-q)+q*(1-a)+q*a*(1-h)+q*a*sum(px)
    assert s.expand(total.subs(px[3],h-sum(px[:3])))==1
    return {'records':[{'record':'physical_ground_preparation_failure','probability':'1-q'},
             {'record':'P0_localization_failure','probability':'q(1-a)'},
             {'record':'common_two_bin_herald_failure','probability':'qa(1-h)'},
             *[{'record':name,'probability':f'qa p{i}'} for i,name in enumerate(('X+','X-','Y+','Y-'))]],
            'definitions':'q=physical prep success; a=P0 success conditional on q; h=common two-bin herald success; sum_i p_i=h',
            'coarse_record_count':7,'exact_probability_sum':'1',
            'single_preparation_unitary_has_no_postselected_pump_branch':True,
            'all_source_copies_reference_occupations_and_environment_complements_retained':True,
            'time_dependent_control_after_initial_preparation':False,
            'terminal_observables_are_fixed_in_time':True}


def source_pins():
    paths={str(path.relative_to(ROOT)) for path in HERE.iterdir() if path.suffix in ('.py','.md','.tex')}
    inherited='model/openwave_energy_resolved_transport/results.json'
    receipt=json.loads((ROOT/inherited).read_text())
    paths.add(inherited)
    for path,expected in receipt['source_sha256'].items():
        if hashlib.sha256((ROOT/path).read_bytes()).hexdigest()!=expected:raise ValueError('changed inherited source: '+path)
        paths.add(path)
    paths.update(import_closure(sorted(paths)))
    return {path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(paths)}


def certificate():
    pins=source_pins();p=parameters();algebra=partial_isometry_algebra()
    inter=interchannel_parameters();inter_algebra=interchannel_algebra()
    pump=operator_certificate()
    assert all(row['full_source_reference_gauge_commutators']==12 and row['full_spin_commutators']==3
               for row in pump['all_five_native_elementary_probe_controls'])
    qah=p['physical_preparation_success_lower']*p['P0_success_lower']*p['herald_success_lower']
    rows=[observation_bounds(t) for t in (0,1,10,1000,10**6)]
    pi_box=pi();pi_upper=s.Rational(pi_box.hi,SCALE)
    cycles=lower(p['omega'][0]*10**6/(2*pi_upper))
    assert cycles>700000
    controls=[]
    for mu in (0,4,7):
        total=(p['source_gap'][0]+mu,p['source_gap'][1]+mu)
        controls.append({'reference_cost_control':str(mu),'total_beat_frequency_interval':encoded_interval(total),
                         'source_gap_after_explicit_reference_subtraction_interval':encoded_interval((total[0]-mu,total[1]-mu)),
                         'changed_reference_cost_is_original_frozen_Hamiltonian':mu==4,
                         'unique_full_reference_vacuum_claimed_for_changed_cost':False})
    result={'schema':'oph.native_free_beat_clock.v1','source_sha256':pins,
        'source':{'Hamiltonian':'unchanged N2 PW119 H0+2^-34 T+4 N_ref',
                  'source_time':'same continuous Hamiltonian parameter used by the native charged readouts',
                  'selected_field':'actual e^c atomic charged pole, complete full-H bins retained',
                  'preparation':'inherited physical ground preparation, complete P0 localization, one declared coherent native-pole unitary, common full-H two-bin herald',
                  'full_reference_and_gauss_constraints_retained':True},
        'native_transition_algebra':algebra,
        'frequency':{'total_beat_interval':encoded_interval(p['omega']),
                     'total_beat_decimal':[str(s.N(x,16)) for x in p['omega']],
                     'source_gap_after_reference_subtraction_interval':encoded_interval(p['source_gap']),
                     'exact_full_H_cross_sector_frequency_support':encoded_interval(p['frequency_support']),
                     'frequency_window_independent_of_initial_P0_spin_state':True,
                     'phase_readout':'atan2(-<Y>,<X>) with the complete four-outcome detector',
                     'phase_unwrapping_or_finite_shot_precision_certified':False},
        'common_herald':{'neutral_reference_vacuum_energy_bin':list(map(str,p['neutral_bin'])),
                         'charged_one_reference_energy_bin':list(map(str,p['charged_bin'])),
                         'one_common_success_flag_preserves_relative_coherence':True,
                         'probability_lower':str(p['herald_success_lower']),
                         'conditional_complete_trace_error_upper':str(p['conditional_herald_trace_error_upper']),
                         'original_input_success_lower':str(lower(qah)),
                         'nonselective_herald_preserves_mean_H_energy':True},
        'time_rows':rows,
        'native_interchannel_beat':{
            'algebra':inter_algebra,
            'frequency_interval':encoded_interval(inter['omega']),
            'frequency_decimal':[str(s.N(x,16)) for x in inter['omega']],
            'exact_full_H_cross_sector_frequency_support':encoded_interval(inter['frequency_support']),
            'L_lowest_pole_pump_weight_interval':encoded_interval(inter['L_pump_weight_interval']),
            'common_herald_success_lower':str(inter['herald_success_lower']),
            'time_rows':[observation_bounds(t,channel='ec_L') for t in (0,1,10,1000,10**6)],
            'common_reference_cost_changes_or_continuing_phase_subtractions_required':False,
            'same_seven_record_law_with_reference_families_coherently_retained':True},
        'disturbance_at_one_million':observation_bounds(10**6,state_trace_error=s.Rational(1,1000),generator_error=p['v']/100),
        'omitted_coherent_preparation_control':'A definite reference-vacuum input has exactly zero X and Y for every time, because H preserves reference occupation.',
        'reference_energy_controls':controls,
        'million_time_atomic_total_beat_cycles_lower':str(cycles),
        'complete_records':complete_records(),
        'scope':'A fixed native observable oscillates during free evolution after one specified coherent preparation. The actual full-H bins give a common frequency window and the source perturbation bounds its waveform. The reference contributes four to the beat energy. This is not passive self-starting from a stationary ground, an autonomous preparation mechanism, frequency locking against changes of H, or a calibrated physical Compton clock.',
        'verdict':'PASS_FREE_NATIVE_BEAT_WITH_FIXED_DETECTOR_AND_COMPLETE_REFERENCE'}
    if source_pins()!=pins:raise ValueError('source changed during native beat replay')
    return result
