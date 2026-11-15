"""Full-CAR compilation, a common kinetic law, and complete mass records."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
import json
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy import linalg as la
import sympy as s

from model.charged_lepton_response import car_koide, relativistic_mass
from model.coherent_defect_wave.mass_instrument import model as apparatus
from model.common_source import native_matter_attachment as matter
from model.common_source import native_record_yukawa_field as field
from model.common_source.chiral_candidate import kinetic
from model.native_flavor_orbit import intervals as iv
from . import analyzer

ROOT=Path(__file__).resolve().parents[2]


def archive_bound(steps,time=F(1),sites=1,per_round_error=F(0)):
    if type(steps) is not int or steps<1 or type(sites) is not int or sites<1:
        raise ValueError('positive integer steps and sites required')
    if type(time) not in (int,F) or type(per_round_error) not in (int,F) or min(time,per_round_error)<0:
        raise ValueError('nonnegative exact clock and implementation error required')
    t=F(time);epsilon=F(per_round_error);leakage=F(36*sites*sites)*t*t/steps
    return {'sites':sites,'steps':steps,'time_exact':str(t),'full_kick_norm_exact':str(6*sites),
            'any_nonvacuum_probability_upper_exact':str(min(F(1),leakage+steps*epsilon)),
            'complete_unhalved_trace_error_exact':str(min(F(2),3*leakage+2*steps*epsilon)),
            'physical_CAR_modes':90*sites,'enlarged_CAR_modes':324000*sites,
            'fresh_ancillary_modes_each_round':323910*sites,
            'total_exported_ancillary_modes':323910*sites*steps,
            'vacuum_postselection':False,'all_quantum_archives_retained':True,
            'bound_record_convention':'vacuum/nonvacuum flags are committed classical records; quantum coherence within every flag sector is retained',
            'state_and_untouched_reference_domain':'every physical Fock state, including a filled negative band',
            'elapsed_native_attempt_clock_identified':False}


def _annihilator(modes,index):
    return np.asarray(apparatus.annihilator(index,modes),complex)


def four_mode_control(effect,steps=4,time=F(1,8),kinetic_cost=F(1,3)):
    """Exact active flag-span reduction, complete numerical sixteen-state CAR."""
    e=float(effect)
    if not 0<e<1:raise ValueError('a strict source effect is required')
    w=np.array([[e,np.sqrt(e*(1-e))],[np.sqrt(e*(1-e)),1-e]])
    c=[_annihilator(4,i) for i in range(4)]
    pair=sum((w[i,j]*(c[i].conj().T@c[2+j].conj().T) for i in range(2) for j in range(2)),np.zeros((16,16),complex))
    interaction=pair+pair.conj().T
    h0=float(kinetic_cost)*sum((a.conj().T@a for a in c),np.zeros((16,16),complex))
    h=h0+interaction;physical=(0,1,4,5);ancillary=(0,2,8,10)
    target=h[np.ix_(physical,physical)]
    physical_c=[_annihilator(2,i) for i in range(2)]
    hp=e*(physical_c[0].conj().T@physical_c[1].conj().T)
    expected=float(kinetic_cost)*sum((a.conj().T@a for a in physical_c),np.zeros((4,4),complex))+hp+hp.conj().T
    if la.norm(target-expected)>2e-14:raise ArithmeticError('occupied-sector compression failed')
    step=float(time)/steps;u=la.expm(-1j*step*h);ks=[]
    for mask in ancillary:
        signs=np.array([(-1)**sum((mask&((1<<p)-1)).bit_count() for p in (0,2) if word>>p&1) for word in physical])
        ks.append(signs[:,None]*u[np.ix_([word|mask for word in physical],physical)])
    gram=sum((k.conj().T@k for k in ks),np.zeros((4,4),complex))
    if la.norm(gram-np.eye(4))>2e-13:raise ArithmeticError('complete ancillary CAR instrument failed')
    states={():np.eye(4,dtype=complex)}
    for _ in range(steps):states={history+(j,):k@v for history,v in states.items() for j,k in enumerate(ks)}
    probabilities={','.join(map(str,key)):float(la.norm(value,'fro')**2/4) for key,value in states.items()}
    coarse={','.join(map(str,key)):0. for key in product((0,1),repeat=steps)}
    coarse_kets={key:[] for key in coarse}
    for history,operator in states.items():
        label=','.join(str(int(j!=0)) for j in history)
        coarse[label]+=float(la.norm(operator,'fro')**2/4)
        coarse_kets[label].append(operator.reshape(-1)/2)
    quantum_gram_error=0.;retained_coherence=0.
    for label,blocks in coarse_kets.items():
        ket=np.concatenate(blocks)
        quantum_gram_error=max(quantum_gram_error,abs(float(np.vdot(ket,ket).real)-coarse[label]))
        weights=np.array([float(np.vdot(v,v).real) for v in blocks])
        retained_coherence=max(retained_coherence,float(np.sqrt(max(0.,weights.sum()**2-np.dot(weights,weights)))))
    zero=states[(0,)*steps];ideal=la.expm(-1j*float(time)*target)
    v=zero.reshape(-1)/2;z=ideal.reshape(-1)/2;weight=float(np.vdot(v,v).real)
    error=float(np.sqrt(max(0.,(1+weight)**2-4*abs(np.vdot(z,v))**2))+1-weight)
    bound=3*float(time)**2/steps
    if abs(sum(probabilities.values())-1)>2e-12 or error>bound+2e-12:
        raise ArithmeticError('complete reference-safe active-span comparison failed')
    return {'effect':e,'full_Fock_dimension':16,'physical_input_dimension':4,
            'every_physical_occupation_included':True,'all_ancillary_outcomes_each_round':4,
            'steps':steps,'complete_history_count':len(probabilities),'complete_original_history_probabilities':probabilities,
            'coarse_vacuum_nonvacuum_history_probabilities':coarse,
            'retained_coarse_quantum_ket_probability_error':quantum_gram_error,
            'maximum_retained_fine_label_offdiagonal_density_Frobenius_norm':retained_coherence,
            'fine_labels_are_archive_diagonal_coordinates_not_production_measurements':True,
            'quantum_state_within_each_coarse_history':'coherent sum of its Kraus operators tensor the corresponding orthogonal ancillary history kets',
            'reference':'maximally entangled physical/reference state; quantum complements retained',
            'kinetic_cost_exact':str(kinetic_cost),'common_H0_commutes_analyzer_but_not_pairing':True,
            'compression_error':float(la.norm(target-expected)),
            'complete_Kraus_Gram_error':float(la.norm(gram-np.eye(4))),
            'complete_archive_unhalved_error':error,'bound':bound,
            'any_nonvacuum_probability':1-weight,'comparison_is_finite_N_spectral_equality':False}


@lru_cache(maxsize=1)
def full_native_dirac():
    e=np.array(analyzer.certificate()['effect_values_k0_k1_k2'])
    rows=[r for r in matter.invariant_higgs_lines() if r['neutral_rank']==2]
    if len(rows)!=1:raise ArithmeticError('unique directed native charged Higgs line required')
    rep=matter.representation();indices=list(rep['selected'])
    weak=rep['operators']['weak_C2'].extract(indices,indices)
    left=s.diag(*(int(weak[i,i]!=0) for i in range(15)))
    d=left*rows[0]['components'][1]*(s.eye(15)-left)
    if s.trace(d.H*d)!=1 or d*d!=s.zeros(15) or d.rank()!=1:raise ArithmeticError('native charged normalization differs')
    dn=np.asarray(d,complex);directed=np.kron(dn,np.diag(e));kernel=directed+directed.T
    expected_masses=np.sort(np.r_[np.zeros(39),np.repeat(e,2)])
    if np.max(abs(np.sort(la.svdvals(kernel))-expected_masses))>2e-14:
        raise ArithmeticError('the native lift changed the three source masses or its zero spectators')
    controls=[];p=(.17,-.23,.31)
    for valley in product((0,1),repeat=3):
        momentum=tuple((-1)**b*x for b,x in zip(valley,p,strict=True))
        report=field.dispersion_certificate(kernel,momentum,speed=1)
        if report['dispersion_spectral_error']>2e-12 or report['Hermiticity_error']>2e-12:
            raise ArithmeticError('full native Dirac rows differ')
        controls.append({'valley':list(valley),'complete_Nambu_dimension':180,
            'spectral_error':report['dispersion_spectral_error'],'positive_charged_rows':12,
            'complete_180_energy_rows':la.eigvalsh(kinetic.bdg_symbol(momentum,kernel,speed=1)).tolist(),
            'spectator_internal_zero_coefficients':39})
    return {'current_support':car_koide.higgs_covariance(),'Clifford':relativistic_mass.clifford_certificate(),
            'all_eight_valleys':controls,'fixed_scalar_background_norm_exact':'1',
            'background_convention':'w_H=1 is the declared neutral background norm; under the usual parametrization w_H=v/sqrt(2), v=sqrt(2) in source units; no empirical or selected interacting VEV',
            'canonical_physical_Higgs_kinetic_residue_derived':False,
            'compiled_Y_identified_with_absolute_Standard_Model_Yukawa':False,
            'same_effect_rest_energies_k0_k1_k2':e.tolist(),
            'complete_45_internal_rest_coefficients':expected_masses.tolist(),
            'mass_rule':'m_j=w_H e_j; no second squaring',
            'common_kinetic':'identical spin Weyl hopping on every internal/family/flag copy',
            'extended_kinetic_vacuum_reducing':True,'extended_kinetic_commutes_full_analyzer':True,
            'ancillary_vacuum_is_extended_ground':False,
            'same_clock_continuum_error':relativistic_mass.continuum_bound('1/64',1,1,1,'1/20000')}


@lru_cache(maxsize=1)
def exact_rest_instrument():
    """All36 records of one Dirac momentum block, with symbolic energy m>0."""
    beta=s.kronecker_product(s.Matrix([[0,1],[1,0]]),s.eye(2));negative=(s.eye(4)-beta)/2
    c=[apparatus.annihilator(i) for i in range(4)]
    b=(c[0],c[1],c[3].H,-c[2].H)
    native_pair=c[0].H*c[3].H-c[1].H*c[2].H
    native_pair+=native_pair.H
    native_electric=-c[0].H*c[0]-c[1].H*c[1]+c[2].H*c[2]+c[3].H*c[3]
    transform=s.zeros(16)
    for mask in range(16):
        column=s.eye(16)[:,12]
        for index in reversed(range(4)):
            if mask>>index&1:column=b[index].H*column
        transform[:,mask]=column
    if transform.H*transform!=s.eye(16) or any(transform.H*b[i]*transform!=c[i] for i in range(4)):
        raise ArithmeticError('literal all-left to charged-Dirac CAR map failed')
    h=sum((beta[i,j]*c[i].H*c[j] for i in range(4) for j in range(4)),s.zeros(16))+2*s.eye(16)
    relative_number=sum((a.H*a for a in c),s.zeros(16))-2*s.eye(16)
    if transform.H*native_pair*transform!=h-2*s.eye(16) or transform.H*native_electric*transform!=-relative_number:
        raise ArithmeticError('same pairing energy or electric charge changed under CAR map')
    rho=s.zeros(16)
    for a in combinations(range(4),2):
        for b in combinations(range(4),2):rho[sum(1<<i for i in a),sum(1<<i for i in b)]=negative.extract(a,b).det()
    if rho*rho!=rho or s.trace(rho)!=1 or h*rho!=s.zeros(16):raise ArithmeticError('same-mass Dirac sea failed')
    energies=apparatus.energy_projectors(h,s.Integer(1));joint=[]
    for j,p in enumerate(energies):
        for charge in range(-2,3):
            q=s.diag(*(int(mask.bit_count()-2==charge) for mask in range(16)));effect=p*q
            if effect!=s.zeros(16):joint.append((j,charge,effect))
    if sum((p for _,_,p in joint),s.zeros(16))!=s.eye(16) or len(joint)!=9:raise ArithmeticError('complete spectral effects differ')
    records=[]
    for initial in (1,0):
        branches=apparatus.reservoir_branches(rho,initial)
        for flag,state in enumerate(branches):
            for j,q,effect in joint:
                probability=s.trace(effect*state)/2
                if probability.is_nonnegative is not True:raise ArithmeticError('negative original probability')
                if probability!=0 and -q-flag+initial!=0:
                    raise ArithmeticError('complete source/reservoir/dual-reference electric charge changed')
                records.append({'setting':'add' if initial else 'remove','reservoir_out':flag,
                    'success':flag!=initial,'relative_charge':q,'energy_in_units_of_same_m':j,
                    'native_electric_charge':-q,
                    'idle_dual_charge_reference_occupation':initial,
                    'probability_exact':str(probability)})
    if sum(s.Rational(r['probability_exact']) for r in records)!=1:raise ArithmeticError('mass readout lost weight')
    return {'source_Fock_states':16,'with_reservoir_states':32,'complete_records':records,
            'exact_all_left_to_Dirac_annihilator_map':['l_up','l_down','r_down_conjugate†','-r_up_conjugate†'],
            'full_16_state_CAR_intertwiner_energy_and_electric_charge_checked':True,
            'relative_charge_convention':'Dirac number relative to the filled sea; native electric charge is its negative',
            'probe_reference':'reservoir has electric charge -1; idle dual-charge reference starts with the same occupation as the reservoir, so every nonzero complete record has zero total electric charge and even total parity',
            'record_count':36,'zero_probability_records':sum(r['probability_exact']=='0' for r in records),
            'source_reference':'filled negative band of same limiting Dirac H; other modes retain their sea',
            'all_three_families_and_all_eight_valleys_use_same_36_record_rule':True,
            'failure_can_have_energy_twice_m':True,'ideal_spectral_projectors_are_declared_detector_service':True}


def ramsey_readout(steps=108000000000000,time=F(1),sites=1):
    bound=archive_bound(steps,time,sites);epsilon=F(bound['complete_unhalved_trace_error_exact'])
    if epsilon>F(1,4) or time<=0:raise ValueError('resolved complete-output error and positive time required')
    delta=4*epsilon/F(time)
    effect_boxes=[iv.I(*map(F,p)) for p in analyzer.rational_effect_intervals()['effect_intervals_k0_k1_k2']]
    if delta>=F(1,20000):raise ValueError('absolute mass error does not resolve the electron denominator')
    rows=[]
    for k,box in enumerate(effect_boxes):
        angle=box*F(time)
        if angle.hi>=1:raise ValueError('fixed small-angle mass readout branch required')
        sine=iv.sine(angle);cosine=1-2*iv.sine(angle/2)**2
        probabilities=[iv.rounded_out(p,38,probability=True) for p in
                       ((1+cosine)/4,(1-cosine)/4,(1+sine)/4,(1-sine)/4)]
        rows.append({'character':k,'ideal_Xplus_Xminus_Yplus_Yminus_intervals':[p.strings() for p in probabilities],
            'compiled_original_probability_intervals':[[str(max(F(0),p.lo-epsilon/2)),str(min(F(1),p.hi+epsilon/2))] for p in probabilities],
            'energy_estimate_interval':[str(box.lo-delta),str(box.hi+delta)]})
    electron=effect_boxes[1];ratios=[]
    empirical=json.loads((ROOT/'model/lepton_cycle/experimental_inputs.json').read_text())
    anchor=F(empirical['electron_MeV']['value'])/1000
    for name,k in (('muon',2),('tau',0)):
        numerator=effect_boxes[k]
        limits=((numerator.lo-delta)/(electron.hi+delta),(numerator.hi+delta)/(electron.lo-delta))
        ratios.append({'label':name,'ratio_interval':[str(x) for x in limits],
                      'electron_anchored_GeV_interval':[str(anchor*x) for x in limits]})
    return {'compiler':bound,'source_time_exact':str(time),'setting_probabilities':['1/2','1/2'],
            'compensated_reference':'coherent source ground/reference blank plus positive-energy charged excitation/dual-charge reference occupied, normalized by sqrt2; total parity even and total electric charge zero',
            'complete_effects':'(I±X)/4 and (I±Y)/4; X,Y vanish on the orthogonal source/reference complement; Y=i(|0><1|-|1><0|) reads positive sin(mt)',
            'phase_readout':'atan2(<Y>,<X>)/t on the positive small-angle branch',
            'complete_channel_error_to_energy_error_exact':str(delta),
            'energy_error_proof':'each Pauli expectation error<=epsilon; phase error<=4epsilon when epsilon<=1/4',
            'certified_denominator_floor_exact':'1/20000','records':rows,'ratio_bounds':ratios,
            'electron_anchor_GeV_exact':str(anchor),
            'anchor_uncertainty_included':False,
            'finite_shot_noise_included':False,'physical_model_error_included':False,
            'precision_is_computed_record_law_not_a_finite_sample_claim':True}


def electron_calibration():
    inputs=json.loads((ROOT/'model/lepton_cycle/experimental_inputs.json').read_text())
    with mp.workdps(60):
        theta=mp.mpf(2)/9;a=mp.sqrt(mp.mpf(5)/108);b=mp.sqrt(mp.mpf(5)/216)
        effects=[(a+2*b*mp.cos(theta+2*mp.pi*k/3))**2 for k in range(3)]
        anchor=mp.mpf(inputs['electron_MeV']['value'])/1000
        scale=anchor/effects[1]
        values=[{'species':name,'character':k,'compiled_source_rest_energy':mp.nstr(effects[k],35),
                 'electron_anchored_GeV':mp.nstr(scale*effects[k],35)} for name,k in (('electron',1),('muon',2),('tau',0))]
        return {'empirical_input':'electron_MeV only; no muon/tau values enter the compiled source',
            'source_energy_unit_GeV':mp.nstr(scale,35),'rows':values,
            'declared_background_norm_GeV_if_unit_is_converted':mp.nstr(scale,35),
            'empirical_weak_scale_or_Higgs_normalization_matched':False,
            'same_compiled_free_source_rest_coefficients':True,
            'observed_interacting_lepton_identification_proved':False,
            'architecture_exposure':'fixed theta2/9 and mode labels inherited from retrospective lepton-cycle architecture; no blind selection claim'}
