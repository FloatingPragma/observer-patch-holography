"""Complete finite phase pointers, confidence sets and exact mass-frame gaps."""
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import math
import mpmath as mp
import numpy as np
from model.particle_joint_mass_identity import joint
from model.native_flavor_orbit import intervals as iv
from model.particle_common_contract.provenance import import_closure

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
M=32;K=6
TAIL=F(K+1,2*K*K)


def bindings():
    paths={str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py','.md')}
    paths.update(('model/particle_joint_mass_identity/results.json','model/particle_joint_mass_identity/PROOF.md',
        'requirements.txt','model/requirements.txt'))
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(import_closure(paths))}


def authenticate(parent):
    for path,digest in parent['source_bindings'].items():
        if hashlib.sha256((ROOT/path).read_bytes()).hexdigest()!=digest:raise ValueError('changed source: '+path)


def number(value):
    value=F(value);return mp.mpf(value.numerator)/value.denominator


def coarse_key(row):
    return (row['source_sector'],row['charge'],row['color_C2'],row['helicity_twice'],row['momentum_z_GeV'])


def range_for(groups):
    high=max(F(r['energy_GeV_interval'][1])-1 for r in groups);width=F(1)
    while width<4*high:width*=2
    return width


def phase_for(square,width):
    # Stable exact mass-frame formula; no binary64 subtraction sqrt(1+m²)-1.
    s=number(square);delta=s/(mp.sqrt(1+s)+1)
    return mp.mpf(1)/4+delta/number(width)


def kernel(phase,size=M):
    nearest=mp.nint(phase*size)
    if phase*size==nearest:
        return [mp.mpc(int(y==int(nearest)%size)) for y in range(size)]
    return [mp.fsum(mp.exp(2j*mp.pi*j*(phase-mp.mpf(y)/size)) for j in range(size))/size for y in range(size)]


def adaptive_kernel(phase,bits=5):
    """One reusable measured qubit; retain every low-bit-first history."""
    if type(bits) is not int or not 1<=bits<=8:raise ValueError('bounded executed pointer requires one through eight bits')
    result=[]
    for y in range(1<<bits):
        amplitude=mp.mpc(1);retained=0
        for j in range(bits):
            bit=(y>>j)&1
            phase_angle=phase*(1<<(bits-1-j))-mp.mpf(retained)/(1<<(j+1))
            amplitude*= (1+(-1)**bit*mp.exp(2j*mp.pi*phase_angle))/2
            retained+=bit*(1<<j)
        if retained!=y:raise ArithmeticError('an adaptive classical bit was discarded')
        result.append(amplitude)
    return result


def confidence_interval(y,width,size=M):
    if size<4*K+1:raise ValueError('pointer requires a guard interval wider than the confidence half-width')
    low=max(F(1,4),F(y-K,size));high=min(F(1,2),F(y+K,size))
    if low>high:return None
    elo=1+width*(low-F(1,4));ehi=1+width*(high-F(1,4))
    lo=iv.sqrt_rational(elo*elo-1,digits=76).lo
    hi=iv.sqrt_rational(ehi*ehi-1,digits=76).hi
    return {'energy_GeV_rational':[str(elo),str(ehi)],'mass_GeV_rational':[str(lo),str(hi)]}


def planning(group,peers,width):
    low,high=map(F,group['energy_GeV_interval']);mlo,mhi=map(F,group['rest_coefficient_GeV_interval'])
    m=(mlo+mhi)/2;dm=m/100 if m>0 else F(1,10**6)
    plus=iv.sqrt_rational(1+(m+dm)**2,digits=76)
    minus=iv.sqrt_rational(1+max(F(),m-dm)**2,digits=76)
    precision=min(plus.lo-high,low-minus.hi) if m>0 else plus.lo-high
    if precision<=0:raise ArithmeticError('mass precision gap not enclosed positively')
    separations=[]
    for other in peers:
        if other['candidate_species']==group['candidate_species']:continue
        a,b=map(F,other['energy_GeV_interval']);gap=max(a-high,low-b)
        if gap<=0:raise ArithmeticError('same-quantum-number source distinction lacks exact energy separation')
        separations.append(gap)
    gap=min(separations) if separations else None
    allowance=min(precision,gap) if gap is not None else precision
    bits=5
    while F(4*K)*width/(1<<bits)>=allowance:bits+=1
    size=1<<bits;exposure=2*mp.pi*(size-1)/number(width)
    return {'mass_goal_half_width_GeV':str(dm),'goal_rule':'one percent of nonzero conditional source mass; 1e-6 GeV at zero',
        'certified_energy_goal_gap_lower_GeV':str(precision),
        'nearest_same_quantum_number_energy_gap_lower_GeV':None if gap is None else str(gap),
        'ideal_pointer_qubits_sufficient':bits,'ideal_pointer_dimension_sufficient':str(size),
        'equivalent_adaptive_quantum_qubits':1,'adaptive_retained_classical_bits':bits,'adaptive_sequential_controlled_evolutions':bits,
        'ideal_maximum_coherent_exposure_GeV_inverse':mp.nstr(exposure,30),
        'ideal_energy_confidence_half_width_GeV':str(F(K)*width/size),
        'confidence_lower_exact':str(1-TAIL),'full_planned_pointer_executed':False,
        'mass_goals_hold_on_good_records':True,'hardware_clock_and_precision_supplied':True}


def calculate():
    pins=bindings();parent=json.loads((ROOT/'model/particle_joint_mass_identity/results.json').read_text());authenticate(parent)
    # Independently replay the same-source charge/color/helicity/momentum commutators.
    f,b=joint.fermion_joint(),joint.boson_joint()
    source=parent['complete_calibrated_joint_instrument'];groups=source['joint_groups'];classes=defaultdict(list)
    for row in groups:classes[coarse_key(row)].append(row)
    result_classes=[];kernels={};budget={};errors={};all_records=[];coverage=[]
    weights={(r['pointer'],r['joint_index']):F(r['probability_exact']) for r in source['complete_records']}
    pointers=list(dict.fromkeys(r['pointer'] for r in source['complete_records']))
    with mp.workdps(100):
        for ci,(key,rows) in enumerate(classes.items()):
            width=range_for(rows);amplitudes=[];intervals=[confidence_interval(y,width) for y in range(M)]
            for row in rows:
                phase=phase_for(F(row['mass_squared_GeV2_exact']),width);a=kernel(phase);kernels[row['joint_index']]=a
                adaptive=adaptive_kernel(phase);adaptive_error=max(abs(x-y) for x,y in zip(a,adaptive,strict=True))
                if adaptive_error>mp.mpf('1e-85'):raise ArithmeticError('adaptive history changes the complete quantum instrument')
                error=abs(mp.fsum(abs(x)**2 for x in a)-1)
                if error>mp.mpf('1e-85'):raise ArithmeticError('finite Fourier pointer is incomplete')
                # This circulant matrix is the whole pointer unitary, not only its blank column.
                matrix=np.array([[complex(a[(y-x)%M]) for x in range(M)] for y in range(M)])
                unitary=float(np.linalg.norm(matrix.conj().T@matrix-np.eye(M)))
                if unitary>2e-13:raise ArithmeticError('complete pointer unitary differs')
                failure=mp.fsum(abs(a[y])**2 for y in range(M) if abs(mp.mpf(y)/M-phase)>mp.mpf(K)/M)
                if failure>number(TAIL)+mp.mpf('1e-70'):raise ArithmeticError('universal confidence bound failed')
                amplitudes.append({'joint_group':row['joint_index'],'phase_decimal':mp.nstr(phase,90),
                    'amplitudes':[[mp.nstr(x.real,90),mp.nstr(x.imag,90)] for x in a],
                    'probabilities':[float(abs(x)**2) for x in a],
                    'adaptive_one_qubit_amplitudes':[[mp.nstr(x.real,90),mp.nstr(x.imag,90)] for x in adaptive],
                    'adaptive_vs_full_pointer_max_amplitude_error_decimal':mp.nstr(adaptive_error,12),
                    'normalization_error_decimal':mp.nstr(error,12),'all_pointer_columns_unitarity_error':unitary})
                coverage.append({'joint_group':row['joint_index'],'actual_ideal_failure_probability':float(failure)})
                if row['candidate_species'][0] not in budget:
                    plan=planning(row,rows,width);mode=max(range(M),key=lambda y:abs(a[y])**2)
                    compatible=[]
                    interval=intervals[mode]
                    if interval:
                        lo,hi=map(F,interval['energy_GeV_rational'])
                        compatible=[p['candidate_species'][0] for p in rows if F(p['energy_GeV_interval'][1])>=lo and F(p['energy_GeV_interval'][0])<=hi]
                    budget[row['candidate_species'][0]]={'source_coefficient_GeV_interval':row['rest_coefficient_GeV_interval'],
                        'calibration_ancestry':next(x for x in parent['conditional_GeV_targets'] if x['species']==row['candidate_species'][0]),
                        'executed_pointer_qubits':5,'executed_range_GeV':str(width),'most_probable_pointer_bin':mode,
                        'most_probable_bin_probability':float(abs(a[mode])**2),'illustrative_modal_record_interval':interval,
                        'modal_record_compatible_source_species':sorted(set(compatible)),'planning':plan}
            entry={'class_index':ci,'source_quantum_numbers':list(key),'joint_group_indices':[r['joint_index'] for r in rows],
                'retained_mode_indices':sum([r['mode_indices'] for r in rows],[]),'range_GeV':str(width),
                'admitted_energy_window_GeV':['1',str(1+width/4)],'phase_window':['1/4','1/2'],
                'pointer_columns':amplitudes,'complete_pointer_confidence_intervals':intervals,
                'maximum_coherent_exposure_GeV_inverse':mp.nstr(2*mp.pi*(M-1)/number(width),30)}
            result_classes.append(entry)
            for pointer in pointers:
                for y in range(M):
                    probability=mp.fsum(number(weights[pointer,r['joint_index']])*abs(kernels[r['joint_index']][y])**2 for r in rows)
                    all_records.append({'probe_pointer':pointer,'coarse_class':ci,'phase_pointer':y,
                        'probability':float(probability),'inconclusive_mass_interval':intervals[y] is None})
        for pointer in pointers:
            for index,kind in [(-1,'vacuum'),(112,'complete_many_excitation_complement')]:
                all_records.append({'probe_pointer':pointer,'coarse_class':kind,'phase_pointer':0,
                    'probability':float(weights[pointer,index]),'inconclusive_mass_interval':True})
        total=math.fsum(r['probability'] for r in all_records)
        if abs(total-1)>3e-15:raise ArithmeticError('original probe failures or complements were lost')
        # A reference-sensitive counterexample to reading exact family labels before the pointer.
        row=next(c for c in result_classes if len(c['joint_group_indices'])==3 and
                 groups[c['joint_group_indices'][0]]['candidate_species']==['neutrino_mass_eigenstate_1'])
        ia,ib=row['joint_group_indices'][:2];aa,ab=kernels[ia],kernels[ib]
        trace_distance=mp.fsum(abs(x*y) for x,y in zip(aa,ab))/2
        if trace_distance<mp.mpf('.49'):raise ArithmeticError('unresolved-family coherence control disappeared')
        control={'family_pair':[groups[i]['candidate_species'][0] for i in [ia,ib]],
            'maximally_entangled_family_reference_half_trace_distance_after_invalid_exact_label_dephasing':float(trace_distance),
            'all_phase_pointer_outputs_summed_without_postselection':True}
    result={'schema':'oph.finite_particle_spectroscopy.v1','source_bindings':pins,'all_checks_pass':True,
        'source_contract':{'source':'unchanged calibrated quadratic two-momentum action and exact native mass frames',
            'momenta_GeV':[[-0.,0.,-1.],[0.,0.,1.]],'one_excitation_dimension':236,
            'full_Fock_dimension_exact':source['source_Hilbert_dimension_exact'],
            'same_original_five_probes':pointers,'charge_color_helicity_and_momentum_measured':True,
            'exact_energy_or_family_label_measured_before_phase_pointer':False,
            'coherent_phase_evolution_and_detector_service_supplied':True,'source_reference_and_clock':source['source_reference_and_clock']},
        'instrument_formula':'K_(a,y) = sum_(g in a) A_y(phi(E_g)) P_g; K_vac=P_vac, K_many=P_many, then compose with each original probe',
        'phase_formula':'phi(E)=1/4+(E-1)/R_a; R_a is the least power of two >= max(1,4 max_g(E_g-1))',
        'finite_pointer_size':M,'finite_pointer_qubits':5,
        'adaptive_pointer':{'exact_all_input_instrument_equivalence':True,'reusable_quantum_qubits':1,
            'retained_original_classical_bits':5,'sequential_controlled_source_evolutions':5,
            'retained_bit_order':'least significant first; y=sum_j b_j 2^j',
            'conditional_Kraus_formula':'K_(b_j|r)=(I+(-1)^b_j exp(-2pi i r/2^(j+1)) U^(2^(n-1-j)))/2',
            'complete_word_identity':'product_j K_(b_j|r_j)=2^-n sum_(k=0)^(2^n-1) exp(-2pi i k y/2^n) U^k',
            'all_original_histories':[{'phase_pointer':y,'bits_low_to_high':[(y>>j)&1 for j in range(5)]} for y in range(M)],
            'same_coherent_exposure_as_parallel_register':True,'complete_source_reference_output_preserved':True,
            'all_112_times32_adaptive_amplitudes_executed':True,
            'known_qubit_recycling':'after readout b, controlled X^b returns the pointer to blank while keeping b in the original classical record',
            'classical_history_erased':False,'postselection':False,'autonomous_blank_reset_derived':False,
            'ideal_angle_and_source_control_services_supplied':True,'large_n_adaptive_precision_plans_executed':False},
        'coarse_classes':result_classes,
        'complete_records':all_records,'record_count':len(all_records),'positive_records':sum(r['probability']>0 for r in all_records),
        'zero_records':sum(r['probability']==0 for r in all_records),'probability_sum':total,
        'uniform_confidence':{'half_width_bins':K,'failure_upper_exact':str(TAIL),'coverage_lower_exact':str(1-TAIL),
            'failure_bound_proof':'|A_y|²<=1/(4 d²), d=circular lattice distance; two tails sum <=(1/K²+1/K)/2',
            'arbitrary_state_and_untouched_reference':True,'hardware_half_diamond_error_adds_to_failure_bound':True,
            'numerical_or_hardware_error_certified':False,'each_group_actual_ideal_failure':coverage},
        'conditional_GeV_menu':budget,'negative_exact_family_lookup_control':control,
        'same_source_replay_errors':{'fermion':f['errors'],'boson':b['errors']},
        'scope':{'calibrated_action_coefficients_are_independent_predictions':False,'interacting_rest_poles_identified':False,
            'neutrino_species_resolved_by_executed_five_qubit_pointer':False,'planned_large_registers_executed':False,
            'conditional_mass_intervals_are_statistical_confidence_sets_for_ideal_detector':True,
            'source_parameter_uncertainty_propagated':False,'native_pulse_compilation_of_phase_detector':False,
            'all_classical_failures_and_full_many_excitation_complement_retained':True,
            'color_fiber_and_within_class_quantum_coherence_preserved':True,
            'complete_instrument_extends_by_identity_on_arbitrary_reference':True}}
    if bindings()!=pins:raise ArithmeticError('source changed during finite spectroscopy')
    return result


def equal(a,b):
    if isinstance(b,dict):return isinstance(a,dict) and a.keys()==b.keys() and all(equal(a[k],v) for k,v in b.items())
    if isinstance(b,list):return isinstance(a,list) and len(a)==len(b) and all(equal(x,y) for x,y in zip(a,b))
    if isinstance(b,float):return isinstance(a,(float,int)) and math.isfinite(a) and abs(a-b)<=1e-14+1e-9*abs(b)
    return type(a)is type(b) and a==b


def main():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True);g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');a=p.parse_args()
    value=calculate()
    if a.write:(HERE/'results.json').write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    elif not equal(json.loads((HERE/'results.json').read_text()),value):raise ValueError('finite spectroscopy receipt differs')
    print('PASS: complete finite phase instrument, confidence sets, same-source GeV menu and precision budgets.')
