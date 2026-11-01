"""Exact complete-instrument error, pulse and spectral probability bounds."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from math import factorial

from model.certified_reference import certificate as cert, continuation as tmpl, interval as iv
from model.full_fock_field.local_simulation import split_bound
from model.spectral_readout import finite, kernels
from model.protected_particle_readout import execute
from . import preparation as p, circuits as c


def local(word,template):
    offset,error,variables=cert.consume(word,0,template)
    if offset!=len(word) or variables:raise ArithmeticError('unbound native clock event')
    return error+cert.controls_bound(word,{})


def field_certificate(duration,steps):
    half=duration/(2*steps)
    plans={execute.key(c.emit(tmpl.swap())):tmpl.swap()}
    for group in c.GEOMETRY.groups:
        for _,_,z in group:
            template=c.controlled_hopping_template(z,half)
            plans[execute.key(c.emit(template))]=template
    template=c.controlled_contact_template(half)
    plans[execute.key(c.emit(template))]=template
    blocks=c.controlled_field_blocks(duration,steps,16)
    counts=Counter(instrument.tape for instrument,_ in blocks)
    error=F();rows=[]
    for tape,count in counts.items():
        word=[dict(kind=kind,angle=angle,targets=list(targets)) for kind,angle,targets in tape]
        value=local(word,plans[tape]);error+=count*value
        rows.append({'word_sha256':cert.word_digest(word),'pulses':len(word),
                     'blocks':count,'complete_block_error':str(value)})
    splitting=split_bound(c.GEOMETRY.sides,duration,steps,refined=True)
    return {'duration':str(duration),'steps':steps,'literal_source_error':str(error),
            'full_Fock_controlled_Strang_error':str(splitting),
            'pulses':sum(len(instrument.tape) for instrument,_ in blocks),
            'local_word_certificates':rows}


def ground_phase_error(duration):
    word=c.ground_phase_word(duration)
    root=iv.sqrt_fraction(F(3))
    desired=-duration*root/8
    a,b=word
    if a['kind']!='I' or b['kind']!='Z' or a['targets']!=[0] or b['targets']!=[0] or a['angle']!=-b['angle']:
        raise ArithmeticError('ground energy phase word changed')
    error=2*(iv.Interval.point(F(a['angle']))-desired).absolute_upper()
    return error+cert.controls_bound(word,{})


def fourier_templates():
    result=[]
    for high in reversed(range(c.CLOCK_BITS)):
        result.append(c.hadamard_template())
        for low in reversed(range(high)):
            result.append(c.controlled_phase_template(F(1,1<<(high-low))))
    for _ in range(c.CLOCK_BITS//2):
        result.append(cert.cnot(0,1)+cert.cnot(1,0)+cert.cnot(0,1))
    return result


@lru_cache(maxsize=1)
def certificate():
    prep=p.certificate();error=F(prep['complete_mixed_ground_preparation_error'])
    charge=local(c.emit(c.charge_template()),c.charge_template());error+=charge
    initial=c.CLOCK_BITS*local(c.emit(c.hadamard_template()),c.hadamard_template());error+=initial
    fields=[];phase=F()
    for j,steps in enumerate(c.STEPS):
        duration=c.TAU*(1<<j);row=field_certificate(duration,steps);fields.append(row)
        error+=F(row['literal_source_error'])+F(row['full_Fock_controlled_Strang_error'])
        phase+=ground_phase_error(duration)
    error+=phase
    fourier=sum((local(c.emit(template),template) for template in fourier_templates()),F());error+=fourier
    if error>=F(3,100):raise ArithmeticError('complete native clock missed three-percent target')
    clock_pulses=sum(len(instrument.tape) for instrument,_ in c.clock_blocks())
    prep_pulses=prep['original_ground_tree_pulses']+prep['selector']['pulses']+prep['original_inverse_frame_pulses']+max(x['all_eight_sites_pulses'] for x in prep['rotations'])
    return {'clock_outcomes':c.M,'time_step':str(c.TAU),
            'largest_coherent_history_time':str((c.M-1)*c.TAU),
            'controlled_evolution_durations':[str(c.TAU*(1<<j)) for j in range(c.CLOCK_BITS)],
            'preparation_complete_error':prep['complete_mixed_ground_preparation_error'],
            'charge_instrument_complete_error':str(charge),
            'clock_initialization_error':str(initial), 'controlled_fields':fields,
            'source_derived_ground_energy_phase_error':str(phase),
            'positive_fourier_error':str(fourier),
            'complete_state_and_total_variation_error':str(error),
            'pulse_counts':{'worst_case_preparation':prep_pulses,
                            'charge_instrument':len(c.emit(c.charge_template())),
                            'clock':clock_pulses,
                            'worst_case_complete_history':prep_pulses+len(c.emit(c.charge_template()))+clock_pulses},
            'complete_native_pointer_slots':6**(1+c.CLOCK_BITS),
            'maximum_simultaneously_active_qubits':18,
            'selector_read_before_blank_spin_embedding':True,
            'all_first_export_failures_retained':True,
            'full_Fock_field_error_preserves_arbitrary_references':True,
            'floating_execution_is_not_the_control_error_proof':True}


def ideal_probabilities():
    """Joint charge/bin law: addition has mass3/4, removal mass1/4."""
    down=[]
    for index in range(c.M):
        value=iv.Interval.point(0)
        for lo,hi,low,high in finite.load_source():
            value+=finite.enclosure(low,high)*finite.kernel_at_energy(c.M,c.TAU,index,(lo,hi))
        down.append(value)
    root=iv.sqrt_fraction(F(3))/16
    lam=(F(root.lo,iv.SCALE),F(root.hi,iv.SCALE))
    rows=[];error=F(certificate()['complete_state_and_total_variation_error'])
    centers=kernels.energy_centers(c.M,c.TAU)
    for tag in range(2):
        for index in range(c.M):
            fixed=finite.kernel_at_energy(c.M,c.TAU,index,lam)/4
            ideal=fixed+down[index]/2 if tag==0 else fixed
            lo=max(F(),F(ideal.lo,iv.SCALE)-error);hi=min(F(1),F(ideal.hi,iv.SCALE)+error)
            rows.append({'charge':'addition' if tag==0 else 'removal', 'charge_record':tag,
                         'index':index,'response_frequency_sign':1 if tag==0 else -1,
                         'excitation_energy_bin_center_interval':centers[index].json(),
                         'ideal_joint_probability_interval':ideal.json(),
                         'ideal_joint_probability_decimal':finite.outward_decimal(ideal),
                         'native_unconditioned_probability_interval':[str(lo),str(hi)]})
    return rows


def finite_shot_certificate():
    # Eight charge/bin events and the coarse abort event. A union bound
    # controls every empirical probability, including zero/abort records.
    n=100000;epsilon=F(1,100);alpha=F(1,1000);events=2*c.M+1
    exponent=2*n*epsilon*epsilon
    lower=sum((exponent**j/F(factorial(j)) for j in range(40)),F())
    if 2*events/lower>=alpha:raise ArithmeticError('finite-shot union bound failed')
    return {'independent_histories':n,'events_including_coarse_abort':events,
            'simultaneous_empirical_probability_error':str(epsilon),
            'failure_probability_upper':str(alpha),
            'inequality':'2K exp(-2 N epsilon^2), bounded by an exact positive exponential series',
            'shots_or_random_measurement_records_generated':False}
