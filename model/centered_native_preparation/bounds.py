"""Exact complete-output errors and finite-resolution charged energy laws."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import json

from ..ground_spectral_source import preparation as frames
from ..native_spectral_detector import bounds as old, circuits as prior_circuits
from ..certified_reference import certificate as cert, interval as iv
from ..spectral_readout import finite
from ..centered_particle_spectrum import packet as spectrum
from ..cube_eta_protection import cube
from . import tree,target,circuits


@lru_cache(maxsize=1)
def preparation_certificate():
    parsed=tree.certificate();frame_error=F();frame_pulses=0
    for site in cube.GEOMETRY.sites:
        plan=frames.inverse_frame_plan(site);word=frames.emit(plan)
        offset,error,variables=cert.consume(word,0,frames.template(plan))
        if variables or offset!=len(word):raise ArithmeticError('unaccounted full-Fock inverse-frame event')
        frame_error+=error+cert.controls_bound(word,{})
        frame_pulses+=len(word)
    source_target=target.summary()
    error=F(parsed['complete_tree_error'])+frame_error+F(source_target['normalized_target_to_exact_ground_error'])
    if error>=F(1,10**7):raise ArithmeticError('native centered preparation exceeds precision domain')
    return {'target':source_target,'Gray_tree':parsed,'inverse_frame':frames.exact_frame_certificate(),
            'inverse_frame_pulses':frame_pulses,'inverse_frame_error':str(frame_error),
            'complete_preparation_error':str(error),'prepared_number':8,
            'ground_is_original_spin_two_N4_reference':False,
            'full_Fock_inverse_frame_executed':True,'postselection_required':False}


def ground_phase_error(duration):
    word=circuits.ground_phase_word(duration);a,b=word
    if a['kind']!='I' or b['kind']!='Z' or a['angle']!=-b['angle']:raise ArithmeticError('ground scalar clock phase changed')
    _,_,ground=target.ground();lo,hi=ground['energy_intervals'][0]
    represented=F(a['angle'])
    return 2*max(abs(represented-duration*lo/2),abs(represented-duration*hi/2))+cert.controls_bound(word,{})


@lru_cache(maxsize=1)
def certificate():
    prep=preparation_certificate();error=F(prep['complete_preparation_error'])
    charge=old.local(prior_circuits.emit(prior_circuits.charge_template()),prior_circuits.charge_template());error+=charge
    initialization=2*old.local(prior_circuits.emit(prior_circuits.hadamard_template()),prior_circuits.hadamard_template());error+=initialization
    fields=[];centering=F();phase=F()
    for j,steps in enumerate(circuits.STEPS):
        duration=circuits.TAU*(1<<j)
        row=old.field_certificate(duration,steps);fields.append(row)
        error+=F(row['literal_source_error'])+F(row['full_Fock_controlled_Strang_error'])
        template=circuits.center_template(duration)
        centering+=16*old.local(prior_circuits.emit(template),template)
        phase+=ground_phase_error(duration)
    fourier=sum((old.local(prior_circuits.emit(t),t) for t in old.fourier_templates()),F())
    error+=centering+phase+fourier
    if error>=F(3,100):raise ArithmeticError('complete centered detector misses three-percent bound')
    pulses=sum(len(instrument.tape) for instrument,_ in circuits.clock_blocks())
    return {'preparation':prep,'charge_instrument_error':str(charge),
            'clock_initialization_error':str(initialization),'controlled_fields':fields,
            'centering_phase_error':str(centering),'certified_ground_energy_phase_error':str(phase),
            'positive_Fourier_error':str(fourier),'complete_state_and_total_variation_error':str(error),
            'clock_pulses':pulses,'complete_pulses':prep['Gray_tree']['pulses']+prep['inverse_frame_pulses']+13+pulses,
            'clock_bins':4,'coherent_history_times':['0','1','2','3'],
            'ground_clock':'K=H-N/2; excitation generator K-E_ground',
            'complete_pointer_slots':216,'nonzero_logical_charge_bin_slots':8,
            'source_control_rounding_is_not_floating_execution_error':True,
            'native_first_failure_outputs_retained':True,
            'continuing_bare_field_repair_included':False,
            'physical_species_or_relativistic_mass_identified':False}


@lru_cache(maxsize=1)
def spectral_probabilities():
    exact=spectrum.certificate.verify(spectrum.load_witness())
    values=[]
    for index in range(4):
        value=iv.Interval.point(0)
        for row in exact['rows']:
            e=tuple(F(x) for x in row['excitation_energy_interval'])
            w=finite.enclosure(*(F(x) for x in row['local_branch_weight_interval']))
            value+=w*finite.kernel_at_energy(4,F(1),index,e)
        values.append(value)
    error=F(certificate()['complete_state_and_total_variation_error'])
    return [{'charge_record':tag,'charge':'addition' if tag==0 else 'removal','clock_index':index,
             'signed_response_frequency_sign':1 if tag==0 else -1,
             'ideal_probability_interval':value.json(),
             'ideal_probability_decimal':finite.outward_decimal(value),
             'native_probability_interval':[str(max(F(),F(value.lo,iv.SCALE)-error)),str(min(F(1),F(value.hi,iv.SCALE)+error))]}
            for tag in range(2) for index,value in enumerate(values)]
