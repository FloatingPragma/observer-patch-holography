"""Replay the full Hamiltonian and admit a fixed, outward spectral enclosure."""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import hashlib
import json
import math
import numpy as np
import mpmath as mp

from model.affine_gauss_spectrum import model as coarse
from .operator import HERE, ROOT, ActualOperator, algebraic_frame
from .roundoff import U, q, evaluation_budget, sqrt_bounds, squared_norm_exact

PARENT='model/affine_gauss_spectrum/results.json'
PARENT_SHA='632a3b6ae85c7f1fb67ce768a78b684972831eb1b194517b1dc9c2cb1affb9ff'
CENTER=F.from_float(1.8738097549661659)
RADIUS=F(1,10**12)


def outward(value,upper=False,places=24):
    value=F(value);scale=10**places
    numerator=-(-value.numerator*scale//value.denominator) if upper else value.numerator*scale//value.denominator
    return str(F(numerator,scale))


def interval(lo,hi):return [outward(lo),outward(hi,True)]


def norm_admission(vector):
    if vector.shape!=(70928,) or vector.dtype!=np.dtype('float64') or not np.isfinite(vector).all():
        raise ValueError('the complete finite binary64 trial is required')
    lo,hi=sqrt_bounds(squared_norm_exact(vector))
    if not F(1,2)<=lo<=hi<=2:
        raise ValueError('trial norm outside the certified evaluation range')
    return lo,hi


@lru_cache(None)
def replay():
    raw=(ROOT/PARENT).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=PARENT_SHA:raise ValueError('coarse source certificate changed')
    parent=json.loads(raw);coarse.validate_pins(parent['input_sha256'])
    if parent!=coarse.certificate():raise ValueError('full-complement parent replay differs')
    vector=np.load(HERE/'trial_state.npy',allow_pickle=False)
    mp.iv.dps=80
    normlo,normhi=norm_admission(vector)
    frame=algebraic_frame();operator=ActualOperator(frame)
    image=operator.apply(vector);energy=float(CENTER)
    residual=image-energy*vector
    _,residual_norm=sqrt_bounds(squared_norm_exact(residual))
    _,image_norm=sqrt_bounds(squared_norm_exact(image))
    budget=evaluation_budget(frame,operator)
    final_round=U*abs(CENTER)*normhi+U*(image_norm+(1+U)*abs(CENTER)*normhi)
    residual_upper=(residual_norm+F(budget['H_apply_norm_error_per_unit_input'])*normhi+final_round)/normlo
    if residual_upper>RADIUS:raise ValueError('actual full-H residual exceeds the published spectral radius')
    lower2=F(parent['spectral']['low_charged_levels'][1]['actual_H_energy_interval'][0])
    if CENTER+RADIUS>=lower2:raise ValueError('spectral window does not isolate the actual charged ground')
    return parent,vector,lower2,{'computed_actual_residual_upper':float(residual_upper),
        'computed_H_evaluation_budget':float(F(budget['H_apply_norm_error_per_unit_input'])),
        'published_residual_upper':str(RADIUS)}


def certificate():
    parent,vector,lower2,_=replay()
    norm2=squared_norm_exact(vector)
    vacuum=squared_norm_exact(vector[::5456])/norm2
    sine=RADIUS/(lower2-CENTER)
    compression_lo,compression_hi=map(F,parent['spectral']['low_charged_levels'][0]['zero_photon_compression_energy_interval'])
    photon=interval(max(F(0),1-vacuum-sine),min(F(1),1-vacuum+sine))
    fidelity=interval(max(F(0),vacuum-sine),min(F(1),vacuum+sine))
    positions=[]
    for y,block in enumerate(vector.reshape(13,5456)):
        probability=squared_norm_exact(block)/norm2
        positions.append({'source_vertex':y,'actual_ground_position_probability':
                          interval(max(F(0),probability-sine),min(F(1),probability+sine))})
    # The pump mode is the normalized zero-photon component of the supplied
    # binary trial. It has13 entries and no measured-data input.
    pump_fidelity_lo,pump_fidelity_hi=map(F,fidelity)
    result={'schema':'oph.affine_gauss_full_charged_spectrum.v1',
        'parent_receipt':PARENT,'parent_receipt_sha256':PARENT_SHA,
        'same_H_state_reference_clock':True,
        'physical_one_charge_dimension':70928,'photon_carrier_dimension':5456,
        'photon_total_occupation_cutoff':3,'additional_state_truncation':False,
        'real_basis_unitary':'exp(-i pi N_ph/2); leaves photon number, positions and dGamma(Omega) unchanged',
        'Peierls_exponential':'exact degree8 polynomial on all active spectator blocks of sizes1,2,3,4',
        'trial_is_untrusted_until_full_residual_passes':True,
        'arithmetic_admission':{'standard':'IEEE754 binary64 round-to-nearest',
            'trial_norm_range':['1/2','2'],'coefficient_enclosures':'outward interval arithmetic at80 decimal digits',
            'evaluation':'full70928 source, all photon occupations, sparse/Horner rounding and underflow bound',
            'maximum_admitted_actual_H_residual':str(RADIUS),
            'numerical_eigensolver_used_in_verifier':False},
        'energy_center_exact_binary_value':str(CENTER),
        'actual_charged_ground_energy_interval':interval(CENTER-RADIUS,CENTER+RADIUS),
        'zero_photon_compression_ground_energy_interval':interval(compression_lo,compression_hi),
        'energy_lowering_from_full_transverse_dressing_interval':interval(
            compression_lo-CENTER-RADIUS,compression_hi-CENTER+RADIUS),
        'charged_ground_is_unique':True,
        'charged_ground_gap_lower':outward(lower2-CENTER-RADIUS),
        'reference_energy_retained':'1',
        'source_energy_after_fixed_reference_subtraction':interval(CENTER-RADIUS-1,CENTER+RADIUS-1),
        'normalized_trial_sine_angle_to_actual_ground_upper':outward(sine,True),
        'actual_ground_nonvacuum_transverse_probability_interval':photon,
        'actual_ground_position_records':positions,
        'complete_ground_preparation':{
            'initial_state':'same actual neutral vacuum',
            'mode':'normalized13-entry zero-photon component of trial_state.npy',
            'pump':'D_w/2; exact compensated physical CAR creation',
            'pump_success_probability':'1/4','pump_failure_probability':'3/4',
            'conditional_ground_window_probability_interval':fidelity,
            'unconditional_ground_window_probability_interval':interval(pump_fidelity_lo/4,pump_fidelity_hi/4),
            'unconditional_energy_complement_probability_interval':interval((1-pump_fidelity_hi)/4,(1-pump_fidelity_lo)/4),
            'exact_complete_law':['3/4','f/4','(1-f)/4'],
            'energy_window':interval(CENTER-RADIUS,CENTER+RADIUS),
            'finite_matrix_preparation_and_spectral_detector_services_declared':True,
            'efficient_native_circuit_claimed':False},
        'scope':{'finite_boundary_reference_problem':True,'photon_dressing_is_nonzero':True,
            'translation_invariant_mass_shell_identified':False,'Lorentz_spin_identified':False,
            'Standard_Model_identity':None,'mass_GeV':None}}
    pins=dict(parent['input_sha256']);pins[PARENT]=PARENT_SHA
    for path in HERE.iterdir():
        if path.suffix in ('.py','.md','.npy'):
            pins[path.relative_to(ROOT).as_posix()]=hashlib.sha256(path.read_bytes()).hexdigest()
    result['input_sha256']=dict(sorted(pins.items()))
    return result
