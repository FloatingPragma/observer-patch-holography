"""Certified detector probabilities, retaining finite energy-filter error."""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import hashlib
import json
from math import isqrt
import numpy as np

from model.affine_gauss_charged_clock import certificate as clock
from model.affine_gauss_full_spectrum import certificate as ground
from model.affine_gauss_full_spectrum.roundoff import squared_norm_exact

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).parent
PARENT = 'model/affine_gauss_charged_clock/results.json'
PARENT_SHA = '7c8f0214f1412e38bdc410c607d0a8f8765df323c3a33a778da923f3892e8613'


def roots(value, places=30):
    """Rational square-root enclosure proved by integer inequalities."""
    value=F(value)
    if value<0:raise ValueError('negative root')
    scale=10**places
    numerator=value.numerator*scale*scale
    k=isqrt(numerator//value.denominator)
    hi=k if k*k*value.denominator==numerator else k+1
    return F(k,scale),F(hi,scale)


def projector_bounds(probability_interval, fidelity_lower):
    """Sharp pure-reference fidelity envelope, rounded outwards exactly."""
    lo,hi=map(F,probability_interval);f=F(fidelity_lower)
    if not 0<=lo<=hi<=1 or not 0<=f<=1:raise ValueError('invalid probability/fidelity')
    if lo<=1-f:lower=F(0)
    else:
        lower=max(F(0),roots(f*lo)[0]-roots((1-f)*(1-lo))[1])**2
    if hi>=f:upper=F(1)
    else:upper=min(F(1),(roots(f*hi)[1]+roots((1-f)*(1-hi))[1])**2)
    return lower,upper


def strings(interval):return [str(F(x)) for x in interval]


@lru_cache(None)
def source():
    raw=(ROOT/PARENT).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=PARENT_SHA:raise ValueError('frozen clock changed')
    value=json.loads(raw);clock.validate_pins(value['input_sha256'])
    if value!=clock.certificate():raise ValueError('fresh clock/source replay differs')
    return value


def certificate():
    parent=source();_,vector,lower2,_=ground.replay()
    norm=squared_norm_exact(vector)
    sine=ground.RADIUS/(lower2-ground.CENTER)
    law=parent['finite_rational_gate_record_bounds']
    fidelity=F(law['accepted_source_ground_fidelity_lower'])
    acceptance=tuple(map(F,law['unconditional_accepted_record_interval']))
    rows=[]
    for vertex,block in enumerate(vector.reshape(13,5456)):
        for occupied,part in [(False,block[:1]),(True,block[1:])]:
            trial=squared_norm_exact(part)/norm
            actual=(max(F(0),trial-sine),min(F(1),trial+sine))
            after=projector_bounds(actual,fidelity)
            rows.append({'source_vertex':vertex,'transverse_nonvacuum':occupied,
                'exact_normalized_trial_probability':str(trial),
                'actual_ground_probability_interval':strings(actual),
                'conditional_post_clock_probability_interval':strings(after),
                'unconditional_accepted_then_observed_interval':strings((acceptance[0]*after[0],acceptance[1]*after[1]))})
    assert sum(F(r['exact_normalized_trial_probability']) for r in rows)==1
    for key in ('actual_ground_probability_interval','conditional_post_clock_probability_interval'):
        assert sum(F(r[key][0]) for r in rows)<=1<=sum(F(r[key][1]) for r in rows)
    g=ground.certificate()
    nonvac=projector_bounds(g['actual_ground_nonvacuum_transverse_probability_interval'],fidelity)
    positions=[{'source_vertex':r['source_vertex'],
                'conditional_post_clock_probability_interval':strings(projector_bounds(r['actual_ground_position_probability'],fidelity))}
               for r in g['actual_ground_position_records']]
    pins=dict(parent['input_sha256']);pins[PARENT]=PARENT_SHA
    pins['model/affine_gauss_full_spectrum/results.json']=hashlib.sha256((ROOT/'model/affine_gauss_full_spectrum/results.json').read_bytes()).hexdigest()
    for p in HERE.iterdir():
        if p.suffix in ('.py','.md'):pins[p.relative_to(ROOT).as_posix()]=hashlib.sha256(p.read_bytes()).hexdigest()
    return {'schema':'oph.affine_gauss_clock_observables.v1',
        'parent_receipt':PARENT,'parent_sha256':PARENT_SHA,
        'same_actual_H_preparation_reference_and_clock':True,
        'accepted_source_ground_fidelity_lower':str(fidelity),
        'unconditional_clock_acceptance_interval':strings(acceptance),
        'joint_post_clock_outcomes':rows,
        'post_clock_position_marginals':positions,
        'conditional_post_clock_nonvacuum_probability_interval':strings(nonvac),
        'unconditional_accepted_nonvacuum_probability_interval':strings((acceptance[0]*nonvac[0],acceptance[1]*nonvac[1])),
        'positive_joint_record_lower_bounds':sum(F(r['unconditional_accepted_then_observed_interval'][0])>0 for r in rows),
        'all_source_and_pointer_sectors_retained':True,
        'all_energy_labels_retained':parent['phase_apparatus']['pointer_dimension'],
        'complete_fine_law':{
            'neutral_pump_branch':'3/4, N=0, all other pointers fixed at zero',
            'rejected_phase_label':'||K_k D_z|vac>||^2/4, each rejected k, spatial pointer fixed',
            'accepted_label_and_joint_readout':'||Pi_(x,b) K_k D_z|vac>||^2/4, each accepted k and all26(x,b)',
            'sum':'1, by sum_k K_k†K_k=I and sum_(x,b)Pi_(x,b)=I on the physical N=1 sector',
            'other_number_sectors':'retained with zero probability for this preparation'},
        'measurement_order':'finite energy clock, then charge position and transverse-vacuum flag',
        'uniform_in_extra_autonomous_waiting_time':True,
        'joint_projectors_commute_with_each_other':True,
        'joint_projectors_commute_with_H_claimed':False,
        'exact_ground_projection_or_stationary_post_clock_state_claimed':False,
        'extra_readout_pointer_dimension':26,
        'readout_service':{'unitary':'sum_j Pi_j tensor S26^j, identity off the physical N=1 sector',
            'duration_source_units':'1','Hermitian_log_norm_upper':'pi',
            'extra_interaction':'exp(-it Hraw) L exp(+it Hraw), exp(-iL)=Uread',
            'complete_propagator':'exp(-iHraw) Uread',
            'record_time':'start of the joint readout interval',
            'autonomous_source_paused':False,'passive_native_local_compilation_claimed':False,
            'net_source_work_absolute_upper':'2 ||Hraw|| on the admitted one-charge sector'},
        'additional_readout_log_capacity':'ln(26)',
        'capacity_scope':'logical pointer only; synthesis and environment workspace separate',
        'scope':{'finite_boundary_source':True,'Standard_Model_species':None,'mass_GeV':None,
                 'full_finite_clock_matrix_numerically_executed':False,
                 'complete_H_residual_and_exact_trial_projector_norms_replayed':True},
        'input_sha256':dict(sorted(pins.items()))}
