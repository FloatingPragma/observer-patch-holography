"""Gap-free cooling bounds for the full spatial gauge-coordinate family.

Exact rational plans do not allocate their exponentially large Hilbert spaces.
The208-state native control evaluates its polynomial on every reducing block.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
import math

import sympy as s

from model.particle_field_preparation import preparation
from model.particle_field_state_response import state_response as thermal

SOURCE_INPUTS = (
    'model/particle_observer_outputs/field_route.py',
    'model/particle_observer_outputs/FIELD_ROUTE.md',
    'model/particle_field_preparation/preparation.py',
    'model/particle_field_preparation/PROOF.md',
    'model/particle_field_state_response/state_response.py',
    'model/particle_field_state_response/PROOF.md',
    'model/particle_ground_spectral_records/records.py',
    'model/particle_smooth_sign_field/PROOF.md',
    'model/thermodynamic_ground/README.md',
    'Lean/Thermodynamics/GibbsEntropyLimit.lean',
    'Lean/Thermodynamics/GapFreeCooling.lean',
    'Lean/Thermodynamics/LowTemperatureControl.lean',
)


def _q(value, name, positive=False):
    return preparation.rational(value, name, positive)


def _signed(value, name):
    if type(value) not in (int, str, F):
        raise ValueError(name+' requires exact rational data')
    return F(value)


def _i(value, name, minimum=0):
    return preparation.integer(value, name, minimum)


def spatial_source_bounds(side=2, link_degree=1, spacing='1', electric='1',
                          magnetic='1/7', regulator='wilson', smoothing='1'):
    """Product-per-link word cutoff, full60CAR/cell, no dense allocation.

    Every cell owns three outgoing links. Open boxes omit interactions with
    support outside the box; periodic boxes have the same conservative bounds.
    """
    side=_i(side,'box side',1); degree=_i(link_degree,'per-link word degree')
    a=_q(spacing,'spacing',True); e=_q(electric,'electric coefficient')
    k=_q(magnetic,'magnetic coefficient'); delta=_q(smoothing,'smoothing',True)
    preparation.native_casimirs()  # Validate the actual trace-metric source.
    if regulator not in ('wilson','smooth'):
        raise ValueError('regulator must be wilson or smooth')
    volume=side**3; links=plaquettes=3*volume
    kinetic=F(360*volume)/a if regulator=='wilson' else F(60*volume)*(2+delta)/a
    bounded=kinetic+k*plaquettes/a
    electric_cap=e*F(29,10)*links*degree**2/a
    bits=60*volume+9*degree*links
    return {'side':side,'cells':volume,'links_upper':links,
            'plaquettes_upper':plaquettes,'CAR_modes':60*volume,
            'per_link_word_degree':degree,'cutoff':'product_per_link_native_words',
            'spacing':a,'electric_coefficient':e,'magnetic_coefficient':k,
            'regulator':regulator,'smoothing':delta if regulator=='smooth' else None,
            'dimension_bits_upper':bits,'log_dimension_upper':F(bits),
            'joint_source_reference_bits_upper':2*bits,
            'dimension_factored_upper':{'Fock_base':2,'Fock_power':60*volume,
                                        'link_base':451,'link_power':degree*links},
            'kinetic_norm_upper':kinetic,'electric_norm_upper':electric_cap,
            'magnetic_norm_upper':k*plaquettes/a,
            'spectral_lower':-bounded,'spectral_upper':electric_cap+bounded,
            'spectral_width_upper':electric_cap+2*bounded,
            'all_occupations_retained':True,'all_link_word_states_retained':True,
            'source_Gauss_singlet_projection':False,
            'same_H_clock':True,'large_matrix_allocated':False}


def polynomial_bound(beta, spectral_lower, spectral_upper, order=8, scale=None):
    """Relative amplitude bound for Q=T_m(beta(U-H)/(2s))**s."""
    b=_q(beta,'inverse temperature',True)
    lo=_signed(spectral_lower,'spectral lower'); hi=_signed(spectral_upper,'spectral upper')
    if hi<lo:raise ValueError('spectral interval is reversed')
    m=_i(order,'Taylor order')
    width=hi-lo
    if scale is None:
        power=1
        while F(power)<b*width:power*=2
    else:
        power=_i(scale,'dyadic scale',1)
        if power & (power-1):raise ValueError('scale must be a power of two')
    y=b*width/(2*power)
    if y>F(1,2):raise ValueError('scale must bound beta times the spectral width')
    tail=y**(m+1)/math.factorial(m+1)/(1-y/F(m+2))
    relative=power*tail
    eta=min(F(1),2*relative)
    return {'beta':b,'spectral_lower':lo,'spectral_upper':hi,
            'spectral_width_upper':width,'Taylor_order':m,'dyadic_scale':power,
            'scaled_argument_upper':y,'single_factor_tail_upper':tail,
            'relative_amplitude_error_upper':relative,
            'canonical_state_half_trace_error_upper':eta,
            'polynomial_degree':m*power,
            'matrix_multiplications_upper':m+(power.bit_length()-1),
            'scalar_shift_changes_clock':False}


def cooling_plan(dimension_bits, spectral_lower, spectral_upper, target_excess='1'):
    """Choose finite rational beta and a factored polynomial without a gap."""
    bits=_i(dimension_bits,'dimension bits',1)
    eps=_q(target_excess,'target energy excess',True)
    lo=_signed(spectral_lower,'spectral lower'); hi=_signed(spectral_upper,'spectral upper')
    if hi<lo:raise ValueError('spectral interval is reversed')
    width=hi-lo; beta=2*bits/eps
    eta_target=min(F(1,4),eps/(2*width)) if width else F(1,4)
    order=0
    while True:
        bound=polynomial_bound(beta,lo,hi,order)
        if bound['canonical_state_half_trace_error_upper']<=eta_target:break
        order+=1
    excess=F(bits)/beta+width*bound['canonical_state_half_trace_error_upper']
    if excess>eps:raise ArithmeticError('source cooling budget failed')
    return {**bound,'dimension_bits_upper':bits,'log_dimension_upper':F(bits),
            'target_energy_excess':eps,'mean_ground_energy_excess_upper':excess,
            'target_met':True,'spectral_gap_assumed':False,
            'state_exactly_stationary':True,'reference_H':0,
            'canonical_reference_retained':True,'postselection':False,
            'efficient_in_volume_claimed':False}


def spatial_plan(**kwargs):
    target=kwargs.pop('target_excess','1')
    source=spatial_source_bounds(**kwargs)
    return {'source':source,'cooling':cooling_plan(source['dimension_bits_upper'],
            source['spectral_lower'],source['spectral_upper'],target)}


def _sym(q):
    q=F(q); return s.Rational(q.numerator,q.denominator)


def _matrix_power_dyadic(matrix,power):
    q=matrix
    while power>1:q=q*q;power//=2
    return s.ImmutableMatrix(q)


def prepared_control(beta='4',order=8,scale=None):
    bound=polynomial_bound(beta,-4,22,order,scale)
    return _prepared_control(bound['beta'],bound['Taylor_order'],bound['dyadic_scale'])


@lru_cache(maxsize=4)
def _prepared_control(beta,order,scale):
    bound=polynomial_bound(beta,-4,22,order,scale)
    phases,h=thermal.real_computation_chart()
    blocks=[]; z=s.S.Zero; energy=s.S.Zero
    for indices in thermal.connected_blocks():
        hb=h.extract(indices,indices); eye=s.eye(len(indices))
        x=_sym(beta)*(22*eye-hb)/(2*scale)
        term=eye; factor=eye
        for j in range(1,order+1):term=term*x/j;factor+=term
        q=_matrix_power_dyadic(factor,scale)
        q2=q*q; z+=s.trace(q2);energy+=s.trace(q2*hb)
        blocks.append((indices,q))
    if z<=0:raise ArithmeticError('zero normalization for a positive polynomial')
    return {'amplitude_blocks':tuple(blocks),'norm_squared':z,
            'computational_basis_phases':phases,'mean_energy':energy/z,
            'bound':bound,'source_dimension':208,
            'preparation_trace_error_upper':bound['canonical_state_half_trace_error_upper'],
            'all_occupations_and_fluxes_retained':True,'reference_retained':True,
            'stationary_source_density':True,'joint_purification_stationary':False,
            'native_source_interval':(-4,22)}


def control_density(beta='4',order=8,scale=None):
    """All208 source states; called only when an operator consumer needs rho."""
    prep=prepared_control(beta,order,scale); phases=prep['computational_basis_phases']
    entries={}
    for indices,q in prep['amplitude_blocks']:
        for (i,j),value in (q*q).todok().items():
            if value:entries[indices[i],indices[j]]=phases[indices[i]]*value*s.conjugate(phases[indices[j]])/prep['norm_squared']
    return s.ImmutableSparseMatrix(208,208,entries)


def ground_inequality_defect(energy_excess,operator_norm_squared='1'):
    return _q(energy_excess,'energy excess')*_q(operator_norm_squared,'squared norm')


def serial(value):
    if isinstance(value,F):return str(value)
    if isinstance(value,dict):return {key:serial(item) for key,item in value.items()}
    if isinstance(value,(tuple,list)):return [serial(item) for item in value]
    return value


def certificate(include_control=True):
    result={'schema':'particle_spatial_ground_state_v1',
            'source_inputs':list(SOURCE_INPUTS),
            'spatial_families':[spatial_plan(side=n,link_degree=1,electric='105',
                                 magnetic='1/7',target_excess=F(1,n)) for n in (2,4,8)],
            'thermodynamic_scope':{'spacing_fixed':True,'per_link_cutoff_fixed':True,
                'Wilson_local_ground_subsequence':True,
                'limit_is_effectively_selected':False,'unique_ground_claimed':False,
                'a_to_zero_limit_proved':False,'physical_particle_masses':[]},
            'physical_Lorentz_or_Einstein_claimed':False}
    if include_control:
        result['native208']=native_control_summary()
    return serial(result)


def outward_dyadic(value,bits=40):
    bits=_i(bits,'rounding bits',1); q=F(value); denominator=1<<bits
    numerator=q.numerator*denominator//q.denominator
    return F(numerator,denominator),F(numerator+1,denominator)


def native_control_summary(beta='4',order=8):
    from model.particle_ground_spectral_records.records import global_ground_interval
    prep=prepared_control(beta,order); bound=prep['bound']
    mean=outward_dyadic(prep['mean_energy']); ground=global_ground_interval()
    epsilon=mean[1]-ground[0]
    if epsilon<0:raise ArithmeticError('energy is below a certified global lower bound')
    return {'source_dimension':208,'beta':bound['beta'],
            'Taylor_order':bound['Taylor_order'],'dyadic_scale':bound['dyadic_scale'],
            'polynomial_degree':bound['polynomial_degree'],
            'canonical_state_half_trace_error_upper':bound['canonical_state_half_trace_error_upper'],
            'mean_energy_interval':mean,'ground_energy_interval':ground,
            'mean_ground_energy_excess_upper':epsilon,
            'gap_free_analytic_excess_upper':F(8)/bound['beta']+26*bound['canonical_state_half_trace_error_upper'],
            'source_density_stationary':True,'reference_retained':True,
            'all_occupations_and_fluxes_retained':True,'same_H_and_clock':True,
            'source_Gauss_singlet_projection':False,
            'postselection':False,'rest_mass_interpretation':False}
