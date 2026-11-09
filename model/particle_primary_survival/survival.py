"""Actual Wilson-loop contrast and suppression of invariant Fock extremes.

All large precision and probability exponents remain factored. This module
does not allocate the spatial source or normalize a rare sector.
"""
from fractions import Fraction as F
from functools import lru_cache

import sympy as s

from model.particle_common_contract import contract, provenance
from model.particle_physical_boundary_state import boundary
from model.particle_primary_field_algebra import algebra
from model.particle_filtered_probe import filter as apparatus
from model.particle_positive_window import window
from model.particle_field_preparation import preparation
from model.particle_local_wilson_density import density
from model.native_quantum_gauss import algebra as native
from model.particle_spatial_ground_state import ground

SOURCE_INPUTS = (
    'model/particle_primary_survival/__init__.py',
    'model/particle_primary_survival/__main__.py',
    'model/particle_primary_survival/survival.py',
    'model/particle_primary_survival/test_survival.py',
    'model/particle_primary_survival/PROOF.md',
    'model/particle_physical_boundary_state/PROOF.md',
    'model/particle_spatial_ground_state/PROOF.md',
    'model/particle_filtered_probe/PROOF.md',
    'model/particle_primary_field_algebra/PROOF.md',
    'model/particle_positive_window/PROOF.md',
    'model/particle_yang_mills_gap_contract/PROOF.md',
)


def dyadic(coefficient, negative_exponent):
    value = ground._q(coefficient,'nonnegative dyadic coefficient')
    exponent = ground._i(negative_exponent,'negative binary exponent',0)
    return {'coefficient':str(value),'base':2,'negative_exponent':exponent}


def dyadic_leq(left,right):
    """Exact comparison without constructing either exponential denominator."""
    a=ground._q(left['coefficient'],'left coefficient')
    b=ground._q(right['coefficient'],'right coefficient')
    p=ground._i(left['negative_exponent'],'left exponent',0)
    q=ground._i(right['negative_exponent'],'right exponent',0)
    if left['base'] != 2 or right['base'] != 2:
        raise ValueError('binary factored values required')
    if a == 0:return True
    if b == 0:return False
    return apparatus._ceil_log2(a/b) <= p-q


@lru_cache(maxsize=1)
def native_witness_control():
    """Independent actual generator-metric contraction, without a sector H."""
    currents=native.native_currents();count=len(currents)
    metric=s.Matrix(count,count,lambda i,j:-s.trace(currents[i]*currents[j]))
    inverse=metric.inv()
    casimir=-sum((inverse[i,j]*currents[i]*currents[j]
                  for i in range(count) for j in range(count)),s.zeros(15))
    casimir=apparatus.clean(casimir)
    inherited={int(row['q6Y']):row for row in preparation.native_casimirs()}
    rows=[]
    for center,indices in native.native_blocks():
        charge=int(-2*s.Rational(center));block=casimir.extract(indices,indices)
        value=s.simplify(s.trace(block)/len(indices))
        if block != value*s.eye(len(indices)) or value != inherited[charge]['trace_metric_C2']:
            raise ArithmeticError('actual native metric and block Casimir disagree')
        rows.append({'integer_Q6_charge':charge,'dimension':len(indices),'trace_metric_Casimir':str(value)})
    beta=density.gamma_data()[3];minus=(s.eye(4)-beta)/2
    trial_occupation=s.kronecker_product(minus,s.eye(15))
    onsite=s.kronecker_product(3*beta,s.eye(15))
    if minus*minus != minus or s.trace(minus) != 2:
        raise ArithmeticError('the native negative-beta subspace changed')
    if any(s.trace(t) != 0 for t in currents):
        raise ArithmeticError('the native fully filled determinant is charged')
    if any(apparatus.clean(trial_occupation*s.kronecker_product(s.eye(4),t)
                           -s.kronecker_product(s.eye(4),t)*trial_occupation) != s.zeros(60)
           for t in currents):
        raise ArithmeticError('the filled negative-beta subspace is not gauge invariant')
    return {'native_generator_count':count,'native_internal_dimension':15,
            'Dirac_CAR_modes_per_cell':60,'all_native_blocks':sorted(rows,key=lambda row:row['integer_Q6_charge']),
            'trace_metric_Gram_invertible':True,
            'sum_squared_block_Casimirs':str(sum(F(row['trace_metric_Casimir'])**2 for row in rows)),
            'trial_occupied_modes_per_cell':int(s.trace(trial_occupation)),
            'trial_onsite_energy_times_spacing':str(s.trace(onsite*trial_occupation)),
            'empty_and_full_onsite_energy':str(s.trace(onsite)),
            'actual_spatial_H_or_state_executed':False,
            'loop_character_span_claimed_invariant_under_H':False}


def extreme_sector_bound(member):
    """Combined original probability of the two exact reducing Fock sectors."""
    contract.verify_member(member)
    source=member['source'];a=F(source['spacing']);volume=source['side']**3
    plaquettes=source['interior_plaquettes'];kappa=F(source['magnetic_coefficient'])
    trial=F(native_witness_control()['trial_onsite_energy_times_spacing'])*volume/a
    sector_floor=-kappa*plaquettes/a
    excess=F(member['state']['target_ground_excess']);separation=sector_floor-trial
    if separation <= 0:raise ArithmeticError('the extreme-sector separation is not positive')
    upper=min(F(1),excess/separation)
    n=source['level']
    coarse=F(7,627*2**(8*n+7))
    if upper > coarse:raise ArithmeticError('the refinement suppression bound changed')
    return contract.serial({'member_id':member['member_id'],**member['ids'],
        'sector':'N_F=0 or N_F=60*V, with the full physical gauge carrier',
        'both_sectors_reduce_actual_H_and_Gauss_projector':True,
        'sector_group':'full native S(U(3) x U(2)); not a compact-simple factor selection',
        'global_ground_energy_upper':trial,'extreme_sector_H_lower':sector_floor,
        'extreme_sector_ground_relative_energy_lower':separation,
        'combined_original_sector_probability_upper':upper,
        'coarse_original_probability_upper':coarse,
        'coarse_formula':'7/(627*2^(8*n+7))',
        'conditional_sector_spectrum_transferred_to_primary':False,
        'full_spatial_sector_spectrum_computed':False,
        'GeV_scale':None})


def factored_apparatus(contrast_coefficient, floor_exponent):
    """Use the inherited complete Ramsey errors at T=t=1/W and R_H<=W."""
    c=ground._q(contrast_coefficient,'strict contrast coefficient',True)
    n=ground._i(floor_exponent,'state floor exponent',0)
    cell_bits=max(1,n+apparatus._ceil_log2(16/c))
    step_bits=max(0,-(-(n+apparatus._ceil_log2(F(10,3)/c))//2))
    target=dyadic(c/4,n)
    quadrature=dyadic(4,cell_bits)
    clocks=dyadic(F(5,6),2*step_bits)
    if not dyadic_leq(quadrature,target) or not dyadic_leq(clocks,target):
        raise ArithmeticError('the exact finite contrast budget does not close')
    return {'quadrature_cells':{'base':2,'exponent':cell_bits},
            'Cayley_subdivisions':{'base':2,'exponent':step_bits},
            'contrast_quadrature_error_upper':quadrature,
            'contrast_clock_error_upper':clocks,
            'total_contrast_error_upper':dyadic(c/2,n),
            'compiler_halfdiamond_error':0,'preparation_halftrace_error':0,
            'apparatus_route':'declared exact algebraic neutral dilation, SELECT, PREP and original-H Cayley services',
            'source_energy_zero_changed':False,'reference_actively_counter_evolved':False,
            'all_selector_probe_Ramsey_outcomes_retained':True,
            'counts_or_spatial_matrices_materialized':False}


def selected_witness(member):
    """A proved nonzero finite primary response, at level-dependent times."""
    contract.verify_member(member)
    source=member['source'];a=F(source['spacing'])
    lo=F(member['preparation_plan']['source']['spectral_lower'])
    hi=F(member['preparation_plan']['source']['spectral_upper']);width=hi-lo
    if not lo <= 0 <= hi or width <= 0:raise ArithmeticError('source interval must contain zero')
    duration=1/width
    probe=algebra.probe(member,half_width=str(duration),kind='plaquette')
    algebra.verify_probe(member,probe)
    if not probe['active']:raise ValueError('the centered plaquette is inactive at this level; use n>=1')
    floor=boundary.common_payload(member)['physical_density_floor_factored']['negative_exponent']
    squared=F(native_witness_control()['sum_squared_block_Casimirs'])
    electric=F(source['electric_coefficient'])
    # Each of the four loop links has C_r; ReTr/15 has coefficient1/30.
    k=2*(2*electric/30)**2*squared
    multiplier=F(23,24)**4
    contrast=F(11,24)*multiplier*k/(a*a*width*width)
    band=(multiplier*k-64)/(2*a*a*width*width)
    finite_band=band-contrast/4
    if finite_band <= 0 or width < 8/a:raise ArithmeticError('finite band constants are inconsistent')
    finite=factored_apparatus(contrast,floor)
    return contract.serial({'member_id':member['member_id'],**member['ids'],
        'probe':probe,'source_spectral_width_upper':width,
        'physical_blank':'empty Fock and constant links; only a matrix-element witness, not the prepared state',
        'blank_commutator_norm_squared':k/(a*a),
        'magnetic_compressed_commutator_annihilates_blank':True,
        'general_compressed_multiplication_commutation_claimed':False,
        'selected_polynomial_is_passive':True,
        'selected_density_floor':dyadic(1,floor),
        'original_second_moment_lower':dyadic(k/(a*a),floor),
        'triangle_multiplier_lower_on_full_finite_spectrum':multiplier,
        'filter_half_width':duration,'Ramsey_source_time':duration,
        'continuous_original_contrast_lower':dyadic(contrast,floor),
        'finite_apparatus_original_contrast_lower':dyadic(contrast/2,floor),
        'original_positive_band':{'lower':dyadic(8/a,(floor+1)//2),'upper':width},
        'original_positive_band_weight_lower':dyadic(band,floor),
        'finite_filter_transfer_band_error_upper':dyadic(contrast/4,floor),
        'finite_apparatus_original_positive_band_weight_lower':dyadic(finite_band,floor),
        'band_projector_domain':'original source energy TRANSFERS E_j-E_i on B sqrt(rho)',
        'apparatus':finite,
        'actual_primary_contrast_lower_bound_derived':True,
        'full_primary_response_or_spatial_matrix_executed':False,
        'time_filter_and_energy_band_depend_on_level':True,
        'fixed_positive_energy_band_survival_proved':False,
        'postselection_normalization_applied':False,
        'physical_mass_or_GeV_assignment':None})


def bind_member(member):
    result=selected_witness(member)
    return {'schema':'oph.primary_finite_survival.v1','result':result,'result_sha256':contract.digest(result)}


def verify_bound(member,value):
    if value != bind_member(member):raise ValueError('the finite survival experiment or its source changed')


def fixed_band_requirement(member,contrast_estimate,contrast_error,**kwargs):
    """Reuse the existing fixed-band admission criterion; no witness is inserted."""
    return window.selected_requirement(member,contrast_estimate,contrast_error,**kwargs)


def passive_polynomial_probabilities(energies,beta='1',upper='4',order=3,power=2):
    """Small exact finite validation of the selected functional form only."""
    values=[ground._signed(e,'energy') for e in energies]
    b=ground._q(beta,'positive beta',True);u=ground._signed(upper,'upper energy')
    m=ground._i(order,'positive Taylor order');scale=ground._i(power,'positive scale')
    if not values or max(values)>u:raise ValueError('a valid source upper energy is required')
    q=[]
    for energy in values:
        x=b*(u-energy)/(2*scale);term=F(1);total=term
        for j in range(1,m+1):term*=x/j;total+=term
        q.append(total**(2*scale))
    total=sum(q)
    return tuple(weight/total for weight in q)


def certificate():
    paths=provenance.definition_inputs(SOURCE_INPUTS);hashes=provenance.definition_hashes(paths)
    result={'schema':'oph.primary_survival_certificate.v1',
            'native_source_witness':native_witness_control(),
            'finite_primary_witnesses':[bind_member(contract.member(n)) for n in (1,2)],
            'extreme_sector_bounds':[extreme_sector_bound(contract.member(n)) for n in (0,1,2,4)],
            'source_sha256':hashes,'uniform_noncollapse_or_mass_prediction':False}
    if provenance.definition_hashes(paths) != hashes:raise ValueError('survival inputs changed during computation')
    return result
