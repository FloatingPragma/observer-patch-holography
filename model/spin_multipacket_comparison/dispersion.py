"""Exact best scalar energy-character defect of the actual local charge inputs."""
from fractions import Fraction as F
import numpy as np

from model.protected_defect_catalog import charge
from model.centered_particle_spectrum.certificate import decimal_bounds


def source_projectors():
    op=charge.operators();identity=np.eye(256,dtype=np.int64)
    group=(0,3,5,6)
    permutations={a:np.array([sum(((word>>j)&1)<<(j^a)for j in range(8))for word in range(256)])for a in group}
    for a,p in permutations.items():
        if not np.array_equal(op['H'][np.ix_(p,p)],op['H']):raise ArithmeticError('cube reflection does not preserve source clock')
        if not np.array_equal(op['twice_charge'][np.ix_(p,p)],op['twice_charge']):raise ArithmeticError('cube reflection mixes the source charge')
        for b in group:
            if not np.array_equal(p[permutations[b]],permutations[a^b]):raise ArithmeticError('cube reflections fail their group law')
    numerators=[sum((-1)**((k&a).bit_count())*identity[permutations[a]]for a in group)for k in range(4)]
    if not np.array_equal(sum(numerators),4*identity):raise ArithmeticError('complete cube-reflection characters lost a sector')
    for i,p in enumerate(numerators):
        for j,q in enumerate(numerators):
            if not np.array_equal(p@q,4*p if i==j else np.zeros_like(p)):
                raise ArithmeticError('source character projectors are not orthogonal')
    return op,numerators


def divide_interval(numerator,denominator):
    lo,hi=denominator
    if lo<=0:raise ValueError('a certified positive denominator is required')
    corners=[a/b for a in numerator for b in denominator]
    return min(corners),max(corners)


def variance_interval(weight,first,second):
    if weight[0]<=0:raise ValueError('a nonzero original character outcome is required')
    squares=[x*x for x in first]
    low_square=F()if first[0]<=0<=first[1]else min(squares)
    low=max(F(),second[0]-max(squares)/weight[0])
    high=max(F(),second[1]-low_square/weight[1])
    return low,high


def build():
    op,projectors=source_projectors();g,error,certificate=charge.ground_input()
    norm=int(g@g);scale=op['scale'];h=op['H']+6*scale*np.eye(256,dtype=np.int64)
    rows=[]
    for q,twice_a in sorted(op['twice_creation'].items()):
        wlo=F();whi=F();rlo=F();rhi=F();characters=[]
        for k,p in enumerate(projectors):
            # P_k A_q g = (p/4)(twice_a/2)g. All256 rows retained.
            image=p.astype(object)@(twice_a.astype(object)@g)
            h_image=h.astype(object)@image
            weight=F(int(image@image),64*norm)
            first=F(int(image@h_image),64*norm*scale)
            second=F(int(h_image@h_image),64*norm*scale*scale)
            wi=(weight-2*error,weight+2*error)
            fi=(first-14*error,first+14*error)
            si=(second-98*error,second+98*error)
            if wi[0]<=0:raise ArithmeticError('a retained source character has no certified weight')
            residual=variance_interval(wi,fi,si)
            if residual[0]<=0:raise ArithmeticError('source conditional energy variance was not certified')
            conditional_mean=divide_interval(fi,wi)
            e0=tuple(F(x)for x in certificate['ground_energy_interval'])
            excitation=(conditional_mean[0]-6-e0[1],conditional_mean[1]-6-e0[0])
            conditional_variance=divide_interval(residual,wi)
            characters.append({'even_cube_reflection_character':k,
                'original_probability_interval':decimal_bounds(*wi,12),
                'conditional_excitation_mean_interval':decimal_bounds(*excitation,12),
                'conditional_energy_variance_interval':decimal_bounds(*conditional_variance,12),
                'unnormalized_optimal_squared_defect_interval':decimal_bounds(*residual,12)})
            wlo+=wi[0];whi+=wi[1];rlo+=residual[0];rhi+=residual[1]
        floor=F(49,1000)if q==1 else F(77,1000)
        if rlo<=floor:raise ArithmeticError('simple original-weight source defect floor failed')
        rows.append({'internal_charge':q,'characters':characters,
            'original_total_probability_interval':decimal_bounds(wlo,whi,12),
            'optimal_unnormalized_squared_dispersion_defect_interval':decimal_bounds(rlo,rhi,12),
            'optimal_conditional_squared_dispersion_defect_interval':decimal_bounds(*divide_interval((rlo,rhi),(wlo,whi)),12),
            'strict_simple_original_weight_defect_lower':str(floor)})
    return {'source_ground_certificate':certificate,
            'source_even_cube_reflection_group':[[0,0,0],[0,1,1],[1,0,1],[1,1,0]],
            'actual_H_charge_reflection_commutation_exact':True,
            'character_projector_numerators_and_completeness_checked':True,
            'four_characters_after_each_of_two_coherent_charge_preparations':8,
            'open_cube_reflection_labels_claimed_as_translation_momenta':False,
            'all_original_records_retained':True,'all_source_basis_rows_retained':256,
            'centered_H_operator_norm_upper':'7','centered_H_definition':'H_B+6I',
            'channels':rows,
            'any_scalar_energy_function_on_these_four_characters_excluded_for_this_full_preparation':True,
            'source_spectral_subbranch_or_distinct_interpolator_excluded':False,
            'finite_selection_field_replaced_by_zero_field_continuum':False,
            'full_native_joint_energy_character_tape_executed':False,
            'physical_dispersion_or_mass_fit_supplied':False}
