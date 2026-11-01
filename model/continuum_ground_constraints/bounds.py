"""Exact bounds, with no numerical ground-state or continuum oracle."""
from fractions import Fraction as F
from math import isqrt


def rational(value, *, positive=False):
    if type(value) not in (int,F) or value<0 or (positive and not value):
        raise ValueError('nonnegative exact rational required')
    return F(value)


def sqrt_upper(value, bits=80):
    value=rational(value)
    if type(bits) is not int or not 1<=bits<=256:
        raise ValueError('integer precision from1 to256 required')
    scale=1<<bits
    root=isqrt(value.numerator*scale*scale//value.denominator)
    if F(root*root,scale*scale)<value:root+=1
    return F(root,scale)


def occupation_defect(ratio, threshold):
    """Bound D for epsilon=U/c and any exact 0<delta<=1."""
    epsilon=rational(ratio);delta=rational(threshold,positive=True)
    if delta>1:raise ValueError('threshold exceeds the proved sine-volume domain')
    return min(F(2),2*delta**3+epsilon/(4*delta))


def best_dyadic(ratio, levels=40):
    epsilon=rational(ratio,positive=True)
    if type(levels) is not int or not 1<=levels<=128:
        raise ValueError('finite calculator requires1..128 levels')
    choices=[(occupation_defect(epsilon,F(1,2**j)),F(1,2**j))
             for j in range(levels+1)]
    defect,delta=min(choices)
    return {'ratio':epsilon,'threshold':delta,'occupation_defect_upper':defect,
            'density_lower':max(F(1,4),1-defect),'density_upper':F(1),
            'per_spin_density_error_upper':min(F(3,8),defect/2),
            'double_occupancy_error_upper':min(F(1,4),8*sqrt_upper(defect)),
            'onsite_density_matrix_trace_norm_error_upper':min(F(2),2*defect+32*sqrt_upper(defect))}


def wick_error_coefficient(order):
    if type(order) is not int or order<0 or order%2 or order>100:
        raise ValueError('even local word length from0 to100 required')
    coefficient=0
    for n in range(2,order+1,2):coefficient=2+(n-1)*coefficient
    return coefficient


def moment_escape(n):
    """Positive probability measures with second-moment loss in a weak limit."""
    if type(n) is not int or n<2:raise ValueError('integer n>=2 required')
    weight=F(1,n*n)
    return {'remote_energy':n,'remote_weight':weight,
            'first_moment':F(1,n),'second_moment':F(1),
            'memory_Gram':1-weight,
            'resolvent_error_at_i_upper':2*weight,
            'weak_limit_second_moment':0,'weak_limit_memory_Gram':0}


def momentum_concentration(n):
    """A shrinking momentum cube may retain a normalized packet completely."""
    if type(n) is not int or n<2:raise ValueError('integer n>=2 required')
    volume=F(1,n**3)
    return {'cube_volume':volume,'integrated_occupation_defect':volume,
            'normalized_packet_density':n**3,
            'packet_defect_expectation':volume*n**3}


def certificate():
    return {'energy_excess_bound':'0<=e_kin-e_free+U*d=e_U-e_free<=U/4',
            'double_occupancy_interval':['0','1/4'],
            'sine_pocket_volume_upper':'delta^3 for0<delta<=1, normalized Brillouin measure',
            'occupation_defect_bound':'D<=min(2,2 delta^3+U/(4c delta))',
            'optimized_bound':'D<=8(U/(24c))^(3/4) when0<U<=24c',
            'density_limit':'n->1, r->1/2 as U/c->0',
            'fixed_local_word_Wick_error':'|omega(A1...An)-omega_free(A1...An)|<=C_n sqrt(D)',
            'Wick_coefficients':{str(n):wick_error_coefficient(n) for n in range(0,12,2)},
            'double_occupancy_limit':'d->1/4 as U/c->0',
            'first_smeared_energy_moment_limit':'<f,h_W f>+(U0/2)||f||^2 when c=v/a,U(a)->U0<infinity',
            'second_moment_equality_passed_to_limit':False,
            'macroscopic_smeared_field_Gaussianity_inferred':False,
            'chemical_potential_shift_added':False,
            'physical_mass_or_species_identified':False,
            'refinement_examples':[
                {key:str(value) for key,value in best_dyadic(8*F(1,2**power)).items()}
                for power in (0,8,16,24,32)],
            'moment_escape_controls':[
                {key:str(value) for key,value in moment_escape(n).items()} for n in (2,8,64)],
            'momentum_concentration_controls':[
                {key:str(value) for key,value in momentum_concentration(n).items()} for n in (2,8,64)]}
