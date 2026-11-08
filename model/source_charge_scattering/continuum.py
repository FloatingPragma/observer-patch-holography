"""Separated bounded charge-pair probabilities on the same continuum source."""
from fractions import Fraction as F
from itertools import product
from model.spin_charge_dispersion import bounds,order


def pair_bounds(radius=F(1)):
    one=bounds.band(order.CERTIFIED_M0,radius)
    kappa=F(one['each_charge_original_joint_band_probability_lower'])
    pair_floor=kappa*kappa/4
    energy=27*F(radius)/4
    ceiling=54*F(radius)/(kappa*kappa)
    retained=pair_floor-energy/ceiling
    if retained!=kappa*kappa/8 or retained<=0:
        raise ArithmeticError('original pair energy-tail accounting differs')
    return {'same_source_order_lower':str(order.CERTIFIED_M0),
        'same_source_single_charge_original_weight_lower':str(kappa),
        'spatial_box_mean_pair_probability_lower':str(kappa*kappa),
        'arbitrarily_far_continuous_shift_pair_probability_lower':str(kappa*kappa/2),
        'arbitrarily_far_rational_shift_pair_probability_lower':str(pair_floor),
        'each_pair_unnormalized_first_energy_moment_upper':str(energy),
        'same_source_pair_upper_energy_cut':str(ceiling),
        'original_pair_weight_below_cut_lower':str(retained),
        'cut_includes_zero_energy':True,
        'each_created_pair_charge':[-2,2],
        'all_preparation_records':[list(record) for record in product((-1,0,1),repeat=2)],
        'bound_valid_at_every_large_separation':False,
        'effective_sufficient_separation_computed':False,
        'clustering_or_Fock_factorization_assumed':False}


def exact_mean_controls():
    # A finite periodic configuration control has only some positive shifts.
    # It audits the inference from an average, not the physical spin source.
    a=(F(1,4),F(),F(),F())
    values=tuple(sum(a[j]*a[(j+k)%4] for j in range(4))/4 for k in range(4))
    mean_a=sum(a)/4
    if sum(values)/4!=mean_a**2 or not any(v==0 for v in values):
        raise ArithmeticError('mean positivity must not assert all-shift positivity')
    # Translation-invariant components retain extra invariant norm. No
    # ergodic factorization is justified by the mean theorem alone.
    components=(F(1,16),F(1,4))
    mean=sum(components)/2;mean_square=sum(x*x for x in components)/2
    if mean_square<=mean**2:raise ArithmeticError('nonclustering control must be strict')
    return {'periodic_four_shift_correlations':[str(v) for v in values],
        'periodic_mean_equals_squared_first_moment':str(mean_a**2),
        'invariant_component_mean_correlation':str(mean_square),
        'squared_mean_in_invariant_component_control':str(mean**2),
        'controls_claimed_as_actual_source_states':False}


def report():
    selected=order.certificate()
    if F(selected['selected_magnetization_strict_lower'])!=order.CERTIFIED_M0:
        raise ArithmeticError('order producer and bounded creator must have one source constant')
    return {'replayed_selected_source_order_certificate':selected,
        'source_pair_bounds':[pair_bounds(r) for r in(F(1),F(1,8))],
        'exact_mean_inference_controls':exact_mean_controls(),
        'same_ordered_ground_clock_spatial_translations_and_bounded_creator':True,
        'complete_nine_branch_infinite_source_preparation_constructed':True,
        'positive_original_weight_at_arbitrarily_large_separation':True,
        'finite_energy_pair_probability_derived':True,
        'zero_energy_pair_weight_excluded':False,
        'asymptotic_particle_states_or_scattering_matrix_constructed':False,
        'physical_pair_mass_or_cross_section':None}
