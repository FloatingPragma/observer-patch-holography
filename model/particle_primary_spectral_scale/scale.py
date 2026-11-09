"""Exact primary product-vacuum diagnostics, without a gap extrapolation.

The interpolation changes the action except at coupling one. Its vacuum is
an expansion anchor, never the primary polynomial preparation.
"""
from fractions import Fraction as F
from functools import lru_cache

import sympy as s

from model.particle_common_contract import contract, provenance
from model.particle_primary_survival import survival
from model.particle_spatial_ground_state import ground
from model.particle_local_wilson_density import density

SOURCE_INPUTS = (
    'model/particle_primary_spectral_scale/__init__.py',
    'model/particle_primary_spectral_scale/__main__.py',
    'model/particle_primary_spectral_scale/scale.py',
    'model/particle_primary_spectral_scale/test_scale.py',
    'model/particle_primary_spectral_scale/PROOF.md',
    'model/particle_common_contract/PROOF.md',
    'model/particle_primary_survival/PROOF.md',
    'model/particle_primary_capacity/PROOF.md',
    'model/particle_primary_locality/PROOF.md',
    'model/particle_primary_mass_readout/PROOF.md',
    'model/particle_field_preparation/PROOF.md',
    'model/particle_spatial_ground_state/PROOF.md',
    'model/source_strong_coupling/PROOF.md',
    'model/strong_coupling_refinement/README.md',
)


def casimir(color_p, color_q, twice_weak_spin, integer_charge):
    """A conservative superset of the actual descended representation labels.

    All singlet characters must have Q6 divisible by six. Nonabelian labels
    need not satisfy the extra quotient congruences for this lower bound.
    """
    p=ground._i(color_p,'SU3 highest weight p',0)
    q=ground._i(color_q,'SU3 highest weight q',0)
    j=ground._i(twice_weak_spin,'twice SU2 spin',0)
    if type(integer_charge) is not int:
        raise ValueError('integer Q6 charge required')
    z=integer_charge
    if p==q==j==0 and z%6:
        raise ValueError('a descended singlet character has Q6 in 6Z')
    return F(p*p+q*q+p*q+3*p+3*q,6)+F(j*(j+2),8)+F(z*z,120)


@lru_cache(maxsize=1)
def vacuum_channels():
    """Actual Dirac matrices and all five native Casimir blocks are retained."""
    native=survival.native_witness_control()
    beta=density.gamma_data()[3]
    negative=(s.eye(4)-beta)/2
    occupied=s.diag(negative,negative)
    vacancy=s.eye(8)-occupied
    spin=[]
    for axis,gamma in enumerate(density.gamma_data()[:3]):
        hop=(-beta+beta*gamma)/2
        h=s.zeros(4).row_join(hop).col_join(hop.H.row_join(s.zeros(4)))
        onsite=s.diag(3*beta,3*beta)
        transition=vacancy*h*occupied
        if (h!=h.H or s.trace(occupied*h)!=0
                or s.simplify(s.trace(transition.H*transition))!=1
                or density.clean(onsite*transition-transition*onsite)!=6*transition):
            raise ArithmeticError('actual Wilson bond does not create the stated particle-hole channel')
        spin.append({'axis':axis,'one_internal_column_transition_norm_squared':'1',
                     'onsite_excitation_energy':'6'})
    rows=[]
    for row in native['all_native_blocks']:
        c=F(row['trace_metric_Casimir']);d=row['dimension']
        hop_gap=6+F(105,2)*c
        loop_gap=210*c
        rows.append({**row,'bond_vacuum_transition_norm_squared':d,
            'bond_H0_excitation':hop_gap,
            'bond_inverse_energy_weight':F(d)/hop_gap,
            'bond_inverse_square_weight':F(d)/(hop_gap*hop_gap),
            'both_loop_orientations_norm_squared':F(2,7**2*30**2),
            'loop_H0_excitation':loop_gap,
            'plaquette_inverse_energy_weight':F(2,7**2*30**2)/loop_gap,
            'plaquette_inverse_square_weight':F(2,7**2*30**2)/(loop_gap*loop_gap)})
    total=lambda name:sum((row[name] for row in rows),F(0))
    return contract.serial({
        'all_native_blocks':rows,'actual_Dirac_bond_checks':spin,
        'full_CAR_modes_per_cell':60,
        'all_word_nontrivial_link_Casimir_lower':F(3,10),
        'dimensionless_electric_link_gap':F(63,4),
        'dimensionless_decoupled_ambient_gap':3,
        'dimensionless_decoupled_physical_gap':6,
        'bond_vacuum_variance':total('bond_vacuum_transition_norm_squared'),
        'plaquette_vacuum_variance':total('both_loop_orientations_norm_squared'),
        'bond_inverse_energy_weight':total('bond_inverse_energy_weight'),
        'bond_inverse_square_weight':total('bond_inverse_square_weight'),
        'plaquette_inverse_energy_weight':total('plaquette_inverse_energy_weight'),
        'plaquette_inverse_square_weight':total('plaquette_inverse_square_weight'),
        'interpolation_vacuum_is_primary_selected_state':False,
        'finite_spatial_H_or_primary_density_executed':False})


def selected_expansion(member):
    """Bind finite perturbation coefficients to the unchanged primary endpoint."""
    contract.verify_member(member)
    source=member['source'];side=source['side'];volume=side**3
    a=F(source['spacing']);bonds=3*side*side*(side-1)
    plaquettes=source['interior_plaquettes'];kappa=F(source['magnetic_coefficient'])
    local=vacuum_channels()
    def sum_terms(suffix):
        return bonds*F(local['bond_'+suffix])+plaquettes*F(local['plaquette_'+suffix])
    strength=30*bonds+kappa*plaquettes
    inverse=sum_terms('inverse_energy_weight')
    inverse_square=sum_terms('inverse_square_weight')
    return contract.serial({'member_id':member['member_id'],**member['ids'],
        'spacing':a,'cells':volume,'actual_interior_bonds':bonds,
        'actual_interior_plaquettes':plaquettes,
        'interpolation':'Hhat(lambda)=a*(H_onsite+H_electric)+lambda*a*(H_hopping+H_magnetic)',
        'primary_endpoint_lambda':'1',
        'anchor_state':'filled two negative-beta native rows per cell, every link constant',
        'anchor_state_is_selected_primary_preparation':False,
        'dimensionless_anchor_ground_energy':-90*volume,
        'dimensionless_anchor_physical_gap':6,
        'actual_source_energy_of_anchor_expectation':-90*volume/a,
        'dimensionless_primary_H_variance_on_anchor':sum_terms('vacuum_variance'),
        'second_order_energy_coefficient':-inverse,
        'first_order_orthogonal_dressing_norm_squared':inverse_square,
        'second_order_energy_coefficient_per_cell':-inverse/volume,
        'first_order_dressing_norm_squared_per_cell':inverse_square/volume,
        'second_order_remainder_bound_at_primary_endpoint':None,
        'dimensionless_global_perturbation_norm_upper':strength,
        'global_Weyl_gap_lower_at_primary_endpoint':6-2*strength,
        'strict_lambda_threshold_for_positive_global_Weyl_bound':F(3)/strength,
        'incident_interaction_norm_upper':F(1272,7),
        'incident_norm_over_anchor_physical_gap':F(212,7),
        'uniform_local_gap_stability_proved':False,
        'primary_gap_or_fixed_band_survival_proved':False,
        'source_energy_zero_or_clock_changed':False,
        'physical_mass_GeV':None})


def fixed_band_scale(member,lower_energy='1/4',upper_energy='4'):
    """Exact dimensionless window and cooling scale, with fixed source units."""
    contract.verify_member(member)
    lower=ground._q(lower_energy,'positive lower source energy',True)
    upper=ground._q(upper_energy,'upper source energy',True)
    if upper<=lower:raise ValueError('ordered nondegenerate source band required')
    a=F(member['source']['spacing']);eps=F(member['state']['target_ground_excess'])
    return contract.serial({'member_id':member['member_id'],**member['ids'],
        'source_band':[lower,upper], 'dimensionless_band':[a*lower,a*upper],
        'Hhat_definition':'a*H on the identical physical carrier',
        'dimensionless_selected_mean_ground_excess_upper':a*eps,
        'source_time_one_in_dimensionless_time':1/a,
        'declared_GeV_scale':None,
        'actual_dimensionless_low_band_weight_computed':False,
        'critical_running_trajectory_inferred':False})


def conditional_gap_exclusion(member,dimensionless_gap_lower,
                              lower_energy='1/4',upper_energy='4',norm_squared='1'):
    """A SUPPLIED full physical gap suppresses fixed-band original transfers.

    The package does not prove the gap for the primary endpoint.
    """
    scale=fixed_band_scale(member,lower_energy,upper_energy)
    gap=ground._q(dimensionless_gap_lower,'supplied positive dimensionless physical gap',True)
    cap=ground._q(norm_squared,'uniform probe norm squared')
    a=F(member['source']['spacing']);eps=F(member['state']['target_ground_excess'])
    admitted=F(upper_energy)<gap/a
    bound=min(cap,cap*a*eps/gap) if admitted else None
    return contract.serial({**scale,'supplied_dimensionless_gap_lower':gap,
        'gap_is_external_unproved_premise_for_primary':True,
        'ground_degeneracy_allowed':True,'probe_norm_squared_upper':cap,
        'band_strictly_below_supplied_gap':admitted,
        'conditional_original_band_weight_upper':bound,
        'success_normalization_applied':False,
        'actual_primary_gap_or_collapse_proved':False})


def polynomial_countercontrol(level):
    """A distinct three-level Q^2 source with rare fixed-energy transfers."""
    n=ground._i(level,'control level',0);a=F(1,2**n);eps=a/16
    energies=(F(0),64/a,64/a+1)
    probabilities=survival.passive_polynomial_probabilities(
        energies,beta=64/a,upper=energies[-1],order=1,power=1)
    excess=sum(p*e for p,e in zip(probabilities,energies))
    if not excess<=eps:raise ArithmeticError('the polynomial control exceeds its cooling budget')
    # A exchanges the two excited states; all other matrix entries vanish.
    positive=probabilities[1];negative=probabilities[2]
    upper=eps*a/64
    if positive>upper:raise ArithmeticError('the full physical-gap bound fails on its control')
    return contract.serial({'comparison_control_only':True,'is_selected_primary':False,
        'spacing':a,'energies':energies,'Q_squared_probabilities':probabilities,
        'mean_ground_excess':excess,'cooling_tolerance':eps,
        'dimensionless_full_gap':64,'positive_transfer_frequency':1,
        'original_positive_fixed_band_weight':positive,
        'original_negative_transfer_weight':negative,
        'complete_filter_rejection_weight':probabilities[0],
        'conditional_gap_original_weight_upper':upper,
        'postselected_positive_fraction':positive/(positive+negative),
        'postselected_fraction_is_original_weight':False})


def certificate():
    paths=provenance.definition_inputs(SOURCE_INPUTS)
    hashes=provenance.definition_hashes(paths)
    members=[contract.member(n) for n in (0,1,2,4)]
    result={'schema':'oph.primary_spectral_scale.v1',
        'actual_local_expansion_channels':vacuum_channels(),
        'selected_source_expansion_diagnostics':[selected_expansion(m) for m in members],
        'actual_fixed_band_scales':[fixed_band_scale(m) for m in members],
        'conditional_gap_controls':[conditional_gap_exclusion(m,'1') for m in members],
        'distinct_polynomial_countercontrols':[polynomial_countercontrol(n) for n in (0,2,5)],
        'primary_dimensionless_gap_proved':False,
        'primary_fixed_band_noncollapse_proved':False,
        'physical_mass_table_GeV':[], 'source_sha256':hashes}
    if provenance.definition_hashes(paths)!=hashes:
        raise ValueError('spectral-scale inputs changed during computation')
    return result
