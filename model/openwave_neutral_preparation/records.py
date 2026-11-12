"""Complete preparation, momentum, energy and native occupation instruments."""
import sympy as s
from model.openwave_neutral_stationary.isolation import data
from model.openwave_neutral_stationary.local_profile import certificate as local_profile
from .filtering import parameters


def certificate():
    p=parameters();d=data();profile=local_profile();eta2=d['eta']**2
    q,h,w=s.symbols('q h w');weights=s.symbols('p0:6')
    per_recipe=[1-q]
    for _ in range(8):
        per_recipe.extend((q*(1-h)/8,q*h*(1-8*w)/8))
        per_recipe.extend(q*h*w*z/8 for _ in range(8) for z in weights)
    assert len(per_recipe)==401
    assert s.expand(sum(per_recipe).subs(weights[-1],1-sum(weights[:-1])))==1
    leak=(p['epsilon']/p['q0'])**2
    return {'preparation_recipe_count':1024,'preparation_recipe_weight':'1/1024',
        'complete_records_per_recipe':401,'complete_original_record_count':1024*401,
        'factorized_record_law':{
            'recipe_and_failure':'(1-q_recipe)/1024',
            'recipe_momentum_energy_complement':'q_recipe (1-h_recipe)/(1024*8)',
            'recipe_momentum_band_profile_complement':'q_recipe h_recipe (1-8w_recipe)/(1024*8)',
            'recipe_momentum_band_cell_occupation':'q_recipe h_recipe w_recipe p_alpha/(1024*8)'},
        'complete_probability_sum':'1 exactly for every recipe, then1 after the uniform mixture',
        'all_eight_momentum_probabilities_conditioned_on_each_recipe_success':'1/8 exactly, by its initial C-eigenvalues and [K,C_i]=0',
        'recipe_success_probability_interval':[str(p['success_lower']),str(p['q0'])],
        'conditional_energy_complement_upper':str(leak),
        'each_band_cell_probability_interval':[str((1-eta2)/8),'1/8'],
        'complete_six_native_local_probabilities':profile['complete_gauge_invariant_native_occupation_records'],
        'measurement_order':'recipe archive retained; complete K accept/fail; all eight momenta; full-H band/complement; eight P_y/profile complement; six local occupation effects',
        'position_and_momentum_joint_sharp_values_claimed':False,
        'archive_not_discarded_to_prove_channel_completeness':True,
        'uniform_reduced_density_stationarity_requires_marginalizing_recipe_archive':True,
        'arbitrary_purification_or_recipe_conditional_quantum_state_stationary':False,
        'exact_source_energy_reference_and_calibration':'unchanged isolated neutral band and actual κ2 ground; no additional empirical input',
        'stationarity_and_preparation_cost_are_distinct_from_physical_particle_identification':True}
