"""Local electric Schur bounds at the actual primary coupling.

The result encloses the full physical ground. It does not use a fitted
perturbation coefficient or normalize a rare apparatus branch.
"""
from fractions import Fraction as F
from functools import lru_cache

import sympy as s

from model.particle_common_contract import contract, provenance
from model.particle_primary_spectral_scale import scale
from model.particle_spatial_ground_state import ground
from model.particle_local_wilson_density import density
from model.spin_packet_readout.bounds import sqrt_upper

SOURCE_INPUTS=(
    'model/particle_primary_endpoint_bounds/__init__.py',
    'model/particle_primary_endpoint_bounds/__main__.py',
    'model/particle_primary_endpoint_bounds/endpoint.py',
    'model/particle_primary_endpoint_bounds/test_endpoint.py',
    'model/particle_primary_endpoint_bounds/PROOF.md',
    'model/particle_primary_spectral_scale/PROOF.md',
    'model/particle_common_contract/PROOF.md',
    'model/particle_primary_locality/PROOF.md',
    'model/particle_primary_survival/PROOF.md',
    'model/particle_primary_variational_readout/PROOF.md',
)


@lru_cache(maxsize=1)
def actual_spin_polynomial():
    """Literal eight-dimensional Dirac matrix, before the native lift."""
    beta=density.gamma_data()[3]
    m,z=s.symbols('m z',real=True)
    expected=(z**4-(2*m*m+1)*z*z+m**4)**2
    rows=[]
    for axis,gamma in enumerate(density.gamma_data()[:3]):
        hop=(-beta+beta*gamma)/2
        h=(m*beta).row_join(hop).col_join(hop.H.row_join(m*beta))
        characteristic=h.charpoly(z)
        polynomial=s.expand(characteristic.as_expr().subs(characteristic.gen,z))
        if s.expand(polynomial-expected)!=0:
            raise ArithmeticError('the actual two-site Wilson characteristic polynomial differs')
        rows.append({'axis':axis,'polynomial':str(expected),
                     'all_eight_Dirac_coordinates_retained':True})
    return {'one_internal_column_checks':rows,
            'full_Fock_lower_shift_formula':'60*(sqrt(m^2+1/4)-m)',
            'native_internal_dimension':15,
            'coordinate_link_removed_only_by_pointwise_unitary_equivalence':True,
            'quantum_link_or_electric_action_removed':False}


def haar_variance_max(dimension, occupation_x, occupation_y):
    """Exact maximal spin-coupled eigenvalue for supplied occupations.

    Allowing maximal spins gives an upper bound even when a particular
    Fock representation lacks that spin. The global maximum is attained.
    """
    d=ground._i(dimension,'native irrep dimension',1)
    x=ground._i(occupation_x,'left active occupation',0)
    y=ground._i(occupation_y,'right active occupation',0)
    if max(x,y)>2*d:raise ValueError('two active Dirac rows permit at most 2d occupations')
    sx=F(min(x,2*d-x),2);sy=F(min(y,2*d-y),2)
    return F(x+y)-F(x*y,d)+F(4,d)*min(sx,sy)*(max(sx,sy)+1)


def local_schur_data(member,retained_onsite_fraction='0'):
    """Keep r H0 explicitly and bound the rest by actual one-link clusters."""
    contract.verify_member(member)
    source=member['source'];r=ground._q(retained_onsite_fraction,'retained H0 fraction')
    if r>=1:raise ValueError('retained H0 fraction must be less than one')
    t=1-r
    maximum_degree=3 if source['side']==2 else 6
    m=t*F(3,maximum_degree)
    electric=t*F(source['electric_coefficient'])
    root_upper=sqrt_upper(m*m+F(1,4))
    root_lower=(m*m+F(1,4))/root_upper
    loss=60*(root_upper-m)
    loss_lower=60*(root_lower-m)
    rows=[]
    for row in scale.vacuum_channels()['all_native_blocks']:
        c=F(row['trace_metric_Casimir']);dimension=row['dimension']
        margin=electric*c/2-loss
        rows.append({'integer_Q6_charge':row['integer_Q6_charge'],
            'dimension':dimension,'trace_metric_Casimir':c,
            'constant_link_Haar_variance_norm_upper':2*(dimension+1),
            'electric_minus_fermion_loss_lower':margin,
            'Schur_loss_upper':F(2*(dimension+1))/margin if margin>0 else None})
    minimum=electric*F(3,20)-loss
    admitted=minimum>0
    return {'retained_H0_fraction':r,'cluster_fraction':t,
        'maximum_interior_vertex_degree':maximum_degree,
        'onsite_coefficient_per_endpoint':m,
        'cluster_electric_coefficient':electric,
        'full_Fock_pointwise_lower_shift_upper':loss,
        'full_Fock_pointwise_lower_shift_lower':loss_lower,
        'all_word_nonconstant_link_Schur_margin_lower':minimum,
        'nonconstant_link_block_positive':admitted,
        'native_channels':rows,
        'per_bond_Schur_loss_upper':sum((row['Schur_loss_upper'] for row in rows),F()) if admitted else None,
        'per_bond_exact_Schur_loss_lower':sum((F(2*(row['dimension']+1))/(electric*F(row['trace_metric_Casimir'])/2-loss_lower)
                                            for row in rows),F()) if admitted else None,
        'pointwise_fermion_bound_is_full_Fock':True,
        'actual_product_link_cutoff_preserves_bound':True}


def retained_penalty_bound(member,retained_onsite_fraction='0'):
    """A full physical ground floor and an ordered second-eigenvalue floor."""
    local=local_schur_data(member,retained_onsite_fraction)
    source=member['source'];a=F(source['spacing']);volume=source['side']**3
    anchor=-90*volume;magnetic=F(source['magnetic_coefficient'])*source['interior_plaquettes']
    if local['nonconstant_link_block_positive']:
        loss=source['interior_hopping_bonds']*local['per_bond_Schur_loss_upper']+magnetic
        lower=(anchor-loss)/a
        second=lower+6*local['retained_H0_fraction']/a
    else:loss=lower=second=None
    return contract.serial({'member_id':member['member_id'],**member['ids'],
        'local_bound':local,'dimensionless_total_loss_upper':loss,
        'physical_ground_lower':lower,
        'physical_second_eigenvalue_lower':second,
        'second_eigenvalue_counts_multiplicity':True,
        'primary_endpoint_lambda':'1','full_physical_carrier':True,
        'anchor_subspace_used_as_invariant':False,
        'branch_isolation_proved':False,
        'selected_polynomial_state_changed':False,
        'source_reference_or_clock_changed':False})


def selected_bound(member):
    """Default exact rational floor, with source-bound scope and comparison."""
    result=retained_penalty_bound(member,'0')
    if result['physical_ground_lower'] is None:
        raise ArithmeticError('the actual primary parameters fail the local electric bound')
    source=member['source'];a=F(source['spacing']);volume=source['side']**3
    inherited=F(member['preparation_plan']['source']['spectral_lower'])
    lower=F(result['physical_ground_lower'])
    loss_lower=F(result['local_bound']['per_bond_exact_Schur_loss_lower'])*source['interior_hopping_bonds']
    loss_lower+=F(source['magnetic_coefficient'])*source['interior_plaquettes']
    method_ceiling=(F(-90*volume)-loss_lower+6)/a
    if lower<inherited:raise ArithmeticError('the local floor does not improve the inherited source floor')
    result.update({
        'schema':'oph.primary_endpoint_ground_floor.v1',
        'physical_ground_upper_from_product_trial':str(F(-90*volume)/a),
        'inherited_full_source_lower':str(inherited),
        'ground_bound_improves_inherited_source_floor':True,
        'bound_route':'actual one-link electric Schur inequality, summed before physical restriction',
        'upper_ceiling_on_this_retained_penalty_method_second_eigenvalue_floors':str(method_ceiling),
        'method_ceiling_is_not_an_upper_bound_on_actual_second_eigenvalue':True,
        'ground_energy_is_not_a_particle_rest_energy':True,
        'physical_mass_table_GeV':[]})
    return result


def bind_member(member):
    result=selected_bound(member)
    return {'schema':'oph.primary_endpoint_bound.v1','result':result,
            'result_sha256':contract.digest(result)}


def verify_bound(member,value):
    if value!=bind_member(member):
        raise ValueError('endpoint ground floor or its actual source changed')


def compare_with_variational(member,upper_receipt):
    """Join independently authenticated full-source bounds, without a pole claim."""
    from model.particle_primary_variational_readout import variational
    variational.verify_bound(member,upper_receipt)
    lower=bind_member(member)
    floor=F(lower['result']['physical_ground_lower'])
    upper=F(upper_receipt['result']['ground_energy_upper'])
    ceiling=F(lower['result']['upper_ceiling_on_this_retained_penalty_method_second_eigenvalue_floors'])
    if floor>upper:raise ArithmeticError('the independently proved ground enclosure is inconsistent')
    trials=[]
    for r in ('0','1/10','1/4','1/2'):
        value=retained_penalty_bound(member,r)
        second=value['physical_second_eigenvalue_lower']
        trials.append({'retained_H0_fraction':r,
            'Schur_Q_block_admitted':value['local_bound']['nonconstant_link_block_positive'],
            'certified_second_eigenvalue_lower':second,
            'ground_isolation_margin':None if second is None else F(second)-upper})
    return contract.serial({'member_id':member['member_id'],**member['ids'],
        'full_physical_ground_energy_interval':[floor,upper],
        'lower_result_sha256':lower['result_sha256'],
        'upper_result_sha256':upper_receipt['result_sha256'],
        'second_level_method_ceiling':ceiling,
        'ceiling_minus_actual_variational_upper':ceiling-upper,
        'entire_retained_penalty_family_fails_with_this_upper':ceiling<=upper,
        'retained_penalty_trials':trials,
        'actual_full_source_gap_ruled_out':False,
        'isolated_branch_or_particle_mass_proved':False})


def certificate():
    from model.particle_primary_variational_readout import variational
    paths=provenance.definition_inputs((*SOURCE_INPUTS,*variational.SOURCE_INPUTS))
    pins=provenance.definition_hashes(paths)
    members=[contract.member(n) for n in (0,1)]
    result={'schema':'oph.primary_endpoint_certificate.v1',
        'actual_spin_polynomial':actual_spin_polynomial(),
        'source_ground_bounds':[bind_member(member) for member in members],
        'retained_penalty_trials':[retained_penalty_bound(m,r)
            for m in members for r in ('0','1/10','1/4','1/2')],
        'authenticated_endpoint_enclosures':[compare_with_variational(m,variational.bind_member(m))
                                             for m in members],
        'endpoint_branch_isolation_proved':False,
        'physical_mass_table_GeV':[], 'source_sha256':pins}
    if provenance.definition_hashes(paths)!=pins:
        raise ValueError('endpoint-bound inputs changed during replay')
    return result
