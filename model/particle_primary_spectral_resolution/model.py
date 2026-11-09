"""Full residual bounds and a complete retained-memory preparation experiment."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product

import sympy as s

from model.particle_common_contract import contract, provenance
from model.particle_primary_odd_residual import model as odd
from model.particle_primary_joint_spectroscopy import spectroscopy
from model.particle_primary_capacity import capacity
from model.spin_packet_readout.bounds import sqrt_upper

SOURCE_INPUTS=(
    'model/particle_primary_spectral_resolution/__init__.py',
    'model/particle_primary_spectral_resolution/__main__.py',
    'model/particle_primary_spectral_resolution/model.py',
    'model/particle_primary_spectral_resolution/test_model.py',
    'model/particle_primary_spectral_resolution/PROOF.md',
    'model/particle_primary_spectral_resolution/README.md',
    'model/particle_common_contract/PROOF.md',
    'model/particle_primary_odd_residual/PROOF.md',
    'model/particle_primary_variational_readout/PROOF.md',
    'model/particle_primary_joint_spectroscopy/PROOF.md',
    'model/particle_primary_capacity/PROOF.md',
    'model/particle_transport/algebra.py',
)

# These coefficients define a physical trial/intervention, never the action.
# They were selected from the source residual quadratic form without mass data.
# Global optimality of this rational choice is not asserted.
HOPPING_WEIGHTS=tuple(F(k,10**6) for k in (17264,24063,14217,18842,31291,28137))
LOOP_WEIGHTS=tuple(F(k,10**6) for k in (5697,9729,4390,6466,13859))


def residual_form(trial,hopping_weights,loop_weights):
    """Full parent-action bound for arbitrary exact rational packet weights."""
    th=tuple(map(F,hopping_weights));tl=tuple(map(F,loop_weights))
    if len(th)!=6 or len(tl)!=5:raise ValueError('six hopping and five loop coefficients required')
    rows=trial['hopping_packets'];n=tuple(F(r['norm_squared']) for r in rows)
    d=tuple(F(r['relative_H0_energy']) for r in rows)
    g=F(trial['loop_packet_norm_squared'])
    A=tuple(tuple(F(v) for v in row) for row in trial['loop_packet_matrix'])
    correction=odd.variational.character_control()['cubic_loop_singlet_matrix']
    C=tuple(tuple(-F(v,420) for v in row) for row in correction)
    moment=odd.odd_parent_moment(th,tl)
    square=sum((v*t*t for v,t in zip(n,th)),F())+g*sum((t*t for t in tl),F())
    q=sum((v*t for v,t in zip(n,th)),F())+g*sum(tl)
    e=(-2*q+sum((v*h*t*t for v,h,t in zip(n,d,th)),F())
       +g*sum((tl[i]*A[i][j]*tl[j] for i,j in product(range(5),repeat=2)),F()))/(1+square)
    correction_vector=[sum((C[i][j]*tl[j] for j in range(5)),F()) for i in range(5)]
    outside=moment['exact_parent_second_action_norm_squared']-q*q-g*sum((x*x for x in correction_vector),F())
    hop_residual=sum((v*(1-(h-e)*t)**2 for v,h,t in zip(n,d,th)),F())
    loop_residual=g*sum(((1-sum((A[i][j]*tl[j] for j in range(5)),F())+e*tl[i])**2 for i in range(5)),F())
    inside=(q+e)**2+hop_residual+loop_residual
    residual2=(inside+outside)/(1+square)
    if outside<0 or residual2<0:raise ArithmeticError('invalid full-space residual form')
    return {'hopping_coefficients':th,'loop_coefficients':tl,
        'packet_norms_squared':n+(g,)*5,'trial_norm_squared':1+square,
        'source_energy_center':-714+e,'relative_Rayleigh_energy':e,
        'full_parent_second_action_norm_squared':moment['exact_parent_second_action_norm_squared'],
        'parent_omitted_norm_squared':outside,'retained_residual_squared':inside,
        'full_residual_squared_upper':residual2,'full_residual_upper':sqrt_upper(residual2),
        'source_action_changed':False,'global_residual_optimum_proved':False,
        'selected_state_expectation_replaced_by_trial':False}


def physical_trial(member):
    odd.branch._member_zero(member)
    inherited=odd.branch.bind_member(member);odd.branch.verify_bound(member,inherited)
    trial=inherited['result']['odd_endpoint_trial']
    result=residual_form(trial,HOPPING_WEIGHTS,LOOP_WEIGHTS)
    center=result['source_energy_center'];radius=result['full_residual_upper']
    lower=F(trial['full_physical_ground_lower']);upper=F(trial['ground_Ritz_upper_for_comparison'])
    if result['full_residual_squared_upper']>=F(49,2) or center-radius<=upper:
        raise ArithmeticError('the prepared scalar-odd trial lost its excited enclosure')
    result.update({'member_id':member['member_id'],**member['ids'],
        'inherited_packet_receipt_sha256':contract.digest(inherited),
        'physical_cubic_irrep':'A1','physical_inversion_parity':-1,'CAR_parity':'even',
        'physical_blank_orthogonal_to_trial':True,'physical_packet_dimension':12,
        'full_ground_lower':lower,'full_ground_upper':upper,
        'some_physical_A1u_energy_interval':[center-radius,center+radius],
        'some_physical_A1u_excitation_interval':[center-radius-upper,center+radius-lower],
        'trial_is_selected_polynomial_state':False,'isolated_eigenvalue_proved':False})
    return result


def prepared_records(member):
    trial=physical_trial(member);center=trial['source_energy_center'];radius=F(23,4)
    variance=trial['full_residual_squared_upper'];weight=1-variance/radius**2
    cap=capacity.member_capacity(member)
    plan=spectroscopy.clock_plan(member,'1/64','1/1024','1/1024')
    loss=F(plan['clock_plan']['complete_clock_tail_upper'])+F(plan['one_complete_clock_error_upper'])
    eta=F(plan['energy_resolution']);reported=[center-radius-eta,center+radius+eta]
    if weight<F(137,529) or loss>F(1,512) or weight-loss<=F(1,4):
        raise ArithmeticError('prepared complete energy-record floor failed')
    if reported[0]<=trial['full_ground_upper']:
        raise ArithmeticError('the reported interval overlaps the ground enclosure')
    eps=F(member['state']['target_ground_excess'])
    return contract.serial({'member_id':member['member_id'],**member['ids'],
        'physical_trial':trial,
        'initial_state_is_selected_primary_preparation':True,
        'apparatus_sequence':['prepare trial on identical idle physical memory with exact Householder',
                              'graded SWAP of source and prepared memory',
                              'complete one-clock energy readout on unchanged source Hamiltonian'],
        'post_intervention_source_state':'pure source-derived scalar-odd trial',
        'memory_retains_displaced_source_and_its_original_reference_correlations':True,
        'source_marginal_preserved_by_preparation_intervention':False,
        'original_reference_discarded':False,'postselection_or_success_normalization':False,
        'preparation_outcome_count':1,'complete_energy_outcomes':plan['clock_plan']['clock_size'],
        'source_Gauss_constraints_preserved':True,'total_fermion_parity_preserved':True,
        'idle_memory_H':0,'additional_memory_physical_capacity_upper':cap['physical_capacity_upper'],
        'additional_memory_qubit_bound':cap['allocated_source_qubit_count_upper'],
        'exact_preparation_primitive':'algebraic Householder in the admitted physical matrix service',
        'source_work_mean_interval':[center-trial['full_ground_upper']-eps,center-trial['full_ground_lower']],
        'autonomous_energy_conserving_controller_derived':False,
        'original_prepared_spectral_window':[center-radius,center+radius],
        'original_prepared_spectral_weight_lower':weight,
        'original_prepared_spectral_weight_simple_lower':F(137,529),
        'original_prepared_excitation_window':[center-radius-trial['full_ground_upper'],
                                               center+radius-trial['full_ground_lower']],
        'complete_clock_plan':plan,'one_clock_event_loss_upper':loss,
        'reported_absolute_source_energy_interval':reported,
        'reported_original_probability_lower':weight-loss,
        'reported_original_probability_simple_lower':F(1,4),
        'reported_interval_lower_minus_ground_upper':reported[0]-trial['full_ground_upper'],
        'finite_source_spectroscopy_preserves_A1u_support':True,
        'all_clock_outcomes_source_memory_and_reference_retained':True,
        'extra_preparation_or_swap_compiler_error_included':False,
        'extra_error_rule':'subtract the separately bounded complete preparation/SWAP output error',
        'full_source_matrix_preparation_tape_or_histogram_executed':False,
        'this_is_the_selected_K00_response':False,
        'selected_K00_residue_attached_to_trial':False,
        'physical_species_or_relativistic_rest_mass_proved':False,'physical_mass_table_GeV':[]})


@lru_cache(maxsize=1)
def swap_reference_control():
    """Exact all-sector graded exchange with a nonclassical retained reference."""
    d=4;parity=(0,0,1,1);target=s.Matrix([s.Rational(3,5),s.Rational(4,5),0,0])
    blank=s.eye(d)[:,0];v=blank-target;U=s.eye(d)-2*v*v.H/(v.H*v)[0]
    swap=s.zeros(d*d)
    for i,j in product(range(d),repeat=2):swap[j*d+i,i*d+j]=(-1)**(parity[i]*parity[j])
    source_ref=s.zeros(d*2,1);source_ref[0]=1/s.sqrt(2);source_ref[5]=s.I/s.sqrt(2)
    initial=s.zeros(d*d*2,1)
    for i,j,r in product(range(d),range(d),range(2)):initial[(i*d+j)*2+r]=source_ref[i*2+r]*blank[j]
    gate=s.kronecker_product(swap,s.eye(2))*s.kronecker_product(s.eye(d),U,s.eye(2))
    final=gate*initial;expected=s.zeros(d*d*2,1)
    for i,j,r in product(range(d),range(d),range(2)):expected[(i*d+j)*2+r]=target[i]*source_ref[j*2+r]
    if final!=expected or gate.H*final!=initial or U.H*U!=s.eye(d):
        raise ArithmeticError('graded preparation lost the retained source-reference state')
    ordinary=s.zeros(d*d)
    for i,j in product(range(d),repeat=2):ordinary[j*d+i,i*d+j]=1
    if ordinary==swap:raise ArithmeticError('graded exchange was replaced by ordinary exchange')
    return {'control_source_dimension':4,'control_reference_dimension':2,
        'exact_recovery_of_entangled_source_reference':True,
        'source_target_CAR_parity':'even','graded_exchange_differs_on_odd_odd_sector':True,
        'primary_spatial_histogram_executed':False}


def bind_member(member):
    swap_reference_control();result=prepared_records(member)
    return {'schema':'oph.primary_spectral_resolution.v1','result':result,'result_sha256':contract.digest(result)}


def verify_bound(member,receipt):
    if receipt!=bind_member(member):raise ValueError('prepared spectral result or physical source changed')


def certificate():
    paths=provenance.definition_inputs((*SOURCE_INPUTS,*odd.SOURCE_INPUTS,*spectroscopy.SOURCE_INPUTS))
    pins=provenance.definition_hashes(paths)
    result={'schema':'oph.primary_spectral_resolution_certificate.v1',
        'selected_member_bound':bind_member(contract.member(0)),
        'graded_reference_control':swap_reference_control(),'source_sha256':pins}
    if provenance.definition_hashes(paths)!=pins:raise ValueError('spectral inputs changed during replay')
    return contract.serial(result)
