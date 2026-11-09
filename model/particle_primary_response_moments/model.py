"""Source-bound original spectral weight from local CAR commutator moments."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product

import sympy as s

from model.particle_common_contract import contract, provenance
from model.particle_local_wilson_density import density
from model.particle_primary_channel_labels import channels
from model.particle_primary_occupation_bound import model as occupation
from model.particle_primary_joint_spectroscopy import spectroscopy
from model.particle_primary_mass_readout import readout
from model.particle_filtered_probe import filter as apparatus

SOURCE_INPUTS = (
    'model/particle_primary_response_moments/__init__.py',
    'model/particle_primary_response_moments/__main__.py',
    'model/particle_primary_response_moments/model.py',
    'model/particle_primary_response_moments/test_model.py',
    'model/particle_primary_response_moments/PROOF.md',
    'model/particle_primary_response_moments/README.md',
    'model/particle_primary_occupation_bound/PROOF.md',
    'model/particle_primary_joint_spectroscopy/PROOF.md',
    'model/particle_primary_mass_readout/PROOF.md',
    'model/particle_primary_channel_labels/PROOF.md',
    'model/particle_common_contract/PROOF.md',
)


def beta_pairs():
    """Two declared matched rows in the actual source Dirac chart."""
    return tuple((s.Matrix([int(j==k)+int(j==k+2) for j in range(4)])/s.sqrt(2),
                  s.Matrix([-int(j==k)+int(j==k+2) for j in range(4)])/s.sqrt(2))
                 for k in range(2))


@lru_cache(maxsize=1)
def source_probe_control():
    """Authenticate native q6 dimension and actual corner Wilson Gram matrices."""
    blocks=channels.native_number_control()['blocks']
    if next(row['internal_dimension'] for row in blocks if row['Q6']==6)!=1:
        raise ArithmeticError('the selected native block is not one dimensional')
    gamma=density.gamma_data();beta=gamma[3]
    spin=tuple(apparatus.clean((-beta+beta*g)/2) for g in gamma[:3])
    for p,m in beta_pairs():
        if beta*p!=p or beta*m!=-m or (p.H*m)[0]!=0:
            raise ArithmeticError('the declared beta probe chart changed')
    rows=[]
    for cell in product(range(2),repeat=3):
        # H_xy=S on each positively oriented source bond, H_yx=S†.
        # Its column at x is S† for the outgoing neighbor, S for incoming.
        local=sum((spin[j]*spin[j].H if cell[j]==0 else spin[j].H*spin[j]
                   for j in range(3)),s.zeros(4))
        for alpha,(p,m) in enumerate(beta_pairs()):
            hp=s.simplify((p.H*local*p)[0]);hm=s.simplify((m.H*local*m)[0])
            c=s.simplify((p.H*local*m)[0]);square=s.simplify(s.conjugate(c)*c)
            if hp!=s.Rational(3,2) or hm!=hp or square!=s.Rational(1,4):
                raise ArithmeticError('the actual corner commutator Gram data changed')
            rows.append({'cell':list(cell),'matched_beta_row':alpha,
                         'neighbor_vector_norm_squared':'3/2','neighbor_overlap_squared':'1/4'})
    return {'native_integer_Q6':6,'native_dimension':1,'matched_rows_per_cell':2,
            'actual_corner_Gram_controls':rows,
            'all_Fock_exact_corner_residual_squared':'3/2+sqrt(2)',
            'all_Fock_corner_residual_squared_rational_upper':'3',
            'five_native_occupation_counts_preserved':True,
            'gauge_invariant':True,'CAR_even':True,
            'individual_probe_is_cubic_scalar_proved':False}


def averaged_instrument(probes):
    """Small exact control of the complete address-plus-failure service."""
    probes=tuple(apparatus.exact_matrix(b) for b in probes)
    if not probes or len(probes)>16:raise ValueError('one through sixteen probes required')
    dimension=probes[0].rows
    if dimension>16 or any(b.shape!=(dimension,dimension) for b in probes):
        raise ValueError('small square probes on one common carrier required')
    count=len(probes)
    success=[apparatus.clean(b/s.sqrt(count)) for b in probes]
    effect=apparatus.clean(sum((b.H*b for b in success),s.zeros(dimension)))
    failure=apparatus.positive_sqrt(s.eye(dimension)-effect)
    rows=tuple(success+[failure])
    if apparatus.clean(sum((b.H*b for b in rows),s.zeros(dimension)))!=s.eye(dimension):
        raise ArithmeticError('complete averaged probe lost probability')
    return rows


def selected_response(member):
    """Regenerate the actual occupation producer before issuing a response bound."""
    contract.verify_member(member)
    if member['source']['side']!=2 or F(member['source']['spacing'])!=1:
        raise ValueError('the selected-state response certificate is proved for member zero only')
    probe=source_probe_control()
    preparation=occupation.bind_member(member)
    occupation.verify_bound(member,preparation)
    inherited=preparation['result']
    defects=F(inherited['defect_count_upper']);weight=F(inherited['average_matched_q6_pair_weight_lower'])
    if defects>=5 or weight<F(11,16):
        raise ArithmeticError('the authenticated selected state lacks the required pair weight')
    moment=F(3);tail=moment/16;band=weight-tail
    epsilon=F(member['state']['target_ground_excess'])
    final=readout.excitation_window((2,10),band,epsilon,F(1,2))
    simple_final=readout.excitation_window((2,10),F(1,2),epsilon,F(1,2))
    clock=spectroscopy.clock_plan(member,'1/16','1/128','1/128')
    clock_tail=2*F(clock['clock_plan']['complete_clock_tail_upper'])
    clock_error=F(clock['two_complete_clock_error_upper'])
    loss=clock_tail+clock_error
    if loss>F(1,32):raise ArithmeticError('the complete two-clock error exceeds the public allowance')
    J=clock['clock_plan']['clock_size']
    return contract.serial({'member_id':member['member_id'],**member['ids'],
        'source_coupling_endpoint':1,'selected_input_state_unchanged':True,
        'source_probe':probe,'probe_address_count':16,'complete_probe_outcomes':17,
        'occupation_producer_receipt_sha256':contract.digest(preparation),
        'occupation_producer_recomputed':True,'selected_total_beta_defect_upper':defects,
        'original_success_weight_lower':weight,'original_success_weight_simple_lower':F(11,16),
        'source_transfer_center':6,'source_transfer_centered_second_moment_upper':moment,
        'original_transfer_interval':[2,10],'outside_transfer_interval_weight_upper':tail,
        'original_transfer_weight_lower':band,'original_transfer_weight_simple_lower':F(1,2),
        'initial_excitation_cut':F(1,2),'selected_ground_excess_upper':epsilon,
        'original_final_excitation':final,
        'original_final_excitation_weight_simple_lower':simple_final['original_final_excitation_weight_lower'],
        'finite_clock_plan':clock,'finite_reported_transfer_interval':[F(15,8),F(81,8)],
        'two_energy_tail_allowance':clock_tail,'complete_two_clock_error_upper':clock_error,
        'complete_clock_and_tail_loss_upper':loss,
        'finite_reported_original_event_weight_lower':max(F(),band-loss),
        'finite_reported_original_event_weight_simple_lower':F(15,32),
        'complete_fine_history_outcomes':17*J*J,
        'probe_service':'exact aggregate address Kraus operators and complementary occupation root',
        'all_failures_and_reference_retained':True,'reference_H_is_zero':True,
        'initial_energy_record_preserves_source_marginal':True,
        'initial_energy_record_preserves_joint_purification':False,
        'additional_probe_compiler_or_sampling_error_included':False,
        'additional_event_error_rule':'subtract any separately certified complete or sampling event error',
        'actual_primary_histogram_executed':False,'efficient_compiler_proved':False,
        'original_weight_conditioned_on_success':False,
        'blank_or_trial_substituted_for_selected_state':False,
        'numerical_ground_energy_subtraction_required':False,
        'individual_address_weight_or_spatial_invariance_proved':False,
        'isolated_eigenbranch_rest_frame_or_species_proved':False,
        'uniform_refinement_survival_proved':False,'physical_mass_table_GeV':[]})


def bind_member(member):
    result=selected_response(member)
    return {'schema':'oph.primary_response_moments.v1','result':result,
            'result_sha256':contract.digest(result)}


def verify_bound(member,value):
    if value!=bind_member(member):raise ValueError('response bound or actual source contract changed')


def certificate():
    paths=provenance.definition_inputs((*SOURCE_INPUTS,*occupation.SOURCE_INPUTS,*spectroscopy.SOURCE_INPUTS))
    if any(not (contract.ROOT/path).is_file() for path in paths):
        raise ValueError('a response evidence path is absent')
    pins=provenance.definition_hashes(paths)
    result={'schema':'oph.primary_response_moments_certificate.v1',
            'selected_member_bound':bind_member(contract.member(0)),
            'original_selected_state_transfer_band_proved':True,
            'resolved_particle_branch_or_mass_proved':False,
            'physical_mass_table_GeV':[],'source_sha256':pins}
    if provenance.definition_hashes(paths)!=pins:raise ValueError('response inputs changed during replay')
    return result
