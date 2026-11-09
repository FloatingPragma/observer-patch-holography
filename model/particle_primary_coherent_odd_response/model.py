"""Original positive-energy response for every fixed flat spatial analyzer."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product

import numpy as np
import sympy as s
from scipy.sparse import csr_matrix, eye

from model.particle_common_contract import contract, provenance
from model.particle_primary_selected_spin_response import model as spin
from . import kernel_reduction

SOURCE_INPUTS=tuple(dict.fromkeys((
 'model/particle_primary_coherent_odd_response/__init__.py',
 'model/particle_primary_coherent_odd_response/__main__.py',
 'model/particle_primary_coherent_odd_response/model.py',
 'model/particle_primary_coherent_odd_response/test_model.py',
 'model/particle_primary_coherent_odd_response/PROOF.md',
 'model/particle_primary_coherent_odd_response/README.md',
 'model/particle_primary_coherent_odd_response/kernel_reduction.py',
 'model/particle_primary_coherent_odd_response/test_kernel_reduction.py',
 'model/particle_primary_coherent_odd_response/KERNEL_PROOF.md',
 'Lean/QFT/GroundRecordTail.lean',
 *spin.SOURCE_INPUTS,
)))


@lru_cache(maxsize=1)
def local_anticommutators():
    """Exact full local CAR effects, including every occupation sector."""
    controls=spin.local_controls();rows=[];total=s.zeros(64)
    for label,family in controls['families'].items():
        for component,A in enumerate(family['operators']):
            effect=spin.clean(A.H*A+A*A.H);values=effect.eigenvals()
            if any(value<0 or value>1 for value in values):
                raise ArithmeticError('coherent signed local probe is not complete')
            total+=effect
            rows.append({'probe_irrep':label,'component':component,
                'anticommutator_spectrum':{str(k):int(v) for k,v in values.items()}})
    lower=spin.clean(total-6*s.eye(64)+3*controls['positive_number'])
    if any(value<0 for value in lower.eigenvals()):
        raise ArithmeticError('signed multiplet occupation floor failed')
    return {'complete_components':6,'rows':rows,'total_success_cap':6,
            'occupation_lower_constant':6,'occupation_lower_slope':3,
            'all_local_occupation_sectors_retained':True}


def cube_stars():
    cells=tuple(product(range(2),repeat=3))
    edges=tuple((x,y) for x in cells for y in cells if x<y and sum(a!=b for a,b in zip(x,y))==1)
    colors=tuple((x,tuple(1-v for v in x)) for x in cells if x[0]==0)
    stars={x:{'matter':{x}|{y if x==z else z for z,y in edges if x in (z,y)},
              'links':{edge for edge in edges if x in edge}} for x in cells}
    for x,y in colors:
        if stars[x]['matter']&stars[y]['matter'] or stars[x]['links']&stars[y]['links']:
            raise ArithmeticError('paired residual stars overlap')
    return {'cells':cells,'edges':edges,'colors':colors,'stars':stars}


@lru_cache(maxsize=2)
def adjoint_residual_control(label):
    """Actual twelve-CAR particle-hole conjugation certifies the adjoint cap."""
    denominator,G=spin.residual_gram(label);Gadj=csr_matrix((4096,4096),dtype=np.int64)
    for polynomial in spin.local_controls()['families'][label]['polynomials']:
        square,coefficients=spin._integer_polynomial(polynomial)
        R=csr_matrix((4096,4096),dtype=np.int64)
        for word,coefficient in coefficients.items():
            for position in range(3):
                shifted=list(word);shifted[position]+=6
                R+=coefficient*spin._creation_sparse(tuple(shifted),12)
        Gadj+=int(denominator*square)*(R@R.T)
    phases=[];targets=[]
    for mask in range(4096):
        target=mask;phase=1
        for i in reversed(range(12)):
            phase*=(-1)**((target&((1<<i)-1)).bit_count());target^=1<<i
        phases.append(phase);targets.append(target)
    J=csr_matrix((np.array(phases,dtype=np.int64),(targets,range(4096))),shape=(4096,4096))
    difference=J@G@J.T-Gadj;difference.eliminate_zeros()
    if difference.nnz or (J@J.T-eye(4096,format='csr',dtype=np.int64)).nnz:
        raise ArithmeticError('full residual particle-hole identity failed')
    inherited=spin.residual_control(label)
    return {'probe_irrep':label,'integer_Gram_denominator':denominator,
        'adjoint_Gram_is_exact_particle_hole_conjugate':True,
        'one_sign_squared_residual_cap':inherited['squared_residual_norm_upper'],
        'both_sign_squared_residual_cap':2*inherited['squared_residual_norm_upper']}


def spatial_menu_control():
    cells=cube_stars()['cells']
    W=s.Matrix([[(-1)**sum(a*b for a,b in zip(q,x)) for x in cells] for q in cells])/s.sqrt(8)
    if spin.clean(W.H*W)!=s.eye(8):raise ArithmeticError('complete Walsh analyzer is not unitary')
    return {'spatial_rows':8,'internal_components':6,'signed_success_records':96,
        'retained_address_failure_records':8,'total_records_before_clocks_and_meters':104,
        'fixed_row_success_effect_relative_to_selected_setting':F(1,8),
        'coherent_address_is_unrecorded_before_success_Walsh_rotation':True,
        'all_spatial_rows_and_original_address_failures_retained':True,
        'complex_Fourier_negative_branch_requires_conjugate_pointer_rotation':True}


def bind_member(member):
    contract.verify_member(member)
    if member['source']['level']!=0:raise ValueError('four-color cube response is certified only at member zero')
    inherited=spin.bind_member(member);spin.verify_bound(member,inherited)
    original=inherited['result'];vacancy=F(original['positive_beta_defect_vacancy_lower'])
    stars=cube_stars();local_anticommutators()
    caps=sum((adjoint_residual_control(k)['both_sign_squared_residual_cap'] for k in ('G1','H')),F())
    degree=F(2*len(stars['edges']),len(stars['cells']))
    moment=F(len(stars['colors']),6)*caps*degree/2
    if moment!=F(110,3):raise ArithmeticError('actual cube residual moment changed')
    radius=F(17,2);negative_cut=F(1,2)
    signed_union=vacancy-moment/radius**2
    downward=F(member['state']['target_ground_excess'])/negative_cut
    positive=signed_union-downward
    plan=spin.odd.spectroscopy.clock_plan(member,F(1,8),'1/16384','1/16384')
    loss=2*F(plan['clock_plan']['complete_clock_tail_upper'])+F(plan['two_complete_clock_error_upper'])
    if loss>F(1,4096) or positive-loss<=F(1,5) or positive/8-loss<=F(1,40):
        raise ArithmeticError('fixed spatial analyzer positive original weight failed')
    channels=[]
    for label,dimension,C,simple in [('G1',2,F(19,3),F(1,16)),('H',4,F(12),F(9,64))]:
        channel_moment=2*C;channel_downward=F(dimension,6)*downward
        channel_positive=F(dimension,6)*vacancy-channel_moment/radius**2-channel_downward
        if channel_positive-loss<=simple:raise ArithmeticError('fixed spatial spin channel lost positive weight')
        channels.append({'probe_irrep':label,'dimension':dimension,
            'signed_effect_upper':F(dimension,6),'signed_original_success_lower':F(dimension,6)*vacancy,
            'signed_centered_second_moment_upper':channel_moment,'negative_transfer_weight_upper':channel_downward,
            'original_positive_creation_weight_lower':channel_positive,
            'reported_original_positive_creation_weight_lower':channel_positive-loss,
            'simple_reported_positive_lower':simple,
            'complete_Walsh_each_positive_spatial_row_lower':(channel_positive-loss)/8,
            'simple_complete_Walsh_each_positive_spatial_row_lower':simple/8})
    result=contract.serial({'member_id':member['member_id'],**member['ids'],
        'selected_spin_receipt_sha256':contract.digest(inherited),
        'source_input_is_original_selected_density':True,'flat_spatial_amplitude_squared':F(1,8),
        'all_flat_complex_phases_allowed':True,'chosen_analyzer_setting_is_not_randomly_postselected':True,
        'complete_signed_success_outcomes':12,'complete_failure_outcomes':1,
        'all_signed_failure_memory_reference_and_clock_outputs_retained':True,
        'creation_native_increment':spin.DELTA,'annihilation_native_increment':tuple(-x for x in spin.DELTA),
        'coherent_native_current_meter':spin.current_meter_control(),
        'color_classes':stars['colors'],'residual_star_color_count':4,
        'product_link_cutoff_preserves_odd_star_support':True,
        'Gauss_projection_is_restriction_after_gauge_invariant_support_identity':True,
        'original_total_signed_success_lower':vacancy,
        'creation_moment_center':9,'annihilation_moment_center':-9,
        'complete_signed_centered_second_moment_upper':moment,
        'centered_signed_union_lower':signed_union,'all_negative_transfer_at_least_half_weight_upper':downward,
        'original_positive_creation_interval':[F(1,2),F(35,2)],
        'original_positive_creation_weight_lower':positive,
        'finite_reported_creation_interval':[F(1,4),F(71,4)],
        'finite_reported_original_creation_weight_lower':positive-loss,
        'simple_reported_original_creation_weight_lower':F(1,5),
        'source_spatial_kernel_reduction':kernel_reduction.control(),
        'complete_Walsh_menu':{**spatial_menu_control(),
            'each_positive_spatial_row_exact_original_weight_lower':positive/8,
            'each_reported_positive_spatial_row_weight_lower':(positive-loss)/8,
            'each_reported_positive_spatial_row_coarse_global_error_lower':positive/8-loss,
            'simple_reported_each_positive_spatial_row_lower':F(1,40),
            'sum_of_eight_simple_disjoint_event_lowers':F(1,5)},
        'complete_two_clock_probability_error_upper':loss,'complete_clock_plan':plan,
        'source_conjugate_apparatus_has_identical_complete_record_law':True,
        'source_conjugation_implies_equal_retained_reference_output_states':False,
        'spin_resolved_channels':channels,
        'each_G1_and_H_multiplet_has_positive_fixed_spatial_weight_proved':True,
        'each_individual_internal_component_has_positive_fixed_spatial_weight_proved':False,
        'fixed_spatial_analyzer_is_conserved_momentum':False,
        'isolated_particle_pole_or_uniform_refinement_survival_proved':False,
        'physical_mass_table_GeV':[]})
    return {'schema':'oph.primary_coherent_odd_response.v1','result':result,'result_sha256':contract.digest(result)}


def verify_bound(member,value):
    if value!=bind_member(member):raise ValueError('coherent odd response or source binding changed')


def certificate():
    paths=provenance.definition_inputs(SOURCE_INPUTS);pins=provenance.definition_hashes(paths)
    value={'schema':'oph.primary_coherent_odd_response_certificate.v1',
        'exact_local_anticommutators':local_anticommutators(),
        'exact_adjoint_residual_controls':[adjoint_residual_control(k) for k in ('G1','H')],
        'members':[bind_member(contract.member(0))],'source_sha256':pins}
    if provenance.definition_hashes(paths)!=pins:raise ValueError('source changed during coherent odd replay')
    return contract.serial(value)
