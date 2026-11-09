"""Exact retained-origin plans on the unchanged finite primary source."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product

import sympy as s

from model.particle_common_contract import contract, provenance
from model.particle_primary_spatial_algebra import spatial
from model.particle_primary_locality import locality
from model.particle_primary_refinement_response import model as response
from model.particle_spatial_ground_state import ground

SOURCE_INPUTS=(
    'model/particle_primary_origin_momentum/__init__.py',
    'model/particle_primary_origin_momentum/__main__.py',
    'model/particle_primary_origin_momentum/model.py',
    'model/particle_primary_origin_momentum/test_model.py',
    'model/particle_primary_origin_momentum/PROOF.md',
    'model/particle_primary_origin_momentum/README.md',
    'model/particle_primary_spatial_algebra/PROOF.md',
    'model/particle_primary_field_algebra/PROOF.md',
    'model/particle_spatial_ground_state/PROOF.md',
    'model/particle_primary_refinement_response/PROOF.md',
    'model/particle_joint_reconstruction/POSITIVE_ENERGY.md',
    'Lean/QFT/GroundRecordTail.lean',
    'Lean/QFT/ColimitSelectedStateGNS.lean',
    'Lean/QFT/PeriodicMomentum.lean',
)


def origin_geometry(level):
    """A nearly full origin cube with a physically diverging boundary collar."""
    n=ground._i(level,'nonnegative source level',0)
    L=2**(2*n+1);a=F(1,2**n)
    b=min((L-1)//2,2**(3*n//2));q=L-2*b
    fraction=F(q**3,L**3)
    return {'source_level':n,'source_side':L,'source_spacing':a,
        'collar_cells':b,'physical_collar':a*b,'origin_side':q,
        'origin_lower_corner':[b]*3,'origin_upper_exclusive':[L-b]*3,
        'origin_count':q**3,'source_cells':L**3,
        'retained_fraction':fraction,'discarded_fraction':1-fraction,
        'classical_origin_register_dimension':q**3,
        'origin_label_bits_sufficient':(q**3-1).bit_length(),
        'origin_distribution':'uniform on the integer cube; one origin per complete history',
        'origin_register_retained':True,'origin_H':0,
        'large_origin_register_allocated':False}


def folner_bound(level,displacement=('1','0','0')):
    geo=origin_geometry(level);shift=spatial.rounded_translation(level,displacement)
    q=geo['origin_side'];m=shift['lattice_steps']
    overlap=1
    for k in m:overlap*=max(0,q-abs(k))
    ratio=2*(1-F(overlap,q**3))
    simple=min(F(2),sum((2*min(F(1),F(abs(k),q)) for k in m),F()))
    if ratio>simple:raise ArithmeticError('cube symmetric-difference estimate failed')
    return {'requested_physical_displacement':shift['requested'],
        'integer_displacement':m,'rounded_physical_displacement':shift['actual'],
        'origin_overlap_count':overlap,
        'symmetric_difference_count':2*(q**3-overlap),
        'unit_norm_expectation_difference_upper':ratio,
        'coordinate_sum_upper':simple,
        'finite_source_translation_symmetry_assumed':False}


def collar_weight(level,full_box_weight):
    w=ground._q(full_box_weight,'original full-box weight')
    if w>1:raise ValueError('original event weight cannot exceed one')
    geo=origin_geometry(level);fraction=geo['retained_fraction']
    return max(F(),(w-geo['discarded_fraction'])/fraction)


def bulk_bound(member,support_radius='1',duration='1/1000',support_cells=1):
    """Uniform origin-conditioned ambient comparison, including arbitrary references."""
    contract.verify_member(member)
    geo=origin_geometry(member['source']['level'])
    R=ground._q(support_radius,'support radius');T=ground._q(duration,'source duration')
    owners=ground._i(support_cells,'seed owner count',1)
    a=geo['source_spacing'];L=geo['source_side']
    r=geo['collar_cells']-(R/a).__ceil__()-2
    inherited=locality.selected_source_bound(member)
    W=F(inherited['family']['weighted_rate_times_spacing_upper'])
    if W!=F(5232,7):raise ArithmeticError('actual Wilson locality coefficient changed')
    cross=F(1272,7)*L*L/a;factor=2*owners*T*cross
    margin=F(r,2)-2*W*T/a
    exact_zero=T==0 and r>=1
    bits=-1
    if not exact_zero and r>=1 and margin>=0 and factor>0:
        bits=max(-1,margin.__floor__()-locality.cone.ceil_log2(factor))
    return {'support_radius':R,'source_time_extent':T,'seed_owner_count':owners,
        'all_origins_support_distance_lower':r,'all_origins_bulk_admitted':r>=1,
        'crossing_interaction_norm_upper':cross,'weighted_rate_times_spacing_upper':W,
        'Duhamel_prefactor':factor,'rational_exponential_margin':margin,
        'operator_error':{'exact_zero':exact_zero,'base':2,
            'negative_exponent':None if exact_zero else bits},
        'one_binary_effect_complete_half_trace_bound':'min(1,sqrt(2*operator_error))',
        'fixed_history_bound':'sum of the per-stage complete bounds, capped at one',
        'origin_conditioned_quantum_reference_retained':True,
        'origin_count_multiplier_in_error':False,
        'finite_physical_norm_embedding_assumed':False}


def same_origin_control():
    """The arithmetic observable mean and independent origins change histories."""
    N0=s.diag(0,1,0,1);N1=s.diag(0,0,1,1);psi=s.Matrix([0,1,0,0])
    choices=[N0,N1];same={};fresh={}
    for r,t in product((0,1),repeat=2):
        same[(r,t)]=s.Rational(1,2)*sum(
            ((A if t else s.eye(4)-A)*(A if r else s.eye(4)-A)*psi).norm()**2
            for A in choices)
        fresh[(r,t)]=s.Rational(1,4)*sum(
            ((B if t else s.eye(4)-B)*(A if r else s.eye(4)-A)*psi).norm()**2
            for A,B in product(choices,repeat=2))
    avg=(N0+N1)/2
    hom=s.diag(N0,N1)
    return {'same_origin':same,'fresh_origins':fresh,
        'homomorphic_second_moment':sum((psi.H*A*A*psi)[0] for A in choices)/2,
        'arithmetic_mean_second_moment':(psi.H*avg*avg*psi)[0],
        'block_product_defect':hom*hom-s.diag(N0*N0,N1*N1),
        'control_is_primary_source':False}


def complete_reference_control():
    """Two measured histories retain the input origin and an entangled reference."""
    identity=s.eye(4);N0=s.diag(0,1,0,1);N1=s.diag(0,0,1,1)
    hadamard=s.Matrix([[1,1],[1,-1]])/s.sqrt(2)
    U=s.kronecker_product(hadamard,s.eye(2))
    psi=s.zeros(16,1)
    for i,p in enumerate((s.Rational(1,2),s.Rational(1,4),s.Rational(1,8),s.Rational(1,8))):
        psi[4*i+i]=s.sqrt(p)
    reference=s.zeros(4);total=s.S(0);origins=[];rows=[]
    for x,A in enumerate((N0,N1)):
        ks=(A,identity-A);effects=s.zeros(4);origin_trace=s.S(0)
        for r,t in product(range(2),repeat=2):
            K=ks[t]*U*ks[r];effects+=K.H*K
            v=s.kronecker_product(K,s.eye(4))*psi/s.sqrt(2)
            density=v*v.H;prob=s.simplify((v.H*v)[0]);total+=prob;origin_trace+=prob
            for i,j in product(range(4),repeat=2):
                reference[i,j]+=sum(density[4*k+i,4*k+j] for k in range(4))
            rows.append({'origin':x,'first':r,'second':t,'probability':str(prob)})
        if effects!=identity:raise ArithmeticError('the complete history is not normalized')
        origins.append(origin_trace)
    expected=s.diag(s.Rational(1,2),s.Rational(1,4),s.Rational(1,8),s.Rational(1,8))
    return {'total_probability':s.simplify(total),'origin_probabilities':origins,
        'reference_marginal_difference':s.simplify(reference-expected),
        'records':rows,'reference_entangled_with_input':True,
        'two_source_times_use_same_controlled_origin':True,
        'control_is_primary_source':False}


@lru_cache(maxsize=32)
def _collective(level):
    member=contract.member(level)
    receipt=response.bind_member(member)
    response.verify_bound(member,receipt)
    return receipt


def selected_bound(member):
    contract.verify_member(member)
    n=member['source']['level'];geo=origin_geometry(n)
    if geo['source_side']!=member['source']['side'] or geo['source_spacing']!=F(member['source']['spacing']):
        raise ArithmeticError('origin geometry does not match the selected source')
    owning=_collective(n);data=owning['result']
    if (data['native_collective_block_dimension'],data['native_collective_Q6_charge'],
        data['collective_normalization_squared'])!=(6,1,42):
        raise ArithmeticError('collective physical seed differs')
    band=data['collective_original_transfer_band']
    protocol={'name':'retained_uniform_origin_complete_history.v1',
        'origin_geometry':contract.serial(geo),'member_id':member['member_id'],
        'one_origin_across_entire_history':True,
        'local_success':'C_x','local_failure':'sqrt(I-C_x^dagger*C_x)',
        'complete_one_stage_record_labels':'(origin x, success or failure)',
        'complete_one_stage_record_count':2*geo['origin_count'],
        'input':'uniform classical origin tensor unchanged canonical source/reference input',
        'clock':'I_origin tensor H_source tensor I_reference',
        'exact_apparatus_service_declared':True,'large_service_executed':False,
        'same_success_measure_as_uniform_address_choice':True,
        'same_complete_failure_law_as_aggregated_failure':False}
    return contract.serial({'member_id':member['member_id'],**member['ids'],
        'observation_protocol':protocol,'observation_protocol_id':contract.digest(protocol),
        'origin_geometry':geo,'fixed_physical_translation_example':folner_bound(n),
        'same_source_bulk_bound':bulk_bound(member),
        'owning_collective_response_sha256':owning['result_sha256'],
        'collective_original_transfer_band':band,
        'full_box_original_weight_lower':data['collective_original_transfer_weight_lower'],
        'origin_original_weight_lower':collar_weight(n,F(data['collective_original_transfer_weight_lower'])),
        'full_box_finite_record_weight_lower':data['finite_record_simple_weight_lower'],
        'origin_finite_record_weight_lower':collar_weight(n,F(data['finite_record_simple_weight_lower'])),
        'finite_record_transfer_band':data['finite_record_transfer_band'],
        'fixed_regular_pointer_histories_have_common_subsequence_limits':True,
        'relative_observation_limiting_state_translation_invariant':True,
        'relative_observation_canonical_H_nonnegative':True,
        'relative_observation_canonical_H_and_P_strongly_commute':True,
        'relative_observation_cyclic_vector_H_and_P_zero':True,
        'fixed_origin_selected_state_translation_invariant_proved':False,
        'nonzero_momentum_response_proved':False,'nonzero_regular_field_proved':False,
        'fixed_energy_survival_proved':False,'dispersion_or_mass_shell_proved':False,
        'large_origin_reference_joint_state_limit_claimed':False,
        'effective_subsequence_or_unique_limit_proved':False,
        'source_input_preparation_changed':False,'source_reference_clock_changed':False,
        'observer_observation_domain_explicitly_extended':True,
        'finite_Fourier_label_identified_as_limiting_momentum':False,
        'physical_mass_table_GeV':[]})


def bind_member(member):
    result=selected_bound(member)
    return {'schema':'oph.primary_origin_momentum.v1','result':result,
            'result_sha256':contract.digest(result)}


def verify_bound(member,receipt):
    if receipt!=bind_member(member):raise ValueError('origin observation result or primary source differs')


def certificate():
    paths=provenance.definition_inputs((*SOURCE_INPUTS,*spatial.SOURCE_INPUTS,*response.SOURCE_INPUTS))
    pins=provenance.definition_hashes(paths)
    common=same_origin_control();complete=complete_reference_control()
    if common['block_product_defect']!=s.zeros(8) or complete['reference_marginal_difference']!=s.zeros(4):
        raise ArithmeticError('exact retained-origin quantum control failed')
    result={'schema':'oph.primary_origin_momentum_certificate.v1',
        'selected_member_bound':bind_member(contract.member(0)),
        'refinement_members':[bind_member(contract.member(n)) for n in (6,13,16)],
        'same_origin_control':{'same_origin':{str(k):str(v) for k,v in common['same_origin'].items()},
            'fresh_origins':{str(k):str(v) for k,v in common['fresh_origins'].items()},
            'homomorphic_second_moment':str(common['homomorphic_second_moment']),
            'arithmetic_mean_second_moment':str(common['arithmetic_mean_second_moment'])},
        'complete_reference_control':{'total_probability':str(complete['total_probability']),
            'origin_probabilities':[str(x) for x in complete['origin_probabilities']],
            'reference_marginal_preserved':True,'records':complete['records']},
        'source_sha256':pins}
    if provenance.definition_hashes(paths)!=pins:raise ValueError('origin inputs changed during replay')
    return result
