"""Polynomial detailed balance and original conserved-sector energy records."""
from fractions import Fraction as F
from functools import lru_cache
from math import prod

import sympy as s

from model.particle_common_contract import contract, provenance
from model.particle_primary_coherent_odd_response import model as coherent
from model.particle_primary_composite_spin import model as composite
from model.particle_primary_gravity_premises import model as modular
from model.particle_spatial_modular_response import regional

SOURCE_INPUTS=tuple(dict.fromkeys((
 'model/particle_primary_original_branch/__init__.py',
 'model/particle_primary_original_branch/__main__.py',
 'model/particle_primary_original_branch/model.py',
 'model/particle_primary_original_branch/test_model.py',
 'model/particle_primary_original_branch/PROOF.md',
 'model/particle_primary_original_branch/README.md',
 *coherent.SOURCE_INPUTS,*composite.SOURCE_INPUTS,*modular.SOURCE_INPUTS,*regional.SOURCE_INPUTS,
)))
DELTA=coherent.spin.DELTA


def dyadic_exponential_upper(exponent,retained_bits=16):
    """e^(-x)<=2^(-floor(x)); use only exact rational comparisons."""
    x=F(exponent)
    if x<0 or type(retained_bits) is not int or retained_bits<0:raise ValueError('nonnegative exponent and integer precision required')
    maximum=x.numerator//x.denominator;bits=min(maximum,retained_bits)
    return {'certified_negative_power_of_two':maximum,'arithmetic_negative_power_of_two':bits,
            'arithmetic_upper':F(1,2**bits),'floating_exponential_evaluated':False}


def paired_window_lower(signed_union,negative_to_positive_ratio):
    union=F(signed_union);ratio=F(negative_to_positive_ratio)
    if union<0 or ratio<0:raise ValueError('positive signed window and ratio bounds required')
    return union/(1+ratio)


@lru_cache(maxsize=1)
def inversion_control():
    gamma=composite.sectors.branch.response.density.gamma_data();beta=gamma[3]
    for g in gamma[:3]:
        hop=(-beta+beta*g)/2
        if beta*hop*beta!=hop.H:raise ArithmeticError('physical beta inversion does not reverse a source hop')
        if (-s.eye(4))*hop*(-s.eye(4))==hop.H:raise ArithmeticError('central spin parity spuriously implements inversion')
    for q in composite.group.elements():
        S=composite.channels.spin_lift(q)
        if beta*S!=S*beta:raise ArithmeticError('separate inversion does not commute with proper source lifts')
    if beta**2!=s.eye(4):raise ArithmeticError('source inversion is not an involution')
    return {'proper_binary_cube_elements':48,'separate_inversion_order':2,
        'source_inversion':'cell x to 1-x, oriented link reversal, Dirac beta lift',
        'full_rotation_inversion_group_order':96,'beta_reverses_actual_Wilson_hop':True,
        'beta_commutes_with_all_proper_source_lifts':True,
        'central_two_pi_CAR_parity_is_spatial_inversion':False,
        'word_projection_and_physical_Gauss_invariance_reused':True}


@lru_cache(maxsize=1)
def fusion_control():
    group=composite.group;tables={};central={}
    minus=group.negative(group.ONE)
    for name,dimension in zip(group.IRREPS,group.DIMENSIONS):
        central[name]=(group.characters(minus)[name]/dimension).integer()
    for probe in ('G1','H'):
        table={}
        for initial in group.IRREPS:
            row=group.decompose(lambda q,initial=initial,probe=probe:group.characters(q)[initial]*group.characters(q)[probe])
            size=sum(group.DIMENSIONS[group.IRREPS.index(label)]*n for label,n in row.items())
            if size!=group.DIMENSIONS[group.IRREPS.index(initial)]*group.DIMENSIONS[group.IRREPS.index(probe)]:
                raise ArithmeticError('literal finite tensor-product dimension differs')
            if any(central[target]!=-central[initial] for target in row):
                raise ArithmeticError('odd probe does not flip the actual central spin lift')
            for parity in ('g','u'):
                table[initial+parity]={target+parity:n for target,n in row.items()}
        count=sum(len(row) for row in table.values())
        if count!=(28 if probe=='G1' else 44):raise ArithmeticError('complete cube transition menu changed')
        tables[probe]={'uniform_operator_irrep':probe+'g','allowed_transition_multiplicities':table,
            'allowed_initial_final_pairs':count,'central_spin_lift_reverses':True,
            'inversion_label_preserved_by_positive_uniform_probe':True}
    return {'actual_source_group_order':96,'separate_source_inversion':inversion_control(),'complete_total_cube_irrep_labels':[name+p for p in ('g','u') for name in group.IRREPS],
        'proper_character_orthogonality_reused':True,'probe_tables':tables,
        'projectors_are_source_total_irrep_projectors':True,'continuous_Lorentz_spin_projectors_claimed':False}


@lru_cache(maxsize=1)
def source_multiplicity_vectors():
    """Sparse actual480-mode matter vectors; every link is the constant function."""
    V=8;omega=sum(composite.OMEGA<<(60*x) for x in range(V))
    local=coherent.spin.local_controls();active=local['active_source_modes'];e=dict(composite.sectors.native_blocks())[6][0]
    vectors=[]
    for label,family in local['families'].items():
        for component,polynomial in enumerate(family['polynomials']):
            three={};seven={}
            for x in range(V):
                for word,coefficient in polynomial.items():
                    ops=tuple((60*x+active[i],True) for i in word)
                    target,phase=composite.apply_operators(omega,ops)
                    if not phase:raise ArithmeticError('uniform actual three-defect component vanished')
                    three[target]=three.get(target,0)+coefficient*phase/s.sqrt(8)
                    for y in range(V):
                        if x==y:continue
                        neutral=((60*y+e,True),(60*y+15+e,True),(60*y+45+e,False),(60*y+30+e,False))
                        target,phase=composite.apply_operators(omega,ops+neutral)
                        if not phase:raise ArithmeticError('uniform actual seven-defect component vanished')
                        seven[target]=seven.get(target,0)+coefficient*phase/s.sqrt(56)
            for state,defects in ((three,3),(seven,7)):
                state={mask:s.simplify(value) for mask,value in state.items() if value}
                if any((mask^omega).bit_count()!=defects for mask in state):raise ArithmeticError('actual full native defect number differs')
                if s.simplify(sum(s.conjugate(v)*v for v in state.values()))!=1:raise ArithmeticError('uniform physical multiplicity vector is not normalized')
            vectors.append((label,component,three,seven))
    for i,(_,_,left3,left7) in enumerate(vectors):
        for j,(_,_,right3,right7) in enumerate(vectors):
            for left,right,same in ((left3,right3,True),(left7,right7,True),(left3,right7,False)):
                inner=s.simplify(sum(s.conjugate(value)*right.get(mask,0) for mask,value in left.items()))
                if inner!=int(same and i==j):raise ArithmeticError('actual source multiplicity vectors are not orthonormal')
    return tuple(vectors)


@lru_cache(maxsize=1)
def physical_multiplicity_witness():
    """Actual local neutral double particle-hole scalar and uniform witnesses."""
    blocks=dict(composite.sectors.native_blocks());e=blocks[6][0]
    modes=(e,15+e,30+e,45+e)
    ops=((modes[0],True),(modes[1],True),(modes[3],False),(modes[2],False))
    target,phase=composite.apply_operators(composite.OMEGA,ops)
    if abs(phase)!=1 or (target^composite.OMEGA).bit_count()!=4:
        raise ArithmeticError('neutral double particle-hole state differs')
    for _,indices in composite.sectors.native_blocks():
        if sum((1 if creation else -1) for i,creation in ops if i%15 in indices):
            raise ArithmeticError('scalar witness changes a conserved native occupation')
    for q in composite.group.elements():
        S=composite.clean(composite.beta_chart().H*composite.channels.spin_lift(q)*composite.beta_chart())
        if S[:2,2:]!=s.zeros(2) or S[2:,:2]!=s.zeros(2):raise ArithmeticError('source rotation changes beta rows')
        if s.simplify(S[:2,:2].det()*s.conjugate(S[2:,2:].det()))!=1:
            raise ArithmeticError('neutral four-defect operator is not a scalar')
        if s.simplify(S[2:,2:].det())!=1:raise ArithmeticError('filled-beta source anchor is not a proper scalar')
    vectors=source_multiplicity_vectors()
    anchor_counts=tuple(16*len(indices) for _,indices in composite.sectors.native_blocks())
    native_masks=[sum(1<<(60*x+15*spin+i) for x in range(8) for spin in range(4) for i in indices) for _,indices in composite.sectors.native_blocks()]
    for _,_,three,seven in vectors:
        for state in (three,seven):
            for mask in state:
                if tuple((mask&block).bit_count()-n for block,n in zip(native_masks,anchor_counts))!=DELTA:
                    raise ArithmeticError('actual full source witness changed its native occupation sector')
    source_states=[]
    for x in range(8):
        for y in range(8):
            if x==y:continue
            # The q6 flipped site is unique; the weak-pair site is distinct.
            source_states.append((x,y))
    if len(set(source_states))!=56:raise ArithmeticError('uniform seven-defect witness normalization differs')
    return {'actual_q6_modes':modes,'neutral_operator_word':ops,'local_anchor_phase':phase,
        'neutral_native_increment':(0,0,0,0,0),'neutral_added_defects':4,
        'neutral_operator_cube_irrep':'A1g','neutral_operator_CAR_parity':0,
        'source_anchor_total_cube_irrep':'A1g','source_anchor_native_occupations':anchor_counts,
        'source_anchor_filled_modes':240,'uniform_three_defect_sites':8,'uniform_seven_defect_site_pairs':56,
        'actual_orthonormal_source_vectors':12,
        'source_matter_sparse_term_counts':[(label,component,len(three),len(seven)) for label,component,three,seven in vectors],
        'three_defect_source_energy_mean':F(-711),'seven_defect_source_energy_mean':F(-699),
        'same_current_and_G1g_Hg_sector_has_two_orthogonal_multiplicity_vectors':True,
        'each_actual_reducing_sector_spectral_diameter_lower':F(12),
        'distinct_means_are_claimed_eigenvalues':False,'isolated_lowest_branch_is_ruled_out':False}


@lru_cache(maxsize=1)
def original_transfer_spread_control():
    """Exact first residual action on a physical anchor, used only in a positive norm."""
    gamma=composite.sectors.branch.response.density.gamma_data();chart=composite.beta_chart()
    for g in gamma[:3]:
        hop=composite.clean(chart.H*((-gamma[3]+gamma[3]*g)/2)*chart)
        if hop[:2,:2]!=-s.eye(2)/2 or hop.H[:2,:2]!=-s.eye(2)/2:
            raise ArithmeticError('actual surviving positive-beta hopping block differs')
    residuals=[]
    for label,family in coherent.spin.local_controls()['families'].items():
        for component,polynomial in enumerate(family['polynomials']):
            state={}
            for word,coefficient in polynomial.items():
                for position in range(3):
                    changed=list(word);changed[position]+=6
                    target,phase=composite.apply_operators(0,tuple((i,True) for i in changed))
                    if not phase:raise ArithmeticError('first residual creation word vanished')
                    state[target]=state.get(target,0)+coefficient*phase
            if any((mask&63).bit_count()!=2 or (mask>>6).bit_count()!=1 for mask in state):
                raise ArithmeticError('residual left its actual two-plus-one occupation profile')
            residuals.append((label,component,state))
    gram=s.Matrix([[s.simplify(sum(s.conjugate(v)*right.get(mask,0) for mask,v in left.items()))
                    for _,_,right in residuals] for _,_,left in residuals])
    if gram!=3*s.eye(6):raise ArithmeticError('canonical vacuum residual Gram is not3I6')
    graph=coherent.cube_stars();profiles=set()
    for x,y in graph['edges']:
        for origin,target in ((x,y),(y,x)):
            profile=tuple(2 if z==origin else 1 if z==target else 0 for z in graph['cells'])
            if profile in profiles:raise ArithmeticError('residual origins have coincident matter profiles')
            profiles.add(profile)
    return {'raw_normalized_component_anchor_success_norm_squared':F(1),
        'canonical_six_component_residual_Gram_diagonal':3,
        'surviving_positive_beta_hopping_image_Gram':F(3,4),
        'actual_each_component_residual_anchor_norm_squared':F(9,4),
        'actual_oriented_cube_residual_profiles':len(profiles),
        'source_energy_commutator_witness_center':F(9),
        'actual_original_mean_transfer_computed':False,
        'all_first_action_link_words_admitted_by_degree_one_cutoff':True,
        'Gauss_projection_preserves_the_first_residual_action':True,
        'residual_anchor_orthogonal_to_created_anchor':True,
        'same_result_for_every_normalized_spatial_coefficient_vector':True,
        'raw_component_centered_moment_lower':'c*((9-center)^2+9/4)',
        'normalized_original_transfer_variance_lower_factored':{'multiplier_of_physical_density_floor':F(9,4)},
        'every_component_original_transfer_measure_is_non_single_line':True,
        'density_inequality_applied_to_source_marginal_positive_norm_only':True,
        'actual_line_locations_or_separation_computed':False,
        'nonzero_variance_excludes_an_isolated_pole':False,
        'current_clock_error_resolves_this_variance_floor':False}


@lru_cache(maxsize=1)
def _source_zero():
    member=contract.member(0)
    response=coherent.bind_member(member);coherent.verify_bound(member,response)
    sectors=composite.sectors.bind_member(member);composite.sectors.verify_bound(member,sectors)
    floor=regional.common_contract_payload(member,1)
    return response,sectors,floor


def bind_member(member):
    contract.verify_member(member)
    if member['source']['level']!=0:raise ValueError('original resolved records are certified only at member zero')
    inherited,sectors,floor_payload=_source_zero();r=inherited['result'];sector=sectors['result']
    slope=modular.selected_modular_slope(member)
    ratio=dyadic_exponential_upper(slope['modular_slope_lower']/2)
    if ratio['certified_negative_power_of_two']!=11135:raise ArithmeticError('selected source modular slope changed')
    clock_loss=F(r['complete_two_clock_probability_error_upper'])
    total=paired_window_lower(F(r['centered_signed_union_lower']),ratio['arithmetic_upper'])-clock_loss
    if total<=F(1,3):raise ArithmeticError('positive original source window lost its detailed-balance improvement')
    channels=[]
    for row,simple in zip(r['spin_resolved_channels'],(F(5,48),F(11,48))):
        label=row['probe_irrep'];signed=F(row['signed_original_success_lower'])-F(row['signed_centered_second_moment_upper'])/F(17,2)**2
        weight=paired_window_lower(signed,ratio['arithmetic_upper'])-clock_loss
        count=fusion_control()['probe_tables'][label]['allowed_initial_final_pairs']
        if weight<=simple:raise ArithmeticError('spin-resolved original detailed-balance window failed')
        channels.append({'probe_irrep':label,'uniform_operator_total_cube_irrep':label+'g',
            'reported_original_selected_phase_weight_lower':weight,'simple_selected_phase_weight_lower':simple,
            'reported_original_complete_Walsh_row_weight_lower':weight/8,'simple_complete_Walsh_row_weight_lower':simple/8,
            'uniform_Walsh_some_allowed_total_cube_transition_weight_lower':weight/(8*count),
            'simple_uniform_Walsh_some_allowed_total_cube_transition_weight_lower':simple/(8*count),
            'particular_robustly_populated_initial_final_pair_identified':False})
    density_floor=floor_payload['physical_density_floor_factored']
    if density_floor['base']!=2 or density_floor['operator']!='P_G':raise ArithmeticError('physical density floor changed')
    actual_floor=F(sector['sector_results'][0]['signed_sectors'][0]['physical_sector_energy_lower'])
    lower=F(-18055527,25000)
    if lower>actual_floor:raise ArithmeticError('short sector floor is not outward')
    mean=F(-711);plan=composite.sectors.spectroscopy.clock_plan(member,'1/64','1/1024','1/1024')
    eta=F(plan['energy_resolution']);loss=F(plan['clock_plan']['complete_clock_tail_upper'])+F(plan['one_complete_clock_error_upper'])
    witness_weight=F(1,2)-loss
    if witness_weight<=F(255,512):raise ArithmeticError('actual complete final-energy witness lost its bound')
    source_window=(lower-eta,2*mean-lower+eta)
    excitation=(source_window[0]-F(sector['global_ground_upper']),source_window[1]-F(sector['global_ground_lower']))
    if excitation[0]<=0:raise ArithmeticError('actual quantum-label window is not above the full ground')
    anchor=physical_multiplicity_witness();exact_branches=[]
    for label,dimension in (('G1',2),('H',4)):
        exact_branches.append({'initial_native_occupations':anchor['source_anchor_native_occupations'],
            'initial_total_cube_irrep':'A1g','final_native_occupations':tuple(n+d for n,d in zip(anchor['source_anchor_native_occupations'],DELTA)),
            'final_total_cube_irrep':label+'g','source_C_conjugate_final_cube_irrep':label+'u',
            'source_C_conjugate_native_increment':tuple(-d for d in DELTA),
            'complete_uniform_Walsh_original_weight_lower_factored':{'physical_density_floor':density_floor,'multiplier':F(dimension,48)*witness_weight},
            'simple_original_weight_strict_lower_factored':{'multiplier_of_physical_density_floor':F(dimension*255,48*512)},
            'original_source_energy_window':source_window,'inferred_excitation_interval':excitation,
            'initial_energy_record_applied':False,'energy_transfer_law_claimed':False,
            'source_input_replaced_by_witness':False,'prepared_reference_density_domination_claimed':False,
            'isolated_pole_or_species_identified':False})
    result=contract.serial({'member_id':member['member_id'],**member['ids'],
        'coherent_response_receipt_sha256':contract.digest(inherited),'sector_threshold_receipt_sha256':contract.digest(sectors),
        'source_modular_slope':slope,'paired_adjoint_ratio_upper':ratio,
        'actual_polynomial_state_replaced_by_Gibbs':False,'negative_probe_is_adjoint_not_source_C_conjugated_creator':True,
        'reported_positive_creation_interval':r['finite_reported_creation_interval'],
        'reported_original_total_selected_phase_weight_lower':total,'simple_total_selected_phase_weight_lower':F(1,3),
        'reported_original_total_complete_Walsh_row_weight_lower':total/8,'simple_total_complete_Walsh_row_weight_lower':F(1,24),
        'complete_two_clock_plan':r['complete_clock_plan'],'complete_two_clock_event_error_upper':clock_loss,
        'spin_resolved_records':channels,'complete_total_cube_current_refinement':fusion_control(),
        'original_energy_law_is_marginal_of_refined_complete_records':True,
        'initial_label_measurement_preserves_source_marginal':True,
        'initial_label_measurement_preserves_source_reference_joint_state':False,
        'explicit_original_final_energy_branches':exact_branches,'final_only_clock_plan':plan,
        'physical_density_floor_payload_sha256':contract.digest(floor_payload),
        'same_label_source_multiplicity_obstruction':anchor,
        'original_transfer_spread':{**original_transfer_spread_control(),'physical_density_floor_factored':density_floor},
        'all_failure_reference_native_current_cube_and_energy_records_retained':True,
        'conserved_spatial_momentum_or_relativistic_dispersion_proved':False,
        'isolated_original_particle_branch_proved':False,'physical_mass_table_GeV':[]})
    return {'schema':'oph.primary_original_branch.v1','result':result,'result_sha256':contract.digest(result)}


def verify_bound(member,value):
    if value!=bind_member(member):raise ValueError('original branch consequence or source binding changed')


def certificate():
    paths=provenance.definition_inputs(SOURCE_INPUTS);pins=provenance.definition_hashes(paths)
    value={'schema':'oph.primary_original_branch_certificate.v1','members':[bind_member(contract.member(0))],
        'source_sha256':pins}
    if provenance.definition_hashes(paths)!=pins:raise ValueError('source changed during original branch replay')
    return contract.serial(value)
