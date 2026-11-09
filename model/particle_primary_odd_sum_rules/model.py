"""Original odd spectral sum rules and nonzero graded source propagation."""
from fractions import Fraction as F
from functools import lru_cache

import sympy as s

from model.particle_common_contract import contract, provenance
from model.particle_primary_coherent_odd_response import model as coherent
from model.particle_primary_composite_spin import model as composite
from model.particle_primary_composite_transport import model as transport
from model.particle_spatial_modular_response import regional
from model.particle_primary_original_branch import model as branch

SOURCE_INPUTS=tuple(dict.fromkeys((
 'model/particle_primary_odd_sum_rules/__init__.py',
 'model/particle_primary_odd_sum_rules/__main__.py',
 'model/particle_primary_odd_sum_rules/model.py',
 'model/particle_primary_odd_sum_rules/test_model.py',
 'model/particle_primary_odd_sum_rules/PROOF.md',
 'model/particle_primary_odd_sum_rules/README.md',
 'model/particle_primary_coherent_odd_response/KERNEL_PROOF.md',
 'model/particle_physical_boundary_state/PROOF.md',
 'model/particle_primary_charged_string/PROOF.md',
 'Lean/Screen/SMStructureAdequacySurface.lean',
 'model/standard_model_source_map/KINETIC_ATTACHMENT.md',
 *coherent.SOURCE_INPUTS,*composite.SOURCE_INPUTS,*transport.SOURCE_INPUTS,*regional.SOURCE_INPUTS,*branch.SOURCE_INPUTS,
)))


def _add(out,state,coefficient=1):
    for mask,value in state.items():out[mask]=s.expand(out.get(mask,0)+coefficient*value)


def _apply_word(state,word):
    out={}
    for mask,value in state.items():
        target,phase=composite.apply_operators(mask,word)
        if phase:out[target]=out.get(target,0)+value*phase
    return out


def _apply_field(state,polynomial,site,adjoint=False):
    active=coherent.spin.local_controls()['active_source_modes'];out={}
    for word,value in polynomial.items():
        order=tuple(reversed(word)) if adjoint else word
        ops=tuple((60*site+active[i],not adjoint) for i in order)
        _add(out,_apply_word(state,ops),s.conjugate(value) if adjoint else value)
    return {mask:s.simplify(value) for mask,value in out.items() if value}


def _apply_neighbor_derivative(state,polynomial,charge,left,right):
    """Coefficient of the actual native R_charge(U)[left,right]."""
    active=coherent.spin.local_controls()['active_source_modes'];blocks=dict(composite.sectors.native_blocks());out={}
    for word,value in polynomial.items():
        for position,i in enumerate(word):
            spin,native=divmod(active[i],15)
            if native!=blocks[charge][right]:continue
            ops=[(60+active[k],True) for k in word]
            ops[position]=(15*spin+blocks[charge][left],True)
            _add(out,_apply_word(state,tuple(ops)),-value/2)
    return {mask:s.simplify(value) for mask,value in out.items() if value}


def _gram_norm(state):
    return s.simplify(sum(s.conjugate(value)*value for value in state.values()))


@lru_cache(maxsize=1)
def neighbor_graded_control():
    """Full120-mode CAR and actual native Haar norm of the physical first action."""
    local=coherent.spin.local_controls();blocks=dict(composite.sectors.native_blocks())
    gamma=composite.sectors.branch.response.density.gamma_data();chart=composite.beta_chart()
    for g in gamma[:3]:
        hop=composite.clean(chart.H*((-gamma[3]+gamma[3]*g)/2)*chart)
        if hop[:2,:2]!=-s.eye(2)/2 or hop.H[:2,:2]!=-s.eye(2)/2:
            raise ArithmeticError('actual positive-beta Wilson block differs')
    omega=composite.OMEGA|(composite.OMEGA<<60);rows=[]
    for label,family in local['families'].items():
        for component,polynomial in enumerate(family['polynomials']):
            witness=_apply_field({omega:s.S.One},polynomial,0)
            if _gram_norm(witness)!=1:raise ArithmeticError('physical trion norm differs')
            coefficients=[];norm=s.S.Zero
            for charge in (-3,6):
                dimension=len(blocks[charge])
                for left in range(dimension):
                    for right in range(dimension):
                        output={}
                        _add(output,_apply_field(_apply_neighbor_derivative(witness,polynomial,charge,left,right),polynomial,0,True))
                        _add(output,_apply_neighbor_derivative(_apply_field(witness,polynomial,0,True),polynomial,charge,left,right))
                        output={mask:s.simplify(value) for mask,value in output.items() if value}
                        for mask in output:
                            if (mask&((1<<30)-1)).bit_count()!=1 or ((mask>>60)&((1<<30)-1)).bit_count()!=2:
                                raise ArithmeticError('graded first action is not a two-particle transfer')
                        square=_gram_norm(output);norm+=square/dimension
                        coefficients.append({'Q6':charge,'left_index':left,'right_index':right,
                            'native_dimension':dimension,'matter_norm_squared':F(square),
                            'nonzero_matter_terms':len(output)})
            norm=F(norm);expected=F(3,8) if label=='G1' else (F(3,4) if component in (0,3) else F(1,4))
            if norm!=expected:raise ArithmeticError('actual Haar/CAR graded propagation norm differs')
            rows.append({'probe_irrep':label,'component':component,'physical_witness_norm_squared':1,
                'first_graded_response_witness_norm_squared':norm,'native_Haar_coefficients':coefficients})
    # Exact exterior chart conversion is part of the inherited source link action.
    permutation=transport.native_chart();transport.nontrivial_link()
    return {'source_matter_modes_in_explicit_control':120,'actual_source_axes_checked':3,
        'source_native_chart_permutation':permutation,'minimum_physical_witness_norm_squared':F(1,4),
        'all_normalized_physical_trion_witnesses':rows,
        'first_action_link_words_admitted_by_actual_degree_one_cutoff':True,
        'first_graded_response_vanishes_for_distinct_nonneighbors':True,
        'physical_witness_substitutes_for_original_density':False,
        'original_signed_first_moment_is_the_positive_norm':False}


@lru_cache(maxsize=1)
def sum_rule_control():
    local=coherent.local_anticommutators();geometry=coherent.kernel_reduction.orbit_control()
    diagonal,neighbor=s.symbols('m0 m1',real=True)
    moments=coherent.kernel_reduction.walsh_eigenresponses((diagonal,neighbor,0,0))
    if moments!=tuple(diagonal+(3-2*k)*neighbor for k in range(4)):
        raise ArithmeticError('nearest-neighbor fermionic sum rule changed')
    return {'actual_local_signed_CAR_effects':local,'source_four_shell_orbits':geometry,
        'positive_measure':'creation transfer plus reflected, transposed adjoint-annihilation transfer',
        'zeroth_spatial_kernel_is_diagonal':True,
        'complete_Walsh_signed_success_effect':'sum_(x,j) {F_xj*,F_xj}/384',
        'all_eight_signed_success_effects_identical_as_operators':True,
        'complete_energy_clocks_preserve_this_original_success_identity':True,
        'first_summed_shells':['m0','m1','0','0'],
        'first_Walsh_moments':['m0+3*m1','m0+m1','m0-m1','m0-3*m1'],
        'first_complete_menu_moment_divisor':48,
        'm1_is_actual_Wilson_dressed_pair_transfer_expectation':True,
        'm1_original_value_or_sign_computed':False,
        'creation_only_success_effects_must_be_spatially_identical':False,
        'retarded_mean_response_is_inferred_from_nonzero_operator':False}


def finite_time_plan(norm,density_exponent):
    R=F(norm);N=int(density_exponent)
    if R<1 or N<0 or N!=density_exponent:raise ValueError('source norm at least one and exact nonnegative exponent required')
    # c=2^-N is never expanded. Every power and rational prefactor is explicit.
    t=F(1,64)/R**2
    nearest=F(1,64)*t*t
    nonneighbor=9*R**4*t**4
    if nonneighbor/nearest!=F(9,64):raise ArithmeticError('strict finite contrast factor differs')
    # With c<=1, Taylor and one-step Cayley errors leave amplitude at least t*c/8.
    relative_remainder=2*R**2*t+R**3*t*t/6
    if relative_remainder>=F(1,8):raise ArithmeticError('positive first-action amplitude lost')
    return {'source_norm_upper':R,'density_floor_negative_exponent':N,
        'time':{'base':2,'negative_exponent':N,'rational_prefactor':t},
        'source_time_formula':'t=c/(64R^2), c=2^(-N)',
        'clock':'one exact Cayley step of the same H, and its inverse',
        'clock_formula':'U_t=(I+i*t*H/2)^(-1)*(I-i*t*H/2)',
        'one_clock_operator_error_upper':'R^3*t^3/12',
        'graded_success_operator':'(F_x* U_t* F_y U_t + U_t* F_y U_t F_x*)/2',
        'complete_Kraus_operators':['K_xy(t)','sqrt(I-K_xy(t)* K_xy(t))'],
        'complete_outcome_count':2,
        'neighbor_original_click_probability_lower':{'base':2,'negative_exponent':3*N,'rational_prefactor':nearest},
        'nonneighbor_original_click_probability_upper':{'base':2,'negative_exponent':4*N,'rational_prefactor':nonneighbor},
        'neighbor_minus_nonneighbor_original_probability_lower':{'base':2,'negative_exponent':3*N,'rational_prefactor':nearest/2},
        'nonneighbor_to_neighbor_bound_ratio_upper':F(9,64),
        'time_and_weight_uniform_in_capacity':False,
        'exponential_size_integers_or_matrices_allocated':False,
        'all_failures_source_and_reference_retained':True,
        'full_original_state_histogram_executed':False}


@lru_cache(maxsize=1)
def gauge_charge_control():
    delta=coherent.spin.DELTA;blocks=composite.sectors.native_blocks();charge=sum(q*n for (q,_),n in zip(blocks,delta))
    if charge:raise ArithmeticError('the actual trion has nonzero global hypercharge')
    kernel=composite.local_multiplet(delta)
    if kernel['kernel_basis'].cols!=6:raise ArithmeticError('actual minimal local singlet kernel changed')
    if len(kernel['gauge_matrices'])!=12 or any(A*kernel['kernel_basis']!=s.zeros(A.rows,6) for A in kernel['gauge_matrices']):
        raise ArithmeticError('actual full gauge-kernel annihilation failed')
    # The native operator commutators are proved by the injective minimal-word map.
    return {'native_charge_order':[q for q,_ in blocks],'trion_native_increment':delta,
        'trion_total_Q6_increment':charge,'actual_local_gauge_kernel_dimension':6,
        'all_twelve_local_gauge_generators_annihilate_trion_kernel':True,
        'trion_is_local_color_and_weak_singlet':True,
        'trion_native_increment_is_nonzero_electric_charge':False,
        'source_C_makes_a_charged_trion_antiparticle':False,
        'same_source_outer_Gauss_boundary_external_charge':0,
        'known_physical_string_endpoint_Q6_increments':[-6,6],
        'known_physical_string_probe_CAR_parity':0,
        'single_gauge_covariant_creator_preserves_original_physical_carrier':False,
        'charged_regional_particle_requires_compensating_dressing_and_separate_spin_dispersion_join':True,
        'neutral_trion_is_identified_as_electron_muon_tau_or_quark':False,
        'gauge_singlet_alone_identifies_or_excludes_a_dressed_neutrino':False}


@lru_cache(maxsize=1)
def chirality_parity_control():
    gamma=composite.sectors.branch.response.density.gamma_data();beta=gamma[3]
    gamma5=composite.clean(gamma[0]*gamma[1]*gamma[2]*beta)
    if gamma5.H!=gamma5 or gamma5**2!=s.eye(4) or beta*gamma5*beta!=-gamma5:
        raise ArithmeticError('native Clifford chirality is not parity exchanged')
    plus=(s.eye(4)+gamma5)/2;minus=(s.eye(4)-gamma5)/2
    if plus.rank()!=2 or minus.rank()!=2 or plus*minus!=s.zeros(4) or beta*plus*beta!=minus:
        raise ArithmeticError('native chirality projectors differ')
    if beta*gamma5==gamma5*beta:raise ArithmeticError('onsite source chirality was spuriously conserved')
    for q in composite.group.elements():
        D=composite.channels.spin_lift(q)
        if D*gamma5!=gamma5*D:raise ArithmeticError('proper source spin lift changes native chirality')
    native=composite.sectors.native.native_currents()
    for G in native:
        for projector in (plus,minus):
            A=s.kronecker_product(projector,s.eye(15));B=s.kronecker_product(s.eye(4),G)
            if A*B!=B*A:raise ArithmeticError('native chirality changed an internal gauge current')
    return {'native_gamma5':[[str(v) for v in row] for row in gamma5.tolist()],
        'native_chiral_projectors':[[[str(v) for v in row] for row in A.tolist()] for A in (plus,minus)],
        'projector_Dirac_ranks':[2,2],'same_internal_representation_in_both_projectors':True,
        'same_native_gauge_currents_in_both_projectors':True,
        'parity_conjugates_internal_gauge_representation':False,
        'onsite_H_conserves_native_gamma5':False,
        'proper_spin_lifts_checked':48,'native_gauge_generators_checked':len(native),
        'actual_source_inversion':branch.inversion_control(),
        'native_gamma5_is_already_proved_physical_Lorentz_chirality':False}


def bind_chirality_member(member):
    contract.verify_member(member);source=member['source'];L=source['side'];B=source['interior_hopping_bonds']
    if B!=3*L*L*(L-1):raise ArithmeticError('actual open source interior graph differs')
    controls=chirality_parity_control()
    result=contract.serial({'member_id':member['member_id'],**member['ids'],
        'source_side':L,'source_link_degree':source['per_link_word_degree'],
        'physical_inversion':'x -> (L-1)-x; reverse interior links; beta on Dirac rows',
        'asymmetric_ambient_ghost_links_are_not_rotated_as_an_invariant_graph':True,
        'one_link_ghost_Gauss_constraints_reduce_every_dangling_link_to_constant':True,
        'interior_graph_inversion_preserves_every_word_cutoff_degree':True,
        'actual_source_polynomial_density_is_parity_invariant':True,
        'native_Clifford_controls':controls,
        'fixed_native_block_probe':'B_(s,x,y)=c_x* (Gamma_s tensor R_r(U_xy)) c_y P0_xy',
        'native_block_probe_norm_upper':'2*d_r',
        'complete_success_records_per_native_block':4*B,'complete_failure_records':1,
        'complete_success_Kraus':'B_(s,x,y)/(2*d_r*sqrt(4*B))',
        'complete_failure_Kraus':'sqrt(I-sum K* K)',
        'complete_original_clock_laws_equal_under_chirality_and_position_reflection':True,
        'source_initial_and_final_energy_bins_are_unchanged_by_parity':True,
        'native_representation_and_all_five_native_currents_unchanged_by_parity':True,
        'controlled_Wilson_dressing_and_controller_must_also_be_parity_conjugated':True,
        'retained_source_reference_output_states_asserted_equal_without_conjugation':False,
        'native_chirality_projectors_define_distinct_particle_species':False,
        'one_handed_physical_low_energy_sector_proved':False,
        'universal_no_go_for_parity_broken_or_reinterpreted_models':False,
        'extra_premise_for_one_handed_identification':['specified parity-breaking source state or selection',
            'different source action','proved parity or current dictionary that changes the attachment'],
        'physical_mass_table_GeV':[]})
    return {'schema':'oph.primary_native_chirality_pairing.v1','result':result,'result_sha256':contract.digest(result)}


def verify_chirality_bound(member,value):
    if value!=bind_chirality_member(member):raise ValueError('native chirality pairing or source binding changed')


@lru_cache(maxsize=1)
def _source():
    member=contract.member(0);value=coherent.bind_member(member);coherent.verify_bound(member,value)
    return value,regional.common_contract_payload(member,1)


def bind_member(member):
    contract.verify_member(member)
    if member['source']['level']!=0:raise ValueError('odd cube sum rules are certified only at member zero')
    inherited,floor=_source();density=floor['physical_density_floor_factored'];N=density['negative_exponent']
    if density['base']!=2 or density['operator']!='P_G':raise ArithmeticError('original physical floor differs')
    source=member['preparation_plan']['source'];R=max(abs(F(source['spectral_lower'])),abs(F(source['spectral_upper'])))
    lower=F(inherited['result']['original_total_signed_success_lower'])
    channels=[]
    for label,d in (('G1',2),('H',4)):
        channels.append({'probe_irrep':label,'dimension':d,
            'each_complete_Walsh_signed_success_lower':F(d,48)*lower,
            'each_complete_Walsh_signed_success_upper':F(d,48),
            'same_original_success_probability_for_all_eight_rows':True})
    result=contract.serial({'member_id':member['member_id'],**member['ids'],
        'coherent_response_receipt_sha256':contract.digest(inherited),
        'physical_floor_payload_sha256':contract.digest(floor),'physical_density_floor_factored':density,
        'spectral_sum_rules':sum_rule_control(),'original_signed_success_channels':channels,
        'actual_first_graded_action':neighbor_graded_control(),
        'original_graded_propagation_coefficient_lower':{'physical_density_floor':density,'multiplier':F(1,16)},
        'finite_complete_graded_probe':finite_time_plan(R,N),'gauge_charge_scope':gauge_charge_control(),
        'native_chirality_pairing':bind_chirality_member(member),
        'graded_probe_net_native_increment':[0]*5,'graded_probe_CAR_parity':0,
        'source_state_reference_and_clock_replaced':False,
        'nonzero_original_four_point_spatial_contrast_proved':True,
        'nonzero_original_odd_two_point_spatial_shell_proved':False,
        'isolated_particle_pole_or_relativistic_dispersion_proved':False,
        'physical_mass_table_GeV':[]})
    return {'schema':'oph.primary_odd_sum_rules.v1','result':result,'result_sha256':contract.digest(result)}


def verify_bound(member,value):
    if value!=bind_member(member):raise ValueError('odd sum rule, original graded record or source binding changed')


def certificate():
    paths=provenance.definition_inputs(SOURCE_INPUTS);pins=provenance.definition_hashes(paths)
    value={'schema':'oph.primary_odd_sum_rules_certificate.v1','members':[bind_member(contract.member(0))],
        'chirality_members':[bind_chirality_member(contract.member(n)) for n in (0,1,2,8)],
        'source_sha256':pins}
    if provenance.definition_hashes(paths)!=pins:raise ValueError('source changed during odd sum-rule replay')
    return contract.serial(value)
