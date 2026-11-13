"""Exact composite support and the missing interacting spectral transfer."""
from fractions import Fraction as F
from model.openwave_capacity_stress.source import enclosure
from model.openwave_hadron_refinement_flow.certified_rest import sqrt_upper
from model.openwave_strange_shell_current.native import contraction,vertex_words
from model.openwave_strange_shell_current.source import lambda_pole,refined_shell
from .native import channel,separated_hopping
import sympy as sy


def matrix_record(A):return [[str(v)for v in row]for row in A.tolist()]


def support_test():
    """Every literal word creates an anti-u beside one of its output u quarks."""
    records=[]
    for coefficient,ops in vertex_words():
        quark=[i for create,i in ops if create and i<48 and (i//3)%2==0]
        antiquark=[i for create,i in ops if create and i>=48]
        assert len(quark)==len(antiquark)==1
        x=quark[0]//24;y=(antiquark[0]-48)//24
        assert x==y
        records.append((x,y))
    # Every charged-pion state contains d + anti-u. Both output u quarks
    # therefore belong to the proton component in this fixed flavour channel.
    return {'literal_words_checked':len(records),'created_antiquark_and_converted_u_site_equal':True,
        'both_output_u_quarks_belong_to_proton_in_p_pi_minus_channel':True,
        'direct_current_to_disjoint_bare_proton_pion_support_exact':'0',
        'holds_on_complete_initial_uds_relative_space':True,
        'vanishing_is_an_on_shell_interacting_decay_amplitude':False}


def run():
    c=channel();M=c['projected_Gram'];G=c['source_Gram'];total=c['total_Gram']
    values=(G.inv()*M).eigenvals()
    expected={sy.Rational(6,5):1,sy.Rational(241,150):2,
              (277-sy.sqrt(54409))/200:2,(277+sy.sqrt(54409))/200:2}
    assert values==expected
    # The exact 2x2 block discriminant is bounded by 237^2; no floating PSD.
    assert 54409<237**2 and 277+237<520
    assert all(v>=0 for v in(M-G/5).eigenvals())
    assert all(v>=0 for v in(13*G/5-M).eigenvals())
    other=total-M
    assert all(v>0 for v in other.eigenvals())
    cab=contraction()['Cabibbo_squared'];pole=lambda_pole()
    hopping=separated_hopping()
    K=hopping['current_then_hopping_Gram_without_t_Cabibbo']
    assert K[1,1]/G[1,1]==sy.Rational(9,4)
    # Full N5 nearest-neighbour kinetic off-diagonal norm <=10t.
    # The actual normalized second Lambda field has ||W psi||=sqrt(6)c.
    kinetic_amplitude=F(3,2000)-sqrt_upper(F(6))/20000
    assert kinetic_amplitude>0
    kinetic_probability=cab*kinetic_amplitude**2/1296
    q=F(pole['certificate']['full_Q_probability_upper']);z=F(pole['certificate']['prepared_trial_pole_weight_lower'])
    # Tight rational bounds use outward square-root arithmetic.
    s=sqrt_upper(F(54409));a=(F(277)-s)/200;b=(F(277)+s)/200
    assert a>F(1,5) and b<F(13,5)
    lower=a*cab*(1-q);upper=b*cab
    complement_lower=(F(6)-b)*cab*(1-q);complement_upper=(F(7)-a)*cab
    p0=lower/1296;p1=upper/1296
    return {'schema':'oph.native_proton_pion_current_support.v1',
        'verdict':'PASS_EXACT_NATIVE_PROTON_PION_CHANNEL_RESPONSE_WITH_EXPLICIT_SCATTERING_BOUNDARY',
        'exact_CAR_channel':{'raw_occupation_union_dimension':c['raw_union_dimension'],
            'proton_native_spin_half_fields':5,'additional_onsite_proton_translate':1,
            'pion_spin_zero_position_fields':4,'canonical_composite_columns':23,
            'independent_composite_rank':len(c['independent_columns']),
            'literal_native_vertex_words':len(vertex_words()),
            'all_CAR_signs_and_Pauli_exchange_retained':True,
            'fitted_mass_eigenvectors_used_to_define_channel':False,
            'projector':'Pi_C=Y(Y*Y)^(-1)Y*, with exact independent native CAR columns Y.',
            'composite_Gram_matrix':matrix_record(c['basis_Gram']),
            'source_Gram_matrix':matrix_record(G),
            'projected_current_Gram_without_Cabibbo':matrix_record(M),
            'orthogonal_compact_other_current_Gram_without_Cabibbo':matrix_record(other),
            'generalized_eigenvalues_exact':{str(k):v for k,v in expected.items()},
            'generalized_eigenvalue_bounds_exact':[str(a),str(b)],
            'Gram_order_exact':'(1/5) G <= W_P* Pi_C W_P / c^2 <= (13/5) G',
            'complete_Pi_C_and_I_minus_Pi_C_instrument':True},
        'complete_Lambda_pole_transfer':{
            'scope':'Fixed total momentum k=0 source fibre, with all relative Q configurations retained; no normalizable Lambda packet executed.',
            'full_ground_GeV_interval':pole['full_ground_GeV_interval'],
            'full_Q_probability_upper':str(q),'Cabibbo_factor_squared_exact':str(cab),
            'proton_pion_compact_response_norm_squared_interval':[str(lower),str(upper)],
            'proton_pion_compact_response_display':enclosure((lower,upper)),
            'other_compact_N5_response_norm_squared_interval':[str(complement_lower),str(complement_upper)],
            'uniform_positive_for_every_compact_Lambda_coordinate':True,
            'proof':'The local vertex preserves occupied spatial span. Pi_C lies in span<=1; noncompact initial Q cannot cancel this channel output. Its full-pole lower bound pays the entire compact-weight deficit.',
            'normalizable_Lambda_packet_executed':False},
        'complete_probe_records':{
            'native_success_operator':'W/36 on the full uds source; zero on its Fock complement.',
            'native_failure_operator':'sqrt(I-W*W/1296), including the full Fock complement.',
            'terminal_outcomes':['compact native proton-pion field space','orthogonal compact N5 output','noncompact N5 output','probe failure'],
            'proton_pion_probability_interval':[str(p0),str(p1)],
            'other_compact_probability_interval':[str(complement_lower/1296),str(complement_upper/1296)],
            'noncompact_probability_interval':['0',str(cab*q)],
            'failure_probability_interval':[str(1-7*cab/1296-cab*q),str(1-6*cab*(1-q)/1296)],
            'original_P7_input_trial_pole_proton_pion_event_probability_lower':str(z*p0/7),
            'trial_preparation_success_exact':'1/7','trial_preparation_failure_exact':'6/7',
            'nonselective_preparation_work_exact':'0 on the specified P7/7 input',
            'probe_and_channel_readout_work':'Tr[H(Phi(rho)-rho)] retained for the complete instrument; not claimed zero.',
            'channel_filter_commutes_with_interacting_H5':False,
            'unobserved_channels_discarded':False},
        'bare_separated_support_control':support_test(),
        'actual_kinetic_substep':{
            'literal_hopping_words':hopping['literal_hopping_words'],
            'target_fields':'Onsite proton at cell0 and onsite pion at cell1, and their exchanged positions; disjoint supports one cell apart.',
            'target_Gram_matrix':matrix_record(hopping['target_Gram']),
            'current_then_hopping_Gram_without_t_Cabibbo':matrix_record(K),
            'actual_input':'The normalized second native Lambda coordinate, spatially localized at cells0/1; not the rest-pole eigenstate.',
            'input_Gram_norm_squared_exact':str(G[1,1]),
            'current_norm_squared_over_Cabibbo_squared_exact':'6',
            'leading_separated_amplitude_squared_over_t_squared_Cabibbo_squared_exact':'9/4',
            'time_source_units':'1/(1000 t)',
            'full_N5_offdiagonal_kinetic_norm_upper':'10 t',
            'separated_amplitude_lower_over_Cabibbo_exact':str(kinetic_amplitude),
            'original_input_probe_then_kinetic_separated_probability_lower':str(kinetic_probability),
            'original_input_probe_then_kinetic_separated_probability_display':enclosure((kinetic_probability,kinetic_probability)),
            'complete_probe_success_probability_exact':str(cab/216),
            'complete_probe_failure_probability_exact':str(1-cab/216),
            'remaining_success_outcome':'Every N5 state orthogonal to the two disjoint fields; retained with probability success minus separated probability.',
            'full_infinite_N5_kinetic_carrier_retained':True,
            'nonzero_finite_time_response_rigorously_certified':True,
            'complete_autonomous_strong_Hamiltonian_used':False,
            'scope':'One genuine complete kinetic substep of the native source; all other strong terms remain required for autonomous scattering. This is not an on-shell amplitude or lifetime.'},
        'same_source_open_energy_shell':refined_shell(),
        'interacting_spectral_test':{
            'required_operator':'Omega_out* W Psi_Lambda, or the corresponding full H5 spectral density/resolvent boundary value at E_Lambda.',
            'compact_channel_projector_equals_H5_spectral_projector':False,
            'full_H5_on_shell_spectral_density_computed':False,
            'Moller_wave_operator_constructed':False,
            'local_response_and_energy_shell_imply_nonzero_on_shell_response':False,
            'precise_remaining_premise':'Construct the outgoing proton-pion wave operator for the same finite-range strong H5, establish its nonzero on-shell overlap with the actual native weak output, then include the full weak Hamiltonian, vacuum and other channels.'},
        'inherited_inputs':'Electron GeV anchor, strong and axial source coefficients, strange coefficients, Cabibbo angle, physical speed and strong length.',
        'new_mass_or_current_calibration':False,'new_independent_match_with_nature':False,
        'physical_Lambda_decay_amplitude_proved':False,'physical_Lambda_lifetime_proved':False,
        'scope':'Exact source composite-channel response and exact bare separated-support control. The channel is a complete native field subspace, not an interacting asymptotic particle or energy-shell projection.'}
