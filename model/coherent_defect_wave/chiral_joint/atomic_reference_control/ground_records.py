"""Complete native-unit-hopping quench records from the atomic reference."""
from itertools import product
import sympy as s
from .lattice import gap_data,hopping_norm,local_norm
from .ground_exchange import second_moment_certificate


def certificate():
    atomic=gap_data();exchange=second_moment_certificate()
    chi_lo,chi_hi=map(s.Rational,exchange['actual_local_ground_spin_coefficients']['chi_interval'])
    tau=hopping_norm();h=atomic['volume']*(local_norm()+3+12)+tau
    assert (tau,h)==(360,936)
    time=1/(2048*h*tau)
    state_error=h*tau*time**2/2
    edge_amplitude_norm=15*time
    probability_error=2*edge_amplitude_norm*state_error+state_error**2
    contrast=time**2*chi_lo/8-2*probability_error
    assert contrast>0
    reference_up=s.diag(0,1,0,1)
    assert reference_up**2==reference_up and reference_up+(s.eye(4)-reference_up)==s.eye(4)
    rows=[]
    for reference_spins in product((0,1),repeat=2):
        aligned=reference_spins[0]==reference_spins[1]
        moment=((15-chi_hi)/16,(15-chi_lo)/16) if aligned else ((15+chi_lo)/16,(15+chi_hi)/16)
        success=((time**2*moment[0]-probability_error)/4,
                 (time**2*moment[1]+probability_error)/4)
        failure=(s.Rational(1,4)-success[1],s.Rational(1,4)-success[0])
        assert 0<success[0]<success[1]<s.Rational(1,4) and 0<failure[0]<failure[1]<s.Rational(1,4)
        rows.append({'retained_reference_spin_record':list(reference_spins),
            'conditional_source_ground_spins':[1-spin for spin in reference_spins],
            'reference_record_probability':'1/4',
            'exact_edge_second_moment_interval':[str(v) for v in moment],
            'original_reference_and_link_success_probability_interval':[str(v) for v in success],
            'original_reference_and_link_failure_probability_interval':[str(v) for v in failure],
            'postmeasurement_states':'retain both complete source/reference states without success renormalization'})
    low=sum(s.Rational(row[name][0]) for row in rows for name in
            ('original_reference_and_link_success_probability_interval','original_reference_and_link_failure_probability_interval'))
    high=sum(s.Rational(row[name][1]) for row in rows for name in
             ('original_reference_and_link_success_probability_interval','original_reference_and_link_failure_probability_interval'))
    assert low<=1<=high
    return {'Hamiltonian':'actual full PW H=H0+T, all coefficients one',
            'hopping_coefficient':'1','atomic_source_ground_reference_is_stationary_for_this_H':False,
            'initial_source_density':'P0/2^8',
            'retained_purification':'product of eight source-ground-spin/reference-spin singlets',
            'reference_at_each_cell':'two gauge-neutral fermionic modes of spin 1/2 and raw L=-1, full four-state Fock space',
            'reference_Hamiltonian':'0','local_total_parity':'even','local_total_raw_L':'0',
            'preparation_controller':'swap full vacuum with normalized local source-ground/reference singlet, identity on full complement',
            'preparation_controller_preserves_Gauss_spin_L_and_parity':True,
            'read_reference_spins_at_cells':[0,1],'read_axis':'z',
            'complete_reference_binary_effects':'label0=n_up; label1=I-n_up on the full four-state reference Fock space',
            'reference_labels_on_prepared_one_particle_subspace':['up','down'],
            'selected_actual_oriented_edge':2,
            'success_effect':'selected edge excited, every other link blank, identity on all matter/scalar/reference occupations',
            'failure_effect':'identity minus success on the complete carrier',
            'source_time':str(time),'centered_full_H_norm_upper_bound':str(h),
            'full_first_residual_norm_upper_bound':str(tau),
            'full_second_state_derivative_norm_upper_bound':str(h*tau),
            'first_order_state_remainder_upper_bound':str(state_error),
            'conditional_link_probability_remainder_upper_bound':str(probability_error),
            'antialigned_minus_aligned_conditional_success_lower_bound':str(contrast),
            'antialigned_minus_aligned_original_joint_success_lower_bound':str(contrast/4),
            'all_eight_reference_and_detector_outcomes':rows,
            'complete_probability_enclosure':[str(low),str(high)],
            'all_local_quantum_ground_occupations_and_higher_order_states_retained':True,
            'unit_hopping_stationary_band_theorem_imported':False,
            'interaction_witness_identified_with_a_particle_mass':False}
