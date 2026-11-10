"""Declared finite programmable pulses with the autonomous clock retained."""
from fractions import Fraction as F
import sympy as s

from .source import admit, constituent_bounds


def exact_noncommuting_pulse_control():
    """A two-state check whose desired intervention does not commute with H."""
    I = s.eye(2); X = s.Matrix([[0, 1], [1, 0]]); Z = s.diag(1, -1)
    free = (I-s.I*Z)/s.sqrt(2)
    pulse = (I-s.I*X)/s.sqrt(2)
    desired = free.H*pulse
    H, L = s.pi*Z/4, s.pi*X/4
    control = L-H
    assert H*desired != desired*H
    assert (free*desired-pulse).applyfunc(s.simplify) == s.zeros(2)
    assert H+control == L
    assert (desired.H*desired-I).applyfunc(s.simplify) == s.zeros(2)
    assert (pulse.H*pulse-I).applyfunc(s.simplify) == s.zeros(2)
    return {'noncommuting_intervention_control_passed': True,
            'source_H': '(pi/4) sigma_z', 'pulse_duration': '1',
            'total_pulse_H': '(pi/4) sigma_x',
            'actual_pulse_propagator': 'exp(-i H) U_desired',
            'native_clock_is_not_silently_paused': True}


def resource_bounds(side, depth, attempts=1):
    if side < 3 or depth < 0 or attempts < 1:
        raise ValueError('finite cubic member, nonnegative depth and positive attempts required')
    parent, _ = admit()
    U = F(parent['intervals']['native_selected_E_interval'][1])
    B = 4*side**3; trials = B*2**depth
    native = F(constituent_bounds()['raw_full_kinematic_native_H_norm_upper'])
    # ||K_a|| <= sqrt(12)/a <4/a. On every full source sector,
    # sum a_j^dagger a_j <=4V. Here a=8pi/(3Nm), pi>3.
    kinetic = F(side, 2)*U*B
    whole_norm = attempts*(trials*(native+120)+kinetic)
    # Preparation tau=1. A detector tau=T/4 has pi/tau=8Omega<=10m.
    preparation_control = whole_norm+F(22, 7)
    detector_control = whole_norm+10*U
    return {'complete_raw_source_H_norm_upper': str(whole_norm),
            'domain': 'all main and retained partner sources, full kinematic complements, every attempt; zero-energy apparatus registers',
            'norm_formula': 'R [4V*2^d (B_native+120)+(N/2) U*4V]',
            'raw_constant_ground_energy_subtraction_is_not_used': True,
            'programmable_Hermitian_log_branch': 'each eigenphase in (-pi,pi], L=i Log(exp(-i H_ext tau) U)/tau',
            'exact_pulse_identity': 'exp(-i tau (H_ext+V_control))=exp(-i H_ext tau) U',
            'preparation_pulse_duration_source_units': '1',
            'each_detector_pulse_duration': 'T/4, T=pi/(2 Omega_a(p))',
            'bare_wait_between_detector_pulses': '3T/4',
            'total_elapsed_time': '1+5T/4',
            'whole_preparation_pulse_control_norm_upper': str(preparation_control),
            'each_detector_pulse_control_norm_upper': str(detector_control),
            'net_bare_source_energy_change_absolute_upper_per_pulse': str(2*whole_norm),
            'net_bare_source_energy_change_absolute_upper_entire_schedule': str(2*whole_norm),
            'sum_absolute_switching_work_upper_three_pulses': str(2*(preparation_control+2*detector_control)),
            'matrix_norm_gate_cost_is_separate_from_work': True,
            'coefficient_error_record_TV_bound': 'sum_j tau_j ||L_implemented,j-L_declared,j||',
            'stronger_programmable_finite_matrix_Hamiltonian_service_declared': True,
            'passive_CE_law_compilation_proved': False,
            'autonomous_local_controller_proved': False,
            'additional_controller_synthesis_workspace_in_capacity': False}
