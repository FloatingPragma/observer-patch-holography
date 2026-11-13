"""Coupling-family certificates reusing the full-complement spectral proof."""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import json
import numpy as np
from model.openwave_coupled_haar_source.source import matrices,hamiltonian,tail_floor
from model.openwave_pure_su2_source.source import coefficients
from model.openwave_colour_string.fixed_witness import spectral_enclosure
from model.openwave_su2_spectral.spectrum import sqrt_upper
from model.openwave_coupled_haar_spectrum.certificate import source_digest,rounded_operator,strings,display,interval_product,ROUND_DENOMINATOR
HERE=Path(__file__).resolve().parent
@lru_cache(None)
def witnesses():
    return json.loads((HERE/'numerical_witnesses.json').read_text())

@lru_cache(None)
def calculate(key):
    packet = witnesses()[key]
    cap, g2, shared = packet['cap'], F(packet['g_squared_exact']), packet['shared_edge']
    assert type(cap) is int and type(shared) is bool
    assert packet['spacing_exact'] == '1'
    basis, (W1, W2) = matrices(cap)
    H = hamiltonian(cap, g2, shared=shared)
    assert packet['exact_source_matrix_sha256'] == source_digest(H)
    assert packet['ordered_spin_network_labels'] == [list(x) for x in basis]
    finite = spectral_enclosure(H, packet['complete_finite_frame'])
    f_intervals = [tuple(map(F, row)) for row in finite['all_eigenvalue_intervals']]
    alpha, beta = coefficients(g2, F(1))
    floor = tail_floor(cap, g2)
    ceiling = f_intervals[3][1]
    assert ceiling < floor
    offblock2 = 4*beta*beta
    shift = offblock2/(floor-ceiling)
    complete = [(max(F(0), lo-shift), hi) for lo, hi in f_intervals[:4]]
    assert complete[1][0] > complete[0][1]

    # Numerical coordinates are proof witnesses; every residual is recomputed.
    D = packet['ground_denominator']
    assert type(D) is int and D == 2**48
    raw = packet['ground_integer_coordinates']
    assert len(raw) == len(basis) and all(type(x) is str and str(int(x)) == x and abs(int(x)) <= D for x in raw)
    u = np.array([int(x) for x in raw], dtype=object)
    norm = sum(int(x)**2 for x in u)
    assert norm > 0
    HD, operator_error = rounded_operator(H)
    w1, error1 = rounded_operator(W1)
    w2, error2 = rounded_operator(W2)
    S = ROUND_DENOMINATOR
    boundary_indices = [i for i, (a, b, _) in enumerate(basis) if a == cap or b == cap]
    low_trials = packet['low_trial_integer_coordinates']
    assert len(low_trials) == 3 and low_trials[0] == raw
    windows = []
    for k, coordinates in enumerate(low_trials):
        assert len(coordinates) == len(basis) and all(type(x) is str and str(int(x)) == x and abs(int(x)) <= D for x in coordinates)
        trial = np.array([int(x) for x in coordinates], dtype=object)
        trial_norm = sum(int(x)**2 for x in trial)
        assert trial_norm > 0
        center = F(packet['complete_finite_frame']['dyadic_eigenvalues'][k])/F(2**40)
        trial_residual = [F(int(x), S)-center*int(y) for x, y in zip(HD@trial, trial)]
        residual_upper = sqrt_upper(sum(x*x for x in trial_residual)/trial_norm)+operator_error
        residual_upper += 2*beta*sqrt_upper(F(sum(int(trial[i])**2 for i in boundary_indices), trial_norm))
        windows.append((center-residual_upper, center+residual_upper))
    separated = all(windows[k][1] < windows[k+1][0] for k in range(2)) and windows[-1][1] < complete[3][0]
    if shared:
        assert separated
    if separated:
        # Each full residual puts a spectral point in its own disjoint window.
        # The full fourth-level floor leaves room for exactly those three points.
        complete[:3] = [(max(old[0], new[0]), min(old[1], new[1])) for old, new in zip(complete[:3], windows)]
        assert all(lo <= hi for lo, hi in complete[:3])
    hu = HD@u
    rayleigh = F(sum(int(x)*int(y) for x, y in zip(u, hu)), S*norm)
    mean_interval = (rayleigh-operator_error, rayleigh+operator_error)
    mu = rayleigh
    residual = [F(int(x), S)-mu*int(y) for x, y in zip(hu, u)]
    finite_residual = sqrt_upper(sum(x*x for x in residual)/norm)+operator_error
    boundary_norm = sqrt_upper(F(sum(int(u[i])**2 for i in boundary_indices), norm))
    full_residual = finite_residual+2*beta*boundary_norm
    separation = complete[1][0]-mu
    assert separation > 0
    distance = full_residual/separation
    assert distance < F(1, 100)
    # Temple's inequality uses the separately certified full first excitation.
    temple = (max(F(0), mean_interval[0]-full_residual**2/(complete[1][0]-mean_interval[1])), mean_interval[1])
    assert temple[0] <= temple[1] < complete[1][0]
    gap = (complete[1][0]-temple[1], complete[1][1]-temple[0])

    first, second = w1@u, w2@u
    m1 = F(sum(int(x)*int(y) for x, y in zip(u, first)), S*norm)
    m2 = F(sum(int(x)*int(y) for x, y in zip(u, second)), S*norm)
    m12 = F(sum(int(x)*int(y) for x, y in zip(first, second)), S*S*norm)
    # Outside W1 and W2 vectors have disjoint loop-spin labels, so this cross
    # product includes the complete, untruncated mixed moment.
    err1, err2 = error1+2*distance, error2+2*distance
    err12 = error1+error2+error1*error2+2*distance
    I1, I2, I12 = (m1-err1, m1+err1), (m2-err2, m2+err2), (m12-err12, m12+err12)
    prod = interval_product(I1, I2)
    covariance = (I12[0]-prod[1], I12[1]-prod[0])
    if shared:
        assert covariance[0] > 0
    else:
        assert covariance[0] <= 0 <= covariance[1]
    moment2 = F(sum(int(x)**2 for x in first), S*S*norm)
    qerror2 = 2*error1+error1*error1
    second_moment = (moment2-qerror2-2*distance,
                     moment2+qerror2+boundary_norm**2+2*distance)
    sq = interval_product(I1, I1)
    variance = (max(F(0), second_moment[0]-sq[1]), min(F(1), second_moment[1]-sq[0]))
    assert variance[0] > 0

    eta = F(1, 2)
    records = []
    centers = []
    for sign1, sign2 in ((-1,-1),(-1,1),(1,-1),(1,1)):
        center = (1+sign1*eta*m1+sign2*eta*m2+sign1*sign2*eta**2*m12)/4
        error = (eta*(err1+err2)+eta**2*err12)/4
        interval = (center-error, center+error)
        assert interval[0] > 0 and interval[1] < 1
        centers.append(center)
        records.append({'outcome': [sign1, sign2], 'probability_interval_exact': strings(interval),
                        'probability_interval': display(interval)})
    assert sum(centers) == 1
    result = {'g_squared_exact': str(g2), 'spacing_exact': '1', 'shared_edge': shared,
              'diagnostic_max_twice_loop_spin': cap, 'retained_complete_dimension': len(basis),
              'complete_finite_spectrum': finite,
              'positive_full_complement_floor_exact': str(floor),
              'full_offblock_norm_squared_upper_exact': str(offblock2),
              'Schur_shift_upper_exact': str(shift),
              'full_native_residual_spectral_windows_exact': [strings(row) for row in windows],
              'first_three_levels_identified_by_disjoint_windows_and_fourth_floor': separated,
              'fourth_level_exclusion_floor_exact': str(complete[3][0]),
              'full_Haar_first_four_energy_intervals_exact': [strings(row) for row in complete],
              'ground_energy_interval_exact': strings(temple), 'ground_energy_interval': display(temple),
              'full_Haar_gap_interval_exact': strings(gap), 'full_Haar_gap_interval': display(gap),
              'full_source_trial_residual_upper_exact': str(full_residual),
              'ground_projector_distance_upper_exact': str(distance),
              'boundary_trial_norm_upper_exact': str(boundary_norm),
              'Wilson_means_intervals_exact': [strings(I1), strings(I2)],
              'Wilson_mixed_moment_interval_exact': strings(I12),
              'connected_Wilson_covariance_interval_exact': strings(covariance),
              'connected_Wilson_covariance_interval': display(covariance),
              'Wilson_variance_interval_exact': strings(variance),
              'nonzero_centered_physical_energy_upper_exact': str(alpha/variance[0]),
              'complete_joint_observer_instrument': {
                  'response_eta_exact': str(eta), 'effects': '(1+s eta W1)(1+t eta W2)/4 for all s,t=+-1',
                  'Kraus_rule': 'positive square root of each joint effect on the full Haar physical space',
                  'complete_effect_sum': 'I', 'effect_floor_exact': '1/16',
                  'external_quantum_reference_retained': True, 'four_outcomes': records,
                  'record_covariance_interval_exact': strings(tuple(eta**2*x for x in covariance)),
                  'nonselective_ground_preparation_energy_gain_interval_exact': ['0', str(alpha/6)]},
              'scope': {'graph_is_full_infinite_Haar_source': True, 'cap_is_only_a_proof_coordinate': True,
                        'physical_GeV_mass': False, 'infinite_spatial_volume': False, 'four_dimensional_continuum': False}}
    internal = {'H': H, 'HD': HD, 'operator_error': operator_error, 'W1': W1, 'W2': W2,
                'w1': w1, 'w2': w2, 'error1': error1, 'error2': error2,
                'u': u, 'norm': norm, 'mu': mu, 'boundary_norm': boundary_norm,
                'basis': basis, 'm1': m1, 'distance': distance, 'ground_interval': temple,
                'packet': packet, 'alpha': alpha, 'beta': beta}
    return result, internal

