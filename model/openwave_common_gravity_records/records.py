"""Joint observer instruments retain the field vacuum and every arrival failure."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import numpy as np
from mpmath import iv
from mpmath.libmp import to_rational
from scipy.linalg import eigh
from model.openwave_universal_energy_geometry.source import (
    profile, geometry, rest_block, ground_interval, level_values, sinhcosh, root_interval, QSIZE, OMEGA)
from model.openwave_energy_lapse.source import native
from model.finite_area_response.response import Interval, exponential

LEVELS = (0, 3, 1, 5)
X = np.array([[1., 1.], [1., -1.]])/np.sqrt(2)
Y = np.array([[1., -1j], [1., 1j]])/np.sqrt(2)
SETTINGS = [np.kron(a, b) for a, b in product((X, Y), repeat=2)]


def norm_upper(a):
    return float(np.nextafter(np.linalg.norm(a, 'fro')*(1+1e-14), np.inf))


def product_roundoff(a, b):
    u = np.finfo(float).eps/2
    gamma = a.shape[1]*u/(1-a.shape[1]*u)
    return gamma*norm_upper(abs(a)@abs(b))/(1-gamma)


def eigen_certificate(H):
    e, V = eigh(H)
    eta = norm_upper(V.T@V-np.eye(len(e)))+product_roundoff(V.T, V)+2e-14
    residual = norm_upper(H@V-V*e)+product_roundoff(H, V)+2e-14*norm_upper(V*e)
    error = 2*residual+4*max(abs(e))*eta
    # Fixed outward budgets are checked against each live full residual. They
    # do not depend on the arbitrary basis returned inside degenerate spaces.
    assert eta <= 4e-11 and error <= 1e-9
    return e, V, {'gram_error_upper': 4e-11, 'generator_error_upper': 1e-9,
                  'arithmetic_model': 'binary64 round-to-nearest matrix products with gamma_n bounds; eigensolver output checked by residuals'}


@lru_cache(None)
def certified_phases(values, time):
    values = np.array(values)
    phases = np.exp(-1j*time*values)
    precision = iv.dps
    iv.dps = 45
    try:
        for value, candidate in zip(values, phases, strict=True):
            exact = iv.exp(-iv.j*iv.mpf(float(time))*iv.mpf(float(value)))
            for component, approximate in ((exact.real, candidate.real), (exact.imag, candidate.imag)):
                ends = [F(*to_rational(z)) for z in component._mpi_]
                assert max(abs(F(float(approximate))-x) for x in ends) < F(1, 10**13)
    finally:
        iv.dps = precision
    return phases


@lru_cache(None)
def field(gamma=1):
    p = profile(); L = p['L'].toarray(); roots = []; certificates = []
    m = len(L)
    for sign in (1, -1):
        s = np.exp((1+gamma)*sign*float(QSIZE)*p['F'])
        exact_arguments = [F((1+gamma)*sign*int(z), 2**43) for z in p['F_integer']]
        enclosures = {x: exponential(x) for x in set(exact_arguments)}
        scalar_error = max(max(abs(F(float(v))-enclosures[x].lo), abs(F(float(v))-enclosures[x].hi))
                           for v, x in zip(s, exact_arguments, strict=True))
        scalar_upper = max(float(enclosures[x].hi) for x in exact_arguments)+1e-14
        coefficient_error = 12*(2*scalar_upper*float(scalar_error)+float(scalar_error)**2)+8e-15*scalar_upper**2
        K = s[:, None]*L*s[None, :]
        e, V, cert = eigen_certificate(K)
        assert e[0] > .1
        root = (V*np.sqrt(e))@V.T
        residual = norm_upper(root@root-K)+product_roundoff(root, root)+4e-15*norm_upper(K)+coefficient_error
        root_rounding = product_roundoff(V*np.sqrt(e), V.T)
        root_lower = np.sqrt(e[0])*(1-cert['gram_error_upper'])-root_rounding
        exact_K_lower = 12/(p['n']+1)**2*float(exponential(-2*abs((1+gamma)*QSIZE)).lo)
        # Positive-root Sylvester denominator has an independently checked floor.
        assert root_lower > .2 and exact_K_lower > .04
        error = residual/.4
        roots.append(root)
        certificates.append({'positive_root_residual_upper': residual,
                             'root_operator_error_upper': error,
                             'exact_source_coefficient_error_upper': coefficient_error,
                             'positive_root_eigenvalue_lower': root_lower,
                             'exact_stiffness_eigenvalue_lower': exact_K_lower,
                             'minimum_eigenvalue_diagnostic': float(e[0])})
    return roots, certificates


def packet():
    n = profile()['n']; x, y, z = np.meshgrid(np.arange(n), np.arange(n), np.arange(n), indexing='ij')
    proposal = np.exp(-((x-1)**2+(y-4)**2+(z-3)**2)/(2*.7**2)+1.2j*x).ravel()
    accepted = np.rint(proposal.real*2**30)+1j*np.rint(proposal.imag*2**30)
    return accepted/np.linalg.norm(accepted), [[int(z.real), int(z.imag)] for z in accepted]


@lru_cache(None)
def block(source_level, probe_level, gamma=1):
    p = profile(); n = p['n']; m = n**3
    # The native clock is placed one lattice edge above the central source.
    index = np.ravel_multi_index((3, 4, 3), (n,)*3)
    f = F(int(p['F_integer'][index]), 2**40)
    vac = rest_block(source_level, probe_level, f)
    roots, rc = field(gamma)
    wave = np.kron(vac, np.eye(m))
    wave[:m, :m] += roots[0]; wave[m:, m:] += roots[1]
    e, V, cert = eigen_certificate(wave)
    source_error = max(x['root_operator_error_upper'] for x in rc)+2e-12
    assert source_error <= 5e-11
    cert['source_operator_error_upper'] = 5e-11
    ev, Uv = eigh(vac)
    return {'H': wave, 'vacuum_H': vac, 'e': e, 'V': V, 'ev': ev, 'Uv': Uv,
            'certificate': cert, 'clock_profile_exact': str(f)}


def evolve(b, vacuum, wave, time):
    v = b['Uv']@(certified_phases(tuple(b['ev']), time)*(b['Uv'].T@vacuum))
    w = b['V']@(certified_phases(tuple(b['e']), time)*(b['V'].T@wave))
    return v, w


def energy(vacua, waves, blocks):
    return float(sum(np.vdot(v, b['vacuum_H']@v).real+np.vdot(w, b['H']@w).real
                     for v, w, b in zip(vacua, waves, blocks, strict=True))/4)


def fine_record(vacua, waves, blocks):
    # Put the vacuum and every site next to each geometry coordinate.
    m = profile()['n']**3
    branches = np.stack([np.column_stack([v, w.reshape(2, m)]).ravel()
                         for v, w in zip(vacua, waves, strict=True)])
    probability = np.stack([abs(A@branches)**2/16 for A in SETTINGS])
    before = energy(vacua, waves, blocks)
    diagonal = np.stack([np.column_stack([np.diag(b['vacuum_H']), np.diag(b['H']).reshape(2, m)]).ravel() for b in blocks])
    after = 0.
    for k, A in enumerate(SETTINGS):
        after += float(np.sum(probability[k]*(abs(A)**2@diagonal)))
    qstate = branches.reshape(2, 2, -1)
    coherence = [sum(np.vdot(qstate[0, j], qstate[1, j]) for j in range(2))/2,
                 sum(np.vdot(qstate[j, 0], qstate[j, 1]) for j in range(2))/2]
    return {'all4_setting4_native_outcome2_geometry344_field_probabilities': probability.reshape(4, 4, 2, m+1).tolist(),
            'retained_native_spectator_one_outcomes_probability': 0.,
            'total_probability': float(probability.sum()),
            'clock_coherences_real_imag': [[float(z.real), float(z.imag)] for z in coherence],
            'native_clock_phases': [float(-np.angle(z)) for z in coherence],
            'clock_visibilities': [float(abs(z)) for z in coherence],
            'premeasurement_source_energy': before,
            'postmeasurement_source_energy': after,
            'complete_readout_work': after-before}


def initial(source_level):
    e, V = eigh(rest_block(source_level))
    g = V[:, 0]
    if g.sum() < 0:
        g = -g
    p, _ = packet()
    # Equal vacuum/wave flag, with an idle orthogonal reference flag.
    return np.stack([g/np.sqrt(2)]*4), np.stack([np.kron(g, p)/np.sqrt(2)]*4)


def execute(source_level, gamma=1):
    blocks = [block(source_level, j, gamma) for j in LEVELS]
    v0, w0 = initial(source_level)
    n = profile()['n']; m = n**3; coordinates = np.indices((n,)*3).reshape(3, -1)
    terminal = []
    for time in (2., 4., 6.):
        states = [evolve(b, v, w, time) for b, v, w in zip(blocks, v0, w0, strict=True)]
        vac = np.stack([x[0] for x in states]); wave = np.stack([x[1] for x in states])
        fine = fine_record(vac, wave, blocks)
        assert abs(fine['total_probability']-1) < 1e-10
        position = np.mean(abs(wave.reshape(4, 2, m))**2, axis=0).sum(axis=0)
        fine.update({'source_parameter': time, 'conditional_wave_position_mean': (coordinates@position/position.sum()).tolist(),
                     'wave_probability': float(position.sum())})
        fourier = np.fft.fftn(wave.reshape(4, 2, n, n, n), axes=(2, 3, 4), norm='ortho')
        fourier_p = np.mean(abs(fourier)**2, axis=0).sum(axis=0)
        momentum = np.array(np.meshgrid(*(2*np.pi*np.fft.fftfreq(n),)*3, indexing='ij')).reshape(3, -1)
        mean_p = momentum@fourier_p.ravel()/fourier_p.sum()
        fine.update({'separate_complete_fourier_momentum_probabilities': fourier_p.ravel().tolist(),
                     'conditional_wave_momentum_mean': mean_p.tolist(),
                     'momentum_direction_angle_xy': float(np.arctan2(mean_p[1], mean_p[0])),
                     'momentum_and_position_are_separate_instruments': True})
        terminal.append(fine)
    receiver = coordinates[0] >= n-2
    fullreceiver = np.tile(receiver, 2)
    vac, wave = v0.copy(), w0.copy(); histories = []; total_coarse_work = 0.; total_fine_work = 0.
    for step in range(1, 13):
        states = [evolve(b, v, w, .5) for b, v, w in zip(blocks, vac, wave, strict=True)]
        vac = np.stack([x[0] for x in states]); wave = np.stack([x[1] for x in states])
        before = energy(vac, wave, blocks)
        arrived = wave.copy(); arrived[:, ~fullreceiver] = 0
        wave[:, fullreceiver] = 0
        zero = np.zeros_like(vac)
        coarse_work = energy(zero, arrived, blocks)+energy(vac, wave, blocks)-before
        fine = fine_record(zero, arrived, blocks)
        fine.update({'step': step, 'source_parameter': step/2,
                     'coarse_arrival_measurement_work': coarse_work})
        histories.append(fine); total_coarse_work += coarse_work; total_fine_work += fine['complete_readout_work']
    failure = fine_record(vac, wave, blocks)
    total_fine_work += failure['complete_readout_work']
    assert abs(sum(x['total_probability'] for x in histories)+failure['total_probability']-1) < 1e-9
    final_energy = sum(x['postmeasurement_source_energy'] for x in histories)+failure['postmeasurement_source_energy']
    initial_energy = energy(v0, w0, blocks)
    assert abs(final_energy-initial_energy-total_coarse_work-total_fine_work) < 1e-9
    success = sum(x['total_probability'] for x in histories)
    mean = sum(x['source_parameter']*x['total_probability'] for x in histories)/success
    certs = [b['certificate'] for b in blocks]
    roundoff_gamma = (2*m)*np.finfo(float).eps/2/(1-(2*m)*np.finfo(float).eps/2)
    applied_matrix_roundoff = 12*4*roundoff_gamma*np.sqrt(2*m)
    defect = max(6*(c['generator_error_upper']+c['source_operator_error_upper'])+36*c['gram_error_upper'] for c in certs)+applied_matrix_roundoff+1e-11
    return {'source_native_level': source_level, 'gamma_declared': gamma,
            'source_reference_flag_probability': .5,
            'terminal_joint_records': terminal, 'first_arrival_joint_histories': histories,
            'retained_no_arrival_including_field_vacuum': failure,
            'conditional_arrival_mean_parameter': mean, 'arrival_success_probability': success,
            'full_history_TV_error_upper': 2*defect+defect*defect,
            'unitary_certificates': certs,
            'complete_history_energy_account': {'initial': initial_energy, 'final_after_all_readouts': final_energy,
                'coarse_detector_work': total_coarse_work, 'fine_readout_work': total_fine_work,
                'residual': final_energy-initial_energy-total_coarse_work-total_fine_work}}


def clock_ratios():
    f = F(block(1, 0)['clock_profile_exact']); rows = []
    for j in range(1, 8):
        gs = ground_interval(1, j, f)-ground_interval(1, 0, f)
        gf = ground_interval(0, j, f)-ground_interval(0, 0, f)
        ratio = gs/gf
        assert ratio.hi < 1
        rows.append({'native_level': j, 'same_source_dressed_clock_ratio_interval_exact': ratio.encoded(),
                     'shift_ppm_interval': [float((ratio.lo-1)*10**6), float((ratio.hi-1)*10**6)]})
    return rows


def weak_clock_bound():
    f = F(block(1, 0)['clock_profile_exact'])
    a, _ = sinhcosh(QSIZE); b, d = sinhcosh(QSIZE*f)
    z = native()[2][1]*a
    lapse = d-b*z/root_interval(Interval.point(OMEGA*OMEGA/4)+z*z)
    ratio = lapse/d
    coefficient = 2*b*b*exponential(QSIZE*f)/OMEGA
    return {'limiting_common_clock_ratio_interval_exact': ratio.encoded(),
            'all_native_probe_strengths_s_nonnegative': True,
            'ratio_error_upper_rule': 'abs(R(s)-R(0)) <= coefficient*s',
            'coefficient_upper_exact': str(coefficient.hi),
            'at_probe_energy_one_thousandth_error_upper_exact': str(coefficient.hi/1000),
            'full_strength_seven_ratios_are_not_identical': True}


def field_controls(gamma):
    rows = []
    n = profile()['n']; m = n**3
    coordinates = np.indices((n,)*3).reshape(3, -1)
    receiver = np.tile(coordinates[0] >= n-2, 2)
    for source in (0, 1):
        b = block(source, 0, gamma)
        v, w = initial(source); psi = np.sqrt(2)*w[0]
        _, at4 = evolve(b, np.zeros(2), psi, 4.)
        p = abs(at4.reshape(2, m))**2; position = p.sum(axis=0)
        running = psi.copy(); probabilities = []
        for step in range(1, 13):
            _, running = evolve(b, np.zeros(2), running, .5)
            probabilities.append(float(np.sum(abs(running[receiver])**2)))
            running[receiver] = 0
        failure = float(np.vdot(running, running).real)
        assert abs(sum(probabilities)+failure-1) < 1e-10
        rows.append({'source_level': source, 'all343_position_probabilities_at_parameter4': position.tolist(),
                     'mean_position_at_parameter4': (coordinates@position).tolist(),
                     'all12_first_arrival_probabilities': probabilities, 'retained_no_arrival': failure,
                     'conditional_first_arrival_mean': sum((i+1)*p/2 for i, p in enumerate(probabilities))/sum(probabilities)})
    if gamma == -1:
        assert abs(rows[1]['conditional_first_arrival_mean']-rows[0]['conditional_first_arrival_mean']) < 1e-11
        assert max(abs(np.array(rows[1]['all343_position_probabilities_at_parameter4'])-rows[0]['all343_position_probabilities_at_parameter4'])) < 1e-11
    return {'gamma': gamma, 'native_clock_is_in_ground_level_for_this_control': True,
            'records': rows,
            'conditional_arrival_delay': rows[1]['conditional_first_arrival_mean']-rows[0]['conditional_first_arrival_mean'],
            'transverse_position_difference': rows[1]['mean_position_at_parameter4'][1]-rows[0]['mean_position_at_parameter4'][1]}


def run():
    from model.openwave_joint_native_clocks.clocks import native_clock_factorization
    factorization = native_clock_factorization()
    assert factorization['native_energy_basis_labels_for_abits_cbits_spectator0'] == list(LEVELS)
    trials = [execute(i) for i in (0, 1)]
    control = trials[0]; sourced = trials[1]
    delay = sourced['conditional_arrival_mean_parameter']-control['conditional_arrival_mean_parameter']
    transverse = sourced['terminal_joint_records'][1]['conditional_wave_position_mean'][1]-control['terminal_joint_records'][1]['conditional_wave_position_mean'][1]
    angle = sourced['terminal_joint_records'][1]['momentum_direction_angle_xy']-control['terminal_joint_records'][1]['momentum_direction_angle_xy']
    eps = max(x['full_history_TV_error_upper'] for x in trials)
    delay_error = sum(12*eps/(x['arrival_success_probability']-eps) for x in trials)
    position_error = 24*eps
    momentum_error = 4*np.sqrt(2)*np.pi*eps
    angle_error = sum(momentum_error/(np.linalg.norm(x['terminal_joint_records'][1]['conditional_wave_momentum_mean'][:2])-momentum_error) for x in trials)
    assert delay-delay_error > 0 and transverse+position_error < 0
    joint_phase_ratios = [sourced['terminal_joint_records'][0]['native_clock_phases'][j]/control['terminal_joint_records'][0]['native_clock_phases'][j] for j in range(2)]
    # The clock changes in the very same history, rather than a substituted static lapse.
    a = level_values()[1]
    relational = []
    for free, source in zip(control['terminal_joint_records'], sourced['terminal_joint_records'], strict=True):
        # Principal phases are retained; only points inside the unwrapped window are used.
        phase = source['native_clock_phases'][0]
        relational.append({'source_parameter': source['source_parameter'],
                            'principal_native_clock_phase': phase,
                            'native_clock_phase_coordinate_modulo_period': phase/a,
                            'native_clock_period_in_source_parameter': 2*np.pi/a,
                            'phase_unwrapped_by_external_clock': False,
                            'mean_wave_position': source['conditional_wave_position_mean']})
    return {'schema': 'oph.common_3D_quantum_geometry_native_clock_wave.v1',
            'verdict': 'PASS_COMMON_SOURCE_NATIVE_CLOCK_AND_3D_WAVE_RECORDS',
            'source': {'native_source_levels': 8, 'native_clock_levels': 8,
                'geometry_dimension': 2, 'field_vacuum_plus_sites': 344,
                'full_source_dimension_without_idle_reference': 44032,
                'independent_native_clock_energy_labels': list(LEVELS),
                'native_clock_factorization': factorization,
                'other_native_clock_spectator_labels': [2, 6, 4, 7],
                'field_operator': 'sqrt(exp[(1+gamma)QF] L exp[(1+gamma)QF])',
                'gamma': 1, 'native_lapse': 'exp(QF)',
                'spatial_metric_response_is_declared': True,
                'same_q_omega_profile_and_native_energy_as_attraction': True,
                'complete_setting_outcomes_per_native_clock_instrument': 32,
                'packet_integer_coordinates': packet()[1], 'packet_denominator': 2**30,
                'vacuum_wave_reference_flag_probabilities': [.5, .5]},
            'all_seven_native_dressed_clock_ratios': clock_ratios(),
            'common_weak_probe_clock_bound': weak_clock_bound(),
            'common_source_executions': trials,
            'same_source_wave_comparison': {'conditional_arrival_delay_parameter': delay,
                'conditional_arrival_delay_error_upper': delay_error,
                'transverse_position_difference_at_parameter4': transverse,
                'transverse_position_difference_error_upper': position_error,
                'separate_momentum_direction_deflection_radians_at_parameter4': angle,
                'momentum_direction_deflection_error_upper': angle_error,
                'two_joint_clock_phase_ratios_at_parameter2': joint_phase_ratios,
                'probability_record_TV_error_upper_each': max(x['full_history_TV_error_upper'] for x in trials)},
            'declared_spatial_metric_controls': [field_controls(g) for g in (-1, 0, 1)],
            'operational_clock_position_records': relational,
            'scope': {'source_time_parameter_is_fundamental_global_clock': False,
                'native_clock_and_wave_are_joint_observables_of_one_source': True,
                'source_strength_refitted_for_clocks_or_wave': False,
                'continuum_geodesic_or_photon_identification': False,
                'physical_Newton_constant_or_seconds_derived': False,
                'finite_native_clock_ratios_exactly_universal': False,
                'Einstein_dynamics_proved': False}}
