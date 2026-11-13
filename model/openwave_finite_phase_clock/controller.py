"""A time-independent finite controller, with complete numerical records.

The exact operator identity is independent of the numerical propagator.
Floating trajectory diagnostics are not outward interval certificates.
"""
from fractions import Fraction as F
import numpy as np
from scipy.sparse import csr_matrix, diags, eye, kron
from scipy.sparse.linalg import expm_multiply
from .source import exchange

SIDE = 256
SCALE = 4096
J = 16
PREPARATION = (48, 64)
READOUT = (192, 208)
TIMES = np.linspace(0., 5.7, 115)


def source(mass, coupling=2., *, twirl=False):
    """H/SCALE; physical source time equals the displayed clock parameter/SCALE."""
    N, A, R, *_ = exchange()
    N, A, R = (np.asarray(z, complex) for z in (N, A, R))
    I = np.eye(2)
    lower = np.array([[0., 1.], [0., 0.]])
    n = np.kron(N, np.eye(4))
    r1 = np.kron(np.eye(4), np.kron(R, I))
    r2 = np.kron(np.eye(8), R)
    b1 = np.kron(A, np.kron(lower, I))+np.kron(A.conj().T, np.kron(lower.T, I))
    b2 = np.kron(A, np.kron(I, lower))+np.kron(A.conj().T, np.kron(I, lower.T))
    aa = np.kron(A, np.eye(4))
    x = np.arange(SIDE)
    prep = (x >= PREPARATION[0]) & (x < PREPARATION[1])
    read = (x >= READOUT[0]) & (x < READOUT[1])
    kinetic = diags([-J*np.ones(SIDE-1), 2*J*np.ones(SIDE), -J*np.ones(SIDE-1)], [-1, 0, 1], format='csr')
    h = float(mass)/SCALE
    pieces = {
        'finite_controller_kinetic': kron(kinetic, eye(16), format='csr'),
        'native_source': kron(eye(SIDE), csr_matrix(h*n), format='csr'),
        'preparation_reference': kron(diags((x < PREPARATION[1])*h), csr_matrix(r1), format='csr'),
        'readout_reference': kron(diags((x >= READOUT[0])*h), csr_matrix(r2), format='csr'),
        'preparation_exchange': kron(diags(prep*coupling), csr_matrix(b1), format='csr'),
        'terminal_exchange': kron(diags(read*2.), csr_matrix(b2), format='csr')}
    H = sum(pieces.values())
    assert (H-H.conj().T).nnz == 0
    clock = np.exp(-((x-20)/5)**2/2)*(1j)**x
    clock[(x < 4) | (x > 36)] = 0
    clock /= np.linalg.norm(clock)
    if twirl:
        seeds = np.zeros((16, 4), complex)
        seeds[:4] = np.eye(4)/2
    else:
        seeds = np.zeros((16, 1), complex)
        seeds[:4, 0] = .5
    initial = np.kron(clock[:, None], seeds)
    return H, pieces, initial, aa, b1, b2, prep, read


def run(mass, coupling=2., *, twirl=False):
    H, pieces, initial, A, B1, B2, prep, read = source(mass, coupling, twirl=twirl)
    out = expm_multiply(-1j*H, initial, start=0., stop=TIMES[-1], num=len(TIMES), traceA=-1j*H.diagonal().sum())
    v = out.reshape(len(TIMES), SIDE, 16, initial.shape[1])
    density = np.sum(abs(v)**2, axis=(2, 3))
    norm_error = float(np.max(abs(density.sum(axis=1)-1)))
    assert norm_error < 1e-12
    coherence = np.einsum('txia,ij,txja->t', v.conj(), A, v)
    h = float(mass)/SCALE
    defect1 = B1@A-A@B1
    defect2 = B2@A-A@B2
    defect = (1j*coupling*np.einsum('txia,ij,txja->t', v[:, prep].conj(), defect1, v[:, prep])+
              2j*np.einsum('txia,ij,txja->t', v[:, read].conj(), defect2, v[:, read]))
    bound = abs(coupling)*density[:, prep].sum(axis=1)+2*density[:, read].sum(axis=1)
    assert np.max(abs(defect)-bound) < 1e-13
    free = (TIMES >= 2.5) & (TIMES <= 4.)
    phases = coherence*np.exp(-1j*h*TIMES)
    drift = float(max(abs(phases[free]-phases[free][0])))
    if not twirl:
        assert min(abs(coherence[free])) > .4
        assert drift < 1e-6
    else:
        assert max(abs(coherence)) < 1e-12
    final = out[-1]
    probability = np.sum(abs(final.reshape(SIDE, 16, initial.shape[1]))**2, axis=2)
    energy = []
    for name, op in pieces.items():
        first = float(np.vdot(initial, op@initial).real)*SCALE
        last = float(np.vdot(final, op@final).real)*SCALE
        after_readout = float(np.dot(op.diagonal().real, probability.ravel()))*SCALE
        energy.append({'component': name, 'initial_native_energy_diagnostic': first,
            'final_native_energy_diagnostic': last,
            'after_complete_occupation_instrument_native_energy_diagnostic': after_readout})
    balance = sum(x['final_native_energy_diagnostic']-x['initial_native_energy_diagnostic'] for x in energy)
    assert abs(balance) < 1e-7
    work = sum(x['after_complete_occupation_instrument_native_energy_diagnostic']-x['final_native_energy_diagnostic'] for x in energy)
    return {'mass_native_units_diagnostic': float(mass), 'preparation_coupling_over_controller_scale': coupling,
        'both_reference_phases_twirled': twirl,
        'all_native_reference_controller_probability_rows': probability.tolist(),
        'probability_row_axes': ['256 controller positions', '4 source occupations times 2 preparation-reference occupations times 2 readout-reference occupations'],
        'complete_record_total_diagnostic': float(probability.sum()),
        'readout_reference_excited_probability_diagnostic': float(probability[:, 1::2].sum()),
        'clock_not_beyond_readout_probability_diagnostic': float(probability[:READOUT[1]].sum()),
        'all_failures_reflections_and_inactive_occupations_retained': True,
        'orthogonal_native_Fock_complement_record_probability_exact': '0',
        'free_interval_native_times': [2.5/SCALE, 4./SCALE],
        'minimum_free_source_coherence_diagnostic': float(min(abs(coherence[free]))),
        'native_frequency_demodulated_drift_diagnostic': drift,
        'maximum_local_clock_contact_defect_diagnostic': float(max(abs(defect[free])))*SCALE,
        'maximum_full_carrier_contact_bound_diagnostic': float(max(bound[free]))*SCALE,
        'source_frequency_control_samples': [{'native_time': float(t)/SCALE,
            'coherence_real': float(z.real), 'coherence_imaginary': float(z.imag),
            'contact_defect_native_units': float(abs(d))*SCALE}
            for t, z, d in zip(TIMES[free][::5], coherence[free][::5], defect[free][::5])],
        'energy_components': energy,
        'closed_hamiltonian_total_energy_change_diagnostic': balance,
        'terminal_complete_occupation_measurement_work_diagnostic': work,
        'maximum_norm_drift_diagnostic': norm_error,
        'numerical_trajectory_is_an_outward_interval_certificate': False}


def certificate(masses):
    cases = [run(float(F(row['source_rest_energy_electron_units_exact'])), strength)
             for row in masses for strength in (2., 3.)]
    twirl = run(1., 2., twirl=True)
    return {'executed_reducing_carrier_dimension': 4096, 'time_independent_autonomous_hamiltonian': True,
        'clock_scale_native_units_exact': str(SCALE), 'controller_hopping_over_scale_exact': str(J),
        'clock_initial_packet': 'Normalized exp(-(x-20)^2/50)*i^x on positions4 through36; zero elsewhere',
        'architecture': 'Open 256-position clock chain, complete native four-occupation pair and two two-level references. Local exchange regions [48,64) and [192,208); reference gaps depend on clock position. All terms are fixed before propagation.',
        'semiboundedness': 'With L=4096 and recorded g=g_prepare/L, the positive controller kinetic and reference/source energies and disjoint exchange regions give H_total >= -L*max(abs(g),2), equivalently H_total/L >= -max(abs(g),2). An additive constant makes the Hamiltonian positive without changing records.',
        'native_complement': 'R=(P_vac+P_j) tensor I_memory reduces the native H and A. The 4096 displayed occupation effects sum to the executed R carrier; an explicit R-perp outcome completes the native identity and has probability zero on this input. No full native Fock enumeration is asserted.',
        'resource_scope': 'The moving finite clock packet and both phaseful references are supplied. Reference-gap switching is controlled by the clock position under one static Hamiltonian; its kinetic energy balance is retained.',
        'exact_local_identity': 'd<A>/dt=i*m*<A>+i*<g_prepare Pi_prepare [B1,A]+g_readout Pi_readout [B2,A]>',
        'exact_full_carrier_remainder': '|d<A>/dt-i*m*<A>| <= |g_prepare| Pr(preparation region)+|g_readout| Pr(readout region)',
        'native_charged_mass_controls': cases, 'both_phase_references_twirled_control': twirl,
        'terminal_detector': 'Complete controller/source/reference occupation instrument; its noncommuting clock readout has displayed mean work, not a claimed energy-compatible zero-resource implementation',
        'lifetime_scope': 'Finite flight only. Reflection, packet dispersion and end-boundary recurrence are retained; no asymptotic linewidth or sustained limit cycle is certified.',
        'literal_openwave_row3_validated': False}
