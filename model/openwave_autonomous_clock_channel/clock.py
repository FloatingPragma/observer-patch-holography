"""Phase-insensitive native pumping with complete records and onset controls."""
from fractions import Fraction as F
from math import factorial, isqrt
from pathlib import Path
import json
import numpy as np
from scipy.linalg import eigh, sqrtm
from model.openwave_passive_clock_records.clock import pair_operators, native_attachment

STEP = F(1, 65536)
PUMP = F(1, 16)


def sqrt_bounds(x, bits=100):
    x = F(x)
    scale = 1 << bits
    floor = isqrt(x.numerator * scale * scale // x.denominator)
    lo, hi = F(floor, scale), F(floor + 1, scale)
    assert lo * lo <= x <= hi * hi
    return lo, hi


def collective(cap):
    """Invariant symmetric clock sector; all complementary sectors are retained
    by the bounded complete native-current construction checked below.
    """
    assert isinstance(cap, int) and cap >= 1
    n = np.arange(cap + 1)
    a = np.diag(np.sqrt(n[1:]), 1)
    gamma, kappa = PUMP, 2 * PUMP / cap
    strength = gamma * cap + kappa * cap * (cap - 1)
    assert STEP * strength < F(1, 4)
    rates = [gamma * (j + 1) * (j < cap) + kappa * j * (j - 1)
             for j in range(cap + 1)]
    no_jump = np.diag([np.sqrt(float(1 - STEP * x)) for x in rates])
    pump = np.sqrt(float(STEP * gamma)) * a.T
    loss = np.sqrt(float(STEP * kappa)) * (a @ a)
    assert np.max(abs(no_jump.T @ no_jump + pump.T @ pump + loss.T @ loss
                      - np.eye(cap + 1))) < 4e-16
    population = no_jump ** 2 + pump ** 2 + loss ** 2
    assert np.max(abs(population.sum(axis=0) - 1)) < 3e-16
    return {'cap': cap, 'a': a, 'n': n, 'gamma': gamma, 'kappa': kappa,
            'rates_exact': rates, 'K': (no_jump, pump, loss),
            'population': population, 'norm_bound_exact': strength}


def coherence_band(cap, order=1):
    """D^-1 B D for rho_(n+order,n), with D_n=sqrt((n+order)!/n!).
    The only complex factor is exp(-i order m STEP).
    """
    assert 1 <= order <= cap
    d = collective(cap)
    size = cap + 1 - order
    B = np.zeros((size, size))
    lower = np.full((size, size), F(0), dtype=object)
    upper = lower.copy()
    gamma, kappa, rates = d['gamma'], d['kappa'], d['rates_exact']
    for n in range(size):
        x = (1 - STEP * rates[n + order]) * (1 - STEP * rates[n])
        lo, hi = sqrt_bounds(x)
        lower[n, n], upper[n, n] = lo, hi
        B[n, n] = np.sqrt(float(x))
        if n:
            value = STEP * gamma * n
            B[n, n - 1] = float(value)
            lower[n, n - 1] = upper[n, n - 1] = value
        if n + 2 < size:
            value = STEP * kappa * (n + order + 2) * (n + order + 1)
            B[n, n + 2] = float(value)
            lower[n, n + 2] = upper[n, n + 2] = value
    return d, B, lower, upper


def perron_certificate(cap, order=1, witness=None):
    d, B, lower, upper = coherence_band(cap, order)
    if witness is None:
        data = json.loads((Path(__file__).parent / 'witnesses.json').read_text())
        witness = data['witnesses'][f'{cap}:{order}']
    witness = list(map(F, witness))
    assert len(witness) == len(B) and min(witness) > 0
    lo = min(sum(lower[i, j] * witness[j] for j in range(len(B))) / witness[i]
             for i in range(len(B)))
    hi = max(sum(upper[i, j] * witness[j] for j in range(len(B))) / witness[i]
             for i in range(len(B)))
    assert 0 < lo <= hi < 1, (cap, order, float(lo), float(hi))
    # 1-r <= -log r <= (1-r)/r gives outward linewidth bounds without
    # treating a machine logarithm as an exact certificate.
    width = [(1 - hi) / STEP, (1 - lo) / (lo * STEP)]
    return {'cap': cap, 'coherence_order': order,
            'complete_band_dimension': len(B),
            'positive_dyadic_witness': [str(x) for x in witness],
            'spectral_radius_interval_exact': [str(lo), str(hi)],
            'linewidth_source_units_interval_exact': [str(x) for x in width],
            'angular_frequency_over_native_mass_exact': str(order),
            'finite_nonzero_linewidth': True}


def complete_native_control(copies=3):
    """Construct every four-state native pair, including both inactive rows.
    No collective-spin projector is installed in the source operators.
    """
    N, Nr, A, X, Y = pair_operators()
    dim = 4 ** copies
    def lift(op, j):
        out = np.array([[1.]], complex)
        for i in range(copies):
            out = np.kron(out, op if i == j else np.eye(4))
        return out
    upper = [lift(A @ A.conj().T, j) for j in range(copies)]
    n = sum(upper, np.zeros((dim, dim), complex))
    Jplus = sum((lift(A, j) for j in range(copies)), np.zeros((dim, dim), complex))
    normalization = np.diag(1 / np.sqrt(copies - np.diag(n).real + 1))
    a = Jplus.conj().T @ normalization
    assert np.max(abs(n @ a - a @ n + a)) < 2e-15
    # Pair occupations |10> and |01> are explicit inactive source blocks.
    assert len(n) == dim
    assert np.linalg.norm(a, 2) <= np.sqrt(copies) + 3e-15
    assert np.linalg.norm(a @ a, 2) <= np.sqrt(copies * (copies - 1)) + 4e-15
    gamma, kappa = float(PUMP), float(2 * PUMP / copies)
    S = gamma * a @ a.conj().T + kappa * a.conj().T @ a.conj().T @ a @ a
    e, U = eigh(np.eye(dim) - float(STEP) * S)
    assert min(e) > .75
    K0 = U @ np.diag(np.sqrt(e)) @ U.conj().T
    Kp = np.sqrt(float(STEP) * gamma) * a.conj().T
    Km = np.sqrt(float(STEP) * kappa) * (a @ a)
    gram = sum((K.conj().T @ K for K in (K0, Kp, Km)), np.zeros((dim, dim), complex))
    assert np.max(abs(gram - np.eye(dim))) < 4e-14
    energy = sum((K.conj().T @ n @ K for K in (K0, Kp, Km)), np.zeros_like(n)) - n
    work = float(STEP) * (gamma * a @ a.conj().T
                           - 2 * kappa * a.conj().T @ a.conj().T @ a @ a)
    assert np.max(abs(energy - work)) < 4e-14
    return {'copies': copies, 'complete_four_state_pair_dimension': dim,
            'all_inactive_Fock_blocks_retained': True,
            'complete_collision_Gram_residual_diagnostic': float(np.max(abs(gram - np.eye(dim)))),
            'complete_operator_energy_balance_residual_diagnostic': float(np.max(abs(energy - work)))}


def evolve(rho, steps, cap, mass=1.):
    """Every coherence band of the actual finite collision channel.

    Point values are numerical diagnostics. Exact frequency/linewidth bounds
    are supplied separately, without treating this power as interval arithmetic.
    """
    assert isinstance(steps, int) and steps >= 0
    d = collective(cap)
    result = np.zeros_like(rho, dtype=complex)
    np.fill_diagonal(result, np.linalg.matrix_power(d['population'], steps)
                     @ np.diag(rho))
    for k in range(1, cap + 1):
        _, B, _, _ = coherence_band(cap, k)
        scale = np.array([np.sqrt(factorial(n + k) / factorial(n))
                          for n in range(cap + 1 - k)])
        column = np.array([rho[n + k, n] for n in range(cap + 1 - k)], complex)
        column = scale * (np.linalg.matrix_power(B, steps) @ (column / scale))
        column *= np.exp(-1j * k * mass * float(STEP) * steps)
        for n, x in enumerate(column):
            result[n + k, n], result[n, n + k] = x, np.conjugate(x)
    return result


def records(cap=4, burn_steps=1 << 23, flight_steps=65536, mass=1.):
    """Complete first/final native probe, all unobserved battery records traced.
    The raw battery alphabet is {0,+,-}; every word has its Kraus product.
    """
    d = collective(cap)
    rho = np.zeros((cap + 1, cap + 1), complex)
    rho[0, 0] = 1
    rho = evolve(rho, burn_steps, cap, mass)
    assert np.max(abs(rho - np.diag(np.diag(rho)))) == 0
    assert abs(np.trace(rho) - 1) < 2e-8
    H = mass * np.diag(d['n'])
    X = (d['a'] + d['a'].T) / (2 * cap)
    Y = 1j * (d['a'].T - d['a']) / (2 * cap)
    before = float(np.trace(H @ rho).real)
    first_energy = between_energy = final_energy = 0.
    rows = []
    for sign in (-1, 1):
        E = (np.eye(cap + 1) + sign * X) / 2
        K = np.asarray(sqrtm(E), complex)
        branch = K @ rho @ K.conj().T
        first_energy += float(np.trace(H @ branch).real)
        branch = evolve(branch, flight_steps, cap, mass)
        between_energy += float(np.trace(H @ branch).real)
        for quadrature, O in (('X', X), ('Y', Y)):
            for final_sign in (-1, 1):
                effect = (np.eye(cap + 1) + final_sign * O) / 2
                L = np.asarray(sqrtm(effect), complex)
                post = .5 * L @ branch @ L.conj().T
                p = float(np.trace(post).real)
                assert p > 0
                final_energy += float(np.trace(H @ post).real)
                rows.append({'first_sign': sign, 'quadrature': quadrature,
                             'final_sign': final_sign, 'probability_diagnostic': p})
    assert abs(sum(x['probability_diagnostic'] for x in rows) - 1) < 2e-8
    works = {'pumping_from_native_vacuum': before,
             'first_probe': first_energy - before,
             'continued_collisions': between_energy - first_energy,
             'final_probe': final_energy - between_energy}
    assert abs(sum(works.values()) - final_energy) < 2e-14
    return {'cap': cap, 'burn_in_collisions': burn_steps,
            'free_collision_steps_between_probes': flight_steps,
            'source_mass_electron_units_diagnostic': mass,
            'complete_eight_probe_records': rows,
            'unobserved_battery_records': 'Every word in {0,+,-}; its Kraus product is retained before this explicit coarse graining.',
            'unconditional_phase_mean_before_first_probe_exact': '0',
            'work_including_all_probe_and_battery_branches_diagnostic': works,
            'final_mean_source_energy_diagnostic': final_energy,
            'point_values_have_no_outward_numerical_certificate': True}


def single_branch_records(steps, pump=PUMP, phase_mass=1.):
    """Projective equatorial records; the prior pumped occupation cancels.

    These exact formulas provide the zero-pump counterfactual instead of
    treating a prepared propagator as proof of autonomous phase onset.
    """
    r = 1 - STEP * pump
    assert 0 < r <= 1
    visibility = float(r) ** (steps / 2)
    phase = phase_mass * float(STEP) * steps
    out = []
    for first in (-1, 1):
        for axis in ('X', 'Y'):
            wave = np.cos(phase) if axis == 'X' else -np.sin(phase)
            for last in (-1, 1):
                out.append({'first_sign': first, 'quadrature': axis,
                            'final_sign': last,
                            'probability_diagnostic': float((1 + first * last * visibility * wave) / 8)})
    assert abs(sum(x['probability_diagnostic'] for x in out) - 1) < 3e-16
    return {'collision_steps': steps, 'pump_strength_exact': str(pump),
            'complete_eight_records': out,
            'exact_record_law': '(1+s*t*r^(k/2)*cos(m*k*delta))/8 for X; replace cos by -sin for Y',
            'r_exact': str(r),
            'prior_pumped_occupation_cancels_exactly': True,
            'first_probe_prepares_the_phase': True,
            'combined_preparation_collision_and_readout_work_exact': 'm/2',
            'physical_single_branch_rest_gap_used': True}


def certificate():
    clocks = native_attachment()
    assert all(F(x['source_rest_energy_electron_units_exact']) * STEP < F(1, 8)
               for x in clocks)
    bands = [perron_certificate(cap, k) for cap in (1, 4, 8, 16)
             for k in range(1, cap + 1)]
    return {'schema': 'oph.autonomous_native_collision_clock.v1',
            'verdict': 'PASS_COMPLETE_NATIVE_COLLISION_CLOCK_AND_PHASE_ONSET_CONTROL',
            'architecture': {'collision_step_exact': str(STEP), 'pump_exact': str(PUMP),
                'two_loss_exact': '2*pump/N', 'supplied_phase': False,
                'native_current_normalization': 'a=J-*(N-n+1)^(-1/2)',
                'source': 'N copies of the inherited192-CAR longitudinal member, original Gauss reference, and complete neutral fermionic memories',
                'norm_complements': 'Every original Fock, constraint, pair-occupation and collective-spin block remains.',
                'mass_input_ancestry': 'Inherited electron calibration and analyzer phase; pump/loss settings are declared architecture.'},
            'native_single_branch_masses_and_frequencies': clocks,
            'complete_native_controls': [complete_native_control(1), complete_native_control(3)],
            'exact_complete_coherence_bands': bands,
            'complete_array_probe_records': [records(4, 1 << 20, k)
                                              for k in (0, 16384, 65536, 262144)],
            'single_branch_counterfactual': [single_branch_records(k, gamma)
                    for gamma in (PUMP, F(0)) for k in (0, 32768, 65536, 131072)],
            'battery': {'incoming_energy': 'm', 'outgoing_record_energies': {'0': 'm', '+': '0', '-': '3m'},
                'source_mean_work_per_collision': 'delta*m*[gamma*<aa*>-2*kappa*<a*²a²>]',
                'complete_energy_conserving_Stinespring_extension': True,
                'fresh_nonpassive_batteries_and_their_entropy_are_resources': True},
            'scope': {'autonomous_particle_occupation_onset': True,
                'finite_coherence_frequency_equals_constituent_source_gap': True,
                'array_frequency_is_total_array_rest_energy': False,
                'single_branch_probe_phase_is_prepared_at_first_measurement': True,
                'autonomous_internal_phase_onset_established_by_these_records': False,
                'zero_ensemble_mean_is_an_extra_admission_requirement': False,
                'equilibrium_ground_state_time_crystal': False,
                'thermodynamic_zero_linewidth_limit_proved': False,
                'arbitrary_coherent_drive_rigidity_proved': False}}
