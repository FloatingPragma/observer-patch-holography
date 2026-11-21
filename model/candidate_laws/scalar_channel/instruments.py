"""Finite scalar preparation and source-on fuel/pointer instruments.

The numerical evolution follows the imported integer-Pauli compiler's even
parity representation. Local ramps, pulse actions and phase settings are
supplied services. Eigensystems are used for simulation and error auditing,
never to initialize the source state. Normal replay writes no artifacts.
"""

import math

import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
from scipy.sparse.linalg import expm_multiply

from .carrier import load_source, chi, require
from .compiler import build, hopping_words

TOL = 2e-11


def _expectation(operator, covariance):
    """Tr(operator covariance), without a cubic matrix product."""
    value = np.einsum('ij,ji->', operator, covariance)
    require(abs(value.imag) <= TOL, 'source energy acquired an imaginary part')
    return float(value.real)


def _energy_stage(name, before, after, work, **services):
    change = after - before
    residual = abs(change - work)
    require(residual <= TOL, 'source energy/work balance failed at ' + name)
    return {'stage': name, 'source_energy_before': before, 'source_energy_after': after,
            'source_energy_change': change, 'external_driving_work': work,
            'work_balance_residual': residual, **services}


def canonical(word):
    p, x, z = word
    q = (x & z).bit_count() % 4
    return (q, x, z), 1j ** ((p - q) % 4)


def compile_matrix(graph, matrix):
    """Imported number/hopping compiler with complete-domain input guards."""
    matrix = np.asarray(matrix, complex)
    n = graph['n']
    require(matrix.shape == (n, n) and np.isfinite(matrix).all(), 'invalid one-body matrix')
    require(np.max(abs(matrix - matrix.conj().T)) <= 1e-13, 'one-body matrix must be Hermitian')
    support = np.eye(n, dtype=bool)
    for i, j in graph['links']:
        support[i, j] = support[j, i] = True
    require(not np.any(matrix[~support] != 0), 'instrument outside fixed local graph')
    terms = {}

    def add(word, coefficient):
        word, phase = canonical(word)
        coefficient = complex(coefficient) * phase
        require(abs(coefficient.imag) <= 1e-13, 'non-Hermitian compiled term')
        terms[word] = terms.get(word, 0.0) + float(coefficient.real)

    for i in range(n):
        if matrix[i, i]:
            add((0, 0, 0), matrix[i, i] / 2)
            add(graph['B'][i], -matrix[i, i] / 2)
    for k, (i, j) in enumerate(graph['links']):
        if matrix[i, j]:
            for coefficient, word in hopping_words(graph, k, complex(matrix[i, j])):
                add(word, coefficient)
    return {word: coefficient for word, coefficient in terms.items() if coefficient}


def decode_paulis(graph, terms):
    """Imported coefficient decoder; its input contains no one-body matrix."""
    require(type(terms) is dict and all(np.isscalar(v) and np.isfinite(v) and np.isreal(v)
                                      for v in terms.values()), 'invalid Pauli coefficients')
    for word in terms:
        require(isinstance(word, tuple) and len(word) == 3 and all(type(v) is int for v in word)
                and 0 <= word[0] < 4 and 0 <= word[1] < 1 << graph['M']
                and 0 <= word[2] < 1 << graph['M']
                and word[0] == (word[1] & word[2]).bit_count() % 4, 'noncanonical Pauli word')
    matrix = np.zeros((graph['n'], graph['n']), complex)
    reconstructed = {}

    def add(word, coefficient):
        reconstructed[word] = reconstructed.get(word, 0.0) + coefficient

    for i, number in enumerate(graph['B']):
        word, phase = canonical(number)
        value = -2 * terms.get(word, 0.0) / phase
        matrix[i, i] = value
        add((0, 0, 0), float(value.real) / 2)
        add(word, -float(value.real) / 2)
    for k, (i, j) in enumerate(graph['links']):
        components = []
        for coefficient in (1 + 0j, 1j):
            basis = hopping_words(graph, k, coefficient)
            norm = sum(value * value for value, word in basis)
            value = sum(value * terms.get(word, 0.0) for value, word in basis) / norm
            for scale, word in basis:
                add(word, value * scale)
            components.append(value)
        matrix[i, j] = components[0] + 1j * components[1]
        matrix[j, i] = matrix[i, j].conjugate()
    residual = max((abs(terms.get(word, 0) - reconstructed.get(word, 0))
                    for word in set(terms) | set(reconstructed)), default=0.0)
    return matrix, float(residual)


def source_matrix(graph, variant='native'):
    require(variant in ('native', 'no_hopping', 'equal_diagonals'), 'undeclared source control')
    matrix = np.zeros((graph['n'], graph['n']), complex)
    diagonal = [1] * graph['E'] + ([-1] * graph['W'] if variant != 'equal_diagonals' else [1] * graph['W'])
    matrix[np.arange(graph['nsource']), np.arange(graph['nsource'])] = diagonal
    if variant != 'no_hopping':
        for (i, j), value in graph['source_coeff'].items():
            matrix[i, j] = matrix[j, i] = value
    return matrix


def prepare_source(steps=16):
    """Execute twelve midpoint star ramps, each of total source action 16.

    256 steps per star are the preparation of the construction manual. The
    shorter default is a different declared preparation with its own audit.
    """
    require(type(steps) is int and 1 <= steps <= 4096, 'steps must be an integer in [1,4096]')
    data = load_source()
    graph = build()
    native = source_matrix(graph)
    decoded, residual = decode_paulis(graph, compile_matrix(graph, native))
    require(residual <= TOL and la.norm(decoded - native, 'fro') <= TOL,
            'preparation source compiler residual')
    E, W = data['E'], data['W']
    incidence = sp.csr_matrix(decoded[:E, E:E + W])
    frame = np.vstack((np.zeros((E, W)), np.eye(W))).astype(complex)
    base = sp.diags(np.r_[np.ones(E), -np.ones(W)], format='csr')
    duration = 16 / steps
    receipts = []
    for center in range(12):
        mask = np.array([int(w[1]) == center for w in data['wedges']], float)
        block = incidence @ sp.diags(mask)
        star = sp.bmat([[sp.csr_matrix((E, E)), block], [block.T, sp.csr_matrix((W, W))]], format='csr')
        for k in range(steps):
            s = (k + .5) / steps
            generator = base + (3 * s * s - 2 * s * s * s) * star
            frame = expm_multiply(-1j * duration * generator, frame,
                                 traceA=-1j * duration * generator.diagonal().sum())
        base += star
        receipts.append({'center': center, 'substeps': steps, 'source_action': 16.0})
    covariance = frame @ frame.conj().T
    isometry_error = float(la.norm(frame.conj().T @ frame - np.eye(W), 2))
    require(isometry_error <= TOL, 'finite preparation lost occupied-frame normalization')
    # This eigensystem only audits the state produced above.
    values, vectors = la.eigh(data['D'])
    ground = vectors[:, values < 0] @ vectors[:, values < 0].conj().T
    error = float(la.norm(covariance - ground, 2))
    circulation_error = max(float(la.norm(covariance @ chi(data, b) - chi(data, b))) for b in range(12))
    require(circulation_error <= TOL, 'finite preparation lost a local circulation')
    if steps == 256:
        require(error < .0058, 'the 256-step preparation exceeds its allowance')
    return frame, {
        'steps_per_star': steps, 'total_source_action_per_star': 16.0,
        'source_action_per_step': duration, 'midpoint_profile': '3*s**2-2*s**3',
        'initially_filled_transit_modes': W, 'initially_empty_seam_modes': E,
        'frame_isometry_error': isometry_error, 'distance_from_exact_sea': error,
        'energy_above_exact_sea': float(np.trace(data['D'] @ (covariance - ground)).real),
        'largest_circulation_occupation_error': circulation_error,
        'manual_preparation_reproduced': steps == 256,
        'exact_sea_eigenvectors_used_to_initialize': False,
        'local_ramp_receipts': receipts,
    }


def pointer_permutation(fuel, word, read_a=True, read_b=True):
    """Two controlled-X writes on the complete dirty two-pointer domain."""
    require(type(fuel) is int and fuel in (0, 1), 'fuel must be a bit')
    require(type(word) is int and 0 <= word < 4, 'pointer word outside two-bit register')
    require(type(read_a) is bool and type(read_b) is bool, 'pointer controls must be booleans')
    return word ^ ((int(read_a) + 2 * int(read_b)) * fuel)


def pointer_law(fuel_probability, word=0, read_a=True, read_b=True):
    require(np.isscalar(fuel_probability) and np.isfinite(fuel_probability)
            and 0 <= fuel_probability <= 1, 'invalid fuel probability')
    result = np.zeros(4)
    for fuel in (0, 1):
        result[pointer_permutation(fuel, word, read_a, read_b)] += fuel_probability if fuel else 1 - fuel_probability
    return result.tolist()


def native_family(graph, source, covariance, target, queries=(0, 4, 8), strength=64.0):
    """Pure even-particle source, both paired apparatus inputs, all outcomes.

    The source covariance must be a projector of even rank. Fuel and parity
    reference start together in 00 or 11, so they add zero or two particles
    to the implemented even-total-parity code. A general covariance alone
    does not certify that preparation or its parity sector.

    The work ledger concerns dGamma(source), with zero bare fuel/reference
    energy. Each pulse switches V on, evolves with source + V, and switches
    V off. Its supplied driving work is <V>before - <V>after. The source is
    on throughout that evolution. This is no controller-energy model.
    """
    ns, fuel, reference = graph['nsource'], graph['F'][0], graph['R'][0]
    source = np.asarray(source, complex)
    require(source.shape == (graph['n'], graph['n']) and np.isfinite(source).all()
            and la.norm(source - source.conj().T, 'fro') <= TOL, 'invalid free source generator')
    require(0 < ns <= fuel < graph['n'] and ns <= reference < graph['n']
            and fuel != reference, 'fuel/reference must be distinct apparatus modes')
    require(not np.any(source[ns:, :] != 0) and not np.any(source[:, ns:] != 0),
            'free generator must act only on source modes')
    covariance = np.asarray(covariance, complex)
    target = np.asarray(target, complex)
    require(covariance.shape == (ns, ns) and np.isfinite(covariance).all()
            and la.norm(covariance - covariance.conj().T, 'fro') <= TOL, 'invalid source covariance')
    eigenvalues = la.eigvalsh(covariance)
    require(eigenvalues[0] >= -TOL and eigenvalues[-1] <= 1 + TOL, 'source covariance outside fermion domain')
    projector_error = float(la.norm(covariance @ covariance - covariance, 'fro'))
    require(projector_error <= TOL, 'source covariance must be a pure Slater projector')
    trace = float(np.trace(covariance).real)
    source_particles = int(np.count_nonzero(eigenvalues > .5))
    trace_error = abs(trace - source_particles)
    # Rounding each Hermitian eigenvalue to {0,1} gives a projector P with
    # ||C-P||_F <= 2 ||C*C-C||_F, hence |Tr(C)-rank(P)| <= 2 sqrt(ns) TOL.
    trace_tolerance = 2 * math.sqrt(ns) * TOL
    require(trace_error <= trace_tolerance and source_particles % 2 == 0,
            'source Slater projector must have even particle number')
    require(target.shape == (ns,) and np.isfinite(target).all()
            and abs(np.vdot(target, target).real - 1) <= TOL, 'target must be a normalized local mode')
    require(type(strength) in (int, float) and np.isfinite(strength) and strength > 0, 'positive finite pulse strength required')
    require(len(queries) > 0 and all(type(q) is int and q >= 0 for q in queries), 'nonnegative local queries required')
    pulse = np.array(source, dtype=complex, copy=True)
    pulse[:ns, fuel] = strength * target
    pulse[fuel, :ns] = strength * target.conj()
    words = compile_matrix(graph, pulse)
    decoded, residual = decode_paulis(graph, words)
    matrix_error = float(la.norm(decoded - pulse, 'fro'))
    require(residual <= TOL and matrix_error <= TOL, 'compiled pulse decoder residual')
    pulse_duration = math.pi / (4 * strength)
    U = la.expm(-1j * pulse_duration * decoded)
    unitary_error = float(la.norm(U.conj().T @ U - np.eye(graph['n']), 'fro'))
    require(unitary_error <= TOL, 'finite pulse unitarity residual')
    converter = decoded - source
    source_after_pulse = U.conj().T @ source @ U
    converter_after_pulse = U.conj().T @ converter @ U
    pulse_services = {'duration': pulse_duration, 'source_strength': 1.0,
                      'source_action': pulse_duration, 'converter_strength': float(strength),
                      'converter_action': math.pi / 4}
    prepared = []
    for initial in (0, 1):
        initial_covariance = np.zeros_like(source)
        initial_covariance[:ns, :ns] = covariance
        initial_covariance[fuel, fuel] = initial_covariance[reference, reference] = initial
        after = U @ initial_covariance @ U.conj().T
        before_energy = _expectation(source, initial_covariance)
        after_energy = _expectation(source, after)
        switch_on = _expectation(converter, initial_covariance)
        switch_off = -_expectation(converter, after)
        stage = _energy_stage('pulse_1', before_energy, after_energy, switch_on + switch_off,
                              switch_on_work=switch_on, switch_off_work=switch_off, **pulse_services)
        prepared.append((after, stage))
    values, vectors = la.eigh(source)
    rows, reference_error, roundoff, work_error = [], 0.0, 0.0, 0.0
    for q in queries:
        free_duration = q / 8
        free = (vectors * np.exp(-1j * free_duration * values)) @ vectors.conj().T
        after_free = []
        for after_first, first_stage in prepared:
            evolved = free @ after_first @ free.conj().T
            energy = _expectation(source, evolved)
            stage = _energy_stage('free_interval', first_stage['source_energy_after'], energy, 0.0,
                                  duration=free_duration, source_strength=1.0,
                                  source_action=free_duration, converter_strength=0.0,
                                  converter_action=0.0)
            after_free.append((evolved, stage))
        for phase in range(4):
            setting = np.ones(graph['n'], complex)
            setting[fuel] = np.exp(-1j * math.pi * phase / 2)
            fuel_row = (U[fuel, :] * setting) @ free @ U
            reference_row = (U[reference, :] * setting) @ free @ U
            for initial in (0, 1):
                raw = float((fuel_row[:ns] @ covariance @ fuel_row[:ns].conj()).real
                            + initial * (abs(fuel_row[fuel]) ** 2 + abs(fuel_row[reference]) ** 2))
                require(-TOL <= raw <= 1 + TOL, 'fuel law outside probability domain')
                probability = min(1.0, max(0.0, raw))
                roundoff = max(roundoff, abs(probability - raw))
                ref = float((reference_row[:ns] @ covariance @ reference_row[:ns].conj()).real
                            + initial * (abs(reference_row[fuel]) ** 2 + abs(reference_row[reference]) ** 2))
                reference_error = max(reference_error, abs(ref - initial))
                evolved, free_stage = after_free[initial]
                phased = setting[:, None] * evolved * setting.conj()[None, :]
                phase_energy = _expectation(source, phased)
                phase_stage = _energy_stage('fuel_phase', free_stage['source_energy_after'], phase_energy, 0.0,
                                            duration=None, source_strength=None, source_action=None,
                                            converter_strength=0.0, converter_action=0.0,
                                            fuel_number_action=math.pi * phase / 2,
                                            supplied_unitary='exp(-i*phase*pi*n_fuel/2)')
                final_energy = _expectation(source_after_pulse, phased)
                switch_on = _expectation(converter, phased)
                switch_off = -_expectation(converter_after_pulse, phased)
                second_stage = _energy_stage('pulse_2', phase_energy, final_energy, switch_on + switch_off,
                                             switch_on_work=switch_on, switch_off_work=switch_off,
                                             **pulse_services)
                # Ignoring the two pointers dephases the fuel covariance.
                # The guarded source block has no fuel entries, so the
                # adjoint pointer channel fixes this observable exactly.
                pointer_stage = _energy_stage('terminal_pointer', final_energy, final_energy, 0.0,
                                              duration=None, source_strength=None, source_action=None,
                                              converter_strength=0.0, converter_action=0.0,
                                              supplied_unitary='CNOT(fuel,A) CNOT(fuel,B)',
                                              source_observable_fixed_by_adjoint_channel=True)
                stages = [prepared[initial][1], free_stage, phase_stage, second_stage, pointer_stage]
                total_work = sum(stage['external_driving_work'] for stage in stages)
                energy_change = final_energy - stages[0]['source_energy_before']
                balance = abs(energy_change - total_work)
                require(balance <= TOL, 'complete source energy/work balance failed')
                work_error = max(work_error, balance, *(stage['work_balance_residual'] for stage in stages))
                rows.append({'query': q, 'phase': phase, 'initial_fuel': initial,
                             'initial_reference': initial, 'raw_final_fuel_probability': raw,
                             'pointer_probabilities': pointer_law(probability),
                             'source_energy_ledger': {
                                 'observable': 'dGamma(free_source)', 'stages': stages,
                                 'net_source_energy_change': energy_change,
                                 'total_external_driving_work': total_work,
                                 'work_balance_residual': balance}})
    require(reference_error <= TOL, 'decoupled parity reference changed')
    return rows, {'pauli_terms': len(words), 'decode_residual': residual,
                  'source_preparation_domain': 'pure even-particle Slater source; paired fuel/reference input',
                  'source_particle_number': source_particles,
                  'source_particle_number_error': trace_error,
                  'source_particle_number_tolerance': trace_tolerance,
                  'source_projector_error': projector_error,
                  'decoded_matrix_error': matrix_error, 'pulse_unitarity_error': unitary_error,
                  'matrix_error_norm': 'Frobenius upper bound on operator norm',
                  'parity_reference_error': reference_error, 'largest_probability_roundoff': roundoff,
                  'largest_source_work_balance_residual': work_error,
                  'source_commutes_with_fuel_phase_and_pointer': True}


def run(extended=False):
    require(type(extended) is bool, 'extended must be a boolean')
    graph, data = build(), load_source()
    frame, preparation = prepare_source(256 if extended else 16)
    covariance = frame @ frame.conj().T
    ns = graph['nsource']
    single = np.zeros(ns); single[0] = 1
    star = np.zeros(ns)
    star[[i for i, edge in enumerate(data['edges']) if 0 in edge]] = 1 / math.sqrt(5)
    a, b, _ = map(int, data['wedges'][0])
    seed = np.zeros(ns, complex)
    seed[data['index'][tuple(sorted((a, b)))]] = 1 / math.sqrt(2)
    seed[data['E']] = 1j / math.sqrt(2)
    require(b == 0, 'complex seed outside active fuel star')
    targets = {'single_seam': single, 'circulation': chi(data), 'five_seam': star, 'complex_seed': seed}
    cases, by_key = [], {}
    for variant in ('native', 'no_hopping', 'equal_diagonals'):
        source = source_matrix(graph, variant)
        decoded, residual = decode_paulis(graph, compile_matrix(graph, source))
        require(residual <= TOL and la.norm(decoded - source, 'fro') <= TOL, 'compiled free source residual')
        for name, target in targets.items():
            records, checks = native_family(graph, decoded, covariance, target)
            cases.append({'variant': variant, 'target': name, 'records': records, 'checks': checks})
            by_key[variant, name] = np.array([row['pointer_probabilities'] for row in records])
    controls = {}
    for variant in ('no_hopping', 'equal_diagonals'):
        difference = max(float(np.max(abs(by_key['native', name] - by_key[variant, name]))) for name in targets)
        require(difference > .01, 'same-apparatus source control has no resolved effect')
        controls[variant] = {'largest_pointer_probability_difference': difference,
                             'same_preparation_targets_pulses_and_readouts': True}
    return {'schema': 'oph.scalar_channel.finite_instruments.v2', 'extended': extended,
            'preparation': preparation,
            'queries': [0, 4, 8], 'source_action_per_local_reference_count': 1 / 8,
            'pulse_strength': 64.0, 'pulse_action': math.pi / 4,
            'pointer_outcome_order': ['00', '10', '01', '11'],
            'cases': cases, 'controls': controls,
            'source_on_during_pulses': True, 'complete_binary_fuel_inputs_retained': True,
            'all_four_pointer_outcomes_retained': True,
            'parity_reference_is_decoupled_and_unread': True,
            'preparation_and_phase_services_are_autonomous': False,
            'source_work_accounting': {
                'observable': 'dGamma(free_source)',
                'bare_fuel_and_reference_energies': 0.0,
                'pulse_work': 'external switch-on <V>before plus switch-off -<V>after',
                'pulse_generator': 'free_source + converter',
                'pulse_source_strength': 1.0, 'pulse_duration': math.pi / 256,
                'phase_and_pointer_durations_declared': False,
                'phase_and_pointer_source_work': 0.0,
                'controller_energy_modeled': False,
                'autonomous_energy_conservation_claimed': False,
                'includes_preparation_ramp_work': False}}
