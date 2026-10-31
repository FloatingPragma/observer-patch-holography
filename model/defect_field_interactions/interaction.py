"""A channel-correct external Higgs drive on the spatial Fock1080 extension.

Every active block is a complete two-mode CAR space, not a restriction to
one or two particles. The common phase drive is specified with the gauge
background held fixed; this is not a passive gauge transformation.
"""
from __future__ import annotations

from functools import lru_cache
from decimal import Decimal, localcontext
import math

import numpy as np

from ..common_source import calibrated_matter as matter
from ..common_source.chiral_candidate import action
from ..defect_field_benchmark.quasiparticles import particle_hole, quasiparticles, validate_chart
from ..observer.source_scalar import compiler

PORTS = 12
LOCAL_MODES = 90
SOURCE_GUARD = 2e-7
OPERATOR_GUARD = 2e-10


def _numeric(value):
    return np.asarray(value.evalf(17) if hasattr(value, 'evalf') else value, complex)


@lru_cache(maxsize=1)
def channel_pairings():
    """Actual all-left invariant coefficients, including conjugate-Higgs flags."""
    coefficients = matter.higgs_components()
    channels = action._representation()[4]
    epsilon = _numeric(action.SPIN_EPSILON)
    norm = float(matter.neutral_higgs_norm())
    return {name: {'pairing': np.kron(epsilon, (norm*_numeric(rows[1])).conj()),
                   'creation_phase_sign': 1 if channels[name]['conjugate_Higgs'] else -1,
                   'conjugate_Higgs': channels[name]['conjugate_Higgs']}
            for name, rows in coefficients.items()}


def higgs_pairing(phase):
    if isinstance(phase, (bool, str)) or not math.isfinite(float(phase)):
        raise ValueError('a finite external Higgs phase required')
    return sum((np.exp(1j*row['creation_phase_sign']*phase)*row['pairing']
                for row in channel_pairings().values()), np.zeros((90, 90), complex))


def nambu(pairing):
    pairing = np.asarray(pairing)
    if pairing.shape != (90, 90):
        raise ValueError('the complete local90 pairing matrix required')
    zero = np.zeros_like(pairing)
    return np.block([[zero, pairing], [pairing.conj().T, zero]])


@lru_cache(maxsize=1)
def paired_chart():
    """Rotate only degenerate energy/charge spaces into42 neutral CAR pairs."""
    chart = quasiparticles()
    base = chart['bdg']
    phase_pairing = sum((1j*row['creation_phase_sign']*row['pairing']
                         for row in channel_pairings().values()), np.zeros((90, 90), complex))
    derivative = nambu(phase_pairing)
    if np.linalg.norm(nambu(higgs_pairing(0))-base) > OPERATOR_GUARD:
        raise ArithmeticError('the external drive does not begin at the same action')
    if np.linalg.norm(base@derivative+derivative@base) > OPERATOR_GUARD:
        raise ArithmeticError('channel-correct phase derivative lost mass anticommutation')
    w = chart['positive_basis']
    anomalous = w.conj().T@derivative@particle_hole(w)
    if np.linalg.norm(w.conj().T@derivative@w) > OPERATOR_GUARD:
        raise ArithmeticError('drive gained a number-preserving derivative')
    rotation = np.eye(90, dtype=complex)
    pairs = []
    for block in chart['projectors']:
        if block['charge6'] >= 0:
            continue
        opposite = next(other for other in chart['projectors']
                        if other['charge6'] == -block['charge6']
                        and abs(other['energy']-block['energy']) < 2e-9)
        left, right = block['indices'], opposite['indices']
        unitary = anomalous[np.ix_(left, right)]/block['energy']
        rotation[np.ix_(right, right)] = unitary.T
        pairs.extend({'left': a, 'right': b, 'energy': block['energy'],
                      'left_charge6': block['charge6']} for a, b in zip(left, right, strict=True))
    basis = w@rotation
    complete = validate_chart(basis, chart['energies'], chart['charge6'])
    desired = np.zeros((90, 90), complex)
    for pair in pairs:
        a, b, energy = pair['left'], pair['right'], pair['energy']
        desired[a, b], desired[b, a] = energy, -energy
    transformed = complete.conj().T@derivative@complete
    target = np.block([[np.zeros((90, 90)), desired], [desired.conj().T, np.zeros((90, 90))]])
    if (len(pairs) != 42 or np.linalg.norm(transformed-target) > OPERATOR_GUARD
            or np.linalg.norm(rotation.conj().T@rotation-np.eye(90)) > 2e-11):
        raise ArithmeticError('the entire phase drive failed its complete42-pair factorization')
    for phase in (.2, math.pi/2, -1.1):
        if np.linalg.norm(nambu(higgs_pairing(phase))-math.cos(phase)*base-math.sin(phase)*derivative) > OPERATOR_GUARD:
            raise ArithmeticError('actual Higgs coefficients disagree with phase interpolation')
    return {'positive_basis': basis, 'full_basis': complete, 'degenerate_rotation': rotation,
            'phase_derivative': derivative, 'pairing_in_paired_chart': desired,
            'pairs': tuple(pairs), 'zero_modes': tuple(map(int, np.flatnonzero(chart['energies'] == 0)))}


def pair_hamiltonian(energy, phase):
    if (isinstance(energy, (bool, str)) or not math.isfinite(float(energy)) or energy <= 0
            or isinstance(phase, (bool, str)) or not math.isfinite(float(phase))):
        raise ValueError('positive finite energy and finite phase required')
    value = energy*math.cos(phase)*np.diag([-1., 0., 0., 1.]).astype(complex)
    value[0, 3] = value[3, 0] = energy*math.sin(phase)
    return value


def pair_word(energy, phase, duration):
    """Source word for exp(-i duration H_pair), retaining its scalar phase."""
    pair_hamiltonian(energy, phase)
    if isinstance(duration, (bool, str)) or not math.isfinite(float(duration)):
        raise ValueError('finite common duration required')
    # R=exp(i phase Y_even/2), Y_even=(X0Y1+Y0X1)/2.
    # R (-Z_even) R^dagger = -cos(phase)Z_even+sin(phase)X_even.
    rotate = (compiler.pauli_rotation('XY', phase/4)
              + compiler.pauli_rotation('YX', phase/4))
    diagonal = [compiler.pulse('Z', duration*energy/2, [0]),
                compiler.pulse('Z', duration*energy/2, [1])]
    return compiler.inverse(rotate)+diagonal+rotate


def pair_gate(energy, phase, duration):
    word = pair_word(energy, phase, duration)
    actual = compiler.execute(word, 2)
    h = pair_hamiltonian(energy, phase)
    even = np.diag([1., 0., 0., 1.])
    target = np.eye(4)+(math.cos(energy*duration)-1)*even-1j*math.sin(energy*duration)*h/energy
    if np.linalg.norm(actual-target, 2) > SOURCE_GUARD:
        raise ArithmeticError('full four-occupation source gate differs from the driven action')
    return actual, word


@lru_cache(maxsize=1)
def fermionic_swap_word():
    """Adjacent CAR permutation: ordinary SWAP and minus on double occupancy."""
    return (compiler.cnot(0, 1)+compiler.cnot(1, 0)+compiler.cnot(0, 1)
            + compiler.compile_diagonal([0., 0., 0., math.pi]))


def routed_pair_word(energy, phase, duration, left, right):
    """Embed a pair gate in the fixed global CAR order, including spectator signs."""
    if type(left) is not int or type(right) is not int or not 0 <= left < right:
        raise ValueError('two increasing CAR mode indices required')
    route = []
    for neighbor in range(right-1, left, -1):
        route += compiler.remap(fermionic_swap_word(), [neighbor, neighbor+1])
    center = compiler.remap(pair_word(energy, phase, duration), [left, left+1])
    return route+center+compiler.inverse(route)


def bare_component_chart(indices):
    """Source-executed QP-to-bare chart for one actual two-mode component.

    In addition to diagonalizing the static mass, a commuting phase fixes
    the driven term to the same +X_even convention as ``pair_hamiltonian``.
    This prevents a chart mismatch when composing with bare propagation.
    """
    if (type(indices) is not tuple or len(indices) != 2
            or not all(type(i) is int for i in indices) or not 0 <= indices[0] < indices[1] < 90):
        raise ValueError('two increasing original internal-mode indices required')
    left, right = indices
    delta = quasiparticles()['pairing']
    if set(np.flatnonzero(delta[left])) != {right} or set(np.flatnonzero(delta[right])) != {left}:
        raise ValueError('a complete isolated two-mode pairing component required')
    z = delta[left, right]
    channels = [row for row in channel_pairings().values() if abs(row['pairing'][left, right]) > 0]
    if len(channels) != 1:
        raise ArithmeticError('two-mode component belongs to multiple Higgs phase conventions')
    sign = channels[0]['creation_phase_sign']
    # D first, then the mass diagonalizer, then the phase of the bare z.
    word = [compiler.pulse('Z', -sign*math.pi/8, [i]) for i in (0, 1)]
    word += compiler.pauli_rotation('XY', math.pi/8)+compiler.pauli_rotation('YX', math.pi/8)
    word += [compiler.pulse('Z', -float(np.angle(z))/4, [i]) for i in (0, 1)]
    actual = compiler.execute(word, 2)
    energy = float(abs(z))
    for phase in (0., .37):
        coefficient = z*np.exp(1j*sign*phase)
        bare = np.zeros((4, 4), complex)
        bare[3, 0], bare[0, 3] = coefficient, coefficient.conjugate()
        if np.linalg.norm(actual@pair_hamiltonian(energy, phase)@actual.conj().T-bare) > SOURCE_GUARD:
            raise ArithmeticError('source chart fails the same bare static/driven action join')
    return {'source_unitary': actual, 'source_word': word, 'energy': energy,
            'creation_phase_sign': sign, 'original_indices': indices}


@lru_cache(maxsize=1)
def pointer_gate():
    """Complete16-column neutral occupation pointer with both bits retained."""
    word = compiler.cnot(0, 2)+compiler.cnot(1, 3)
    actual = compiler.execute(word, 4)
    target = np.zeros((16, 16), complex)
    for index in range(16):
        target[index ^ ((index & 3) << 2), index] = 1
    if np.linalg.norm(actual-target, 2) > SOURCE_GUARD:
        raise ArithmeticError('occupation pointer source truth table failed')
    return actual, word


def observed_pair(energy, phase, duration, initial=0):
    if type(initial) is not int or initial not in range(4):
        raise ValueError('one of all four initial occupation states required')
    gate, word = pair_gate(energy, phase, duration)
    initial_with_pointer = np.zeros(16, complex)
    initial_with_pointer[:4] = gate[:, initial]
    pointer, pointer_word = pointer_gate()
    output = pointer@initial_with_pointer
    amplitudes = output.reshape(4, 4)
    probabilities = np.sum(abs(amplitudes)**2, axis=1)
    expected = abs(gate[:, initial])**2
    if (abs(float(probabilities.sum())-1) > SOURCE_GUARD
            or np.max(abs(probabilities-expected)) > SOURCE_GUARD):
        raise ArithmeticError('complete driven occupation readout lost probability')
    return {'amplitudes': amplitudes, 'probabilities': probabilities,
            'word': word+pointer_word, 'source_gate': gate}


def extended_number_law(duration, phase, current=0):
    """All number outcomes for the old neutral two-particle initial state.

    The electron-like left mode is occupied at port0 and its opposite mode
    at port1. Their two local blocks have odd occupation. Every other
    driven block starts in its vacuum. No same-port pair outcome is cut.
    """
    if type(current) is not int or current not in (-1, 0, 1):
        raise ValueError('one actual signed native clock record required')
    effective = float(duration)*(1+current/4)
    pairs = paired_chart()['pairs']
    first = quasiparticles()['ramsey_indices'][0]
    blocked = next(i for i, pair in enumerate(pairs) if pair['left'] == first)
    factors = []
    for site in range(PORTS):
        for index, pair in enumerate(pairs):
            if site in (0, 1) and index == blocked:
                continue
            chance = math.sin(pair['energy']*effective)**2*math.sin(phase)**2
            factors.append(chance)
    # Decimal exponents preserve tails far below double underflow. The
    # inputs remain explicitly numerical trigonometric probabilities.
    with localcontext() as context:
        context.prec = 70
        probabilities = [Decimal(1)]
        for chance in factors:
            p = Decimal(str(chance))
            updated = [Decimal(0)]*(len(probabilities)+1)
            for index, value in enumerate(probabilities):
                updated[index] += (1-p)*value
                updated[index+1] += p*value
            probabilities = updated
        if abs(sum(probabilities)-1) > Decimal('1e-60'):
            raise ArithmeticError('full number distribution lost probability')
        mean = sum((2+2*i)*p for i, p in enumerate(probabilities))
    numbers = tuple(2+2*i for i in range(len(probabilities)))
    return {'particle_numbers': numbers, 'probabilities_decimal': tuple(map(str, probabilities)),
            'factor_probabilities': tuple(factors),
            'probability_staying_N2_decimal': str(probabilities[0]),
            'expected_particle_number': float(mean),
            'max_pair_occupation': 2, 'odd_blocks': 2, 'driven_vacuum_blocks': PORTS*42-2,
            'effective_duration': effective}


def repair_echo(energy, duration=1.):
    """Actual source drive -> full-domain repair -> inverse-drive echo.

    Two additional particles at ports outside the repaired seam are exact
    identity spectators. Thus vacuum/pair occupations below can equally
    witness total N=2 versus N=4, starting in the old hard-core N=2 sector.
    """
    from ..defect_field_bridge import full_domain

    gate, word = pair_gate(energy, math.pi/2, duration)
    initial = np.zeros(16, complex)
    initial[:4] = gate[:, 0]
    reverse = np.kron(np.eye(4), gate.conj().T)
    chance = math.sin(energy*duration)**2
    expected = 1-2*chance*(1-chance)
    results = {}
    for law in ('C', 'E'):
        isometry, _, resources = full_domain.source_isometry(2, 2, law)
        after = (isometry@initial).reshape(32, 16)
        echo = (reverse@after.T).T
        probability = float(np.sum(abs(echo[:, 0])**2))
        if abs(probability-expected) > 4*SOURCE_GUARD:
            raise ArithmeticError('complete source repair echo differs from recorded-collision law')
        results[law] = {'all_record_field_amplitudes': echo,
                        'observed_return_probability': probability,
                        'ideal_return_probability': expected,
                        'repair_resources': resources}
    return results


def full_repair_echo(duration=.002, phase=math.pi/2, current=0):
    """Exact projector formula on the full spatial state after pair driving.

    Both seam endpoints initially are empty. Two particles may occupy
    other ports. Before the repair, each seam site has even occupation;
    its only noncollision possibility is its full local vacuum. Nothing
    assumes the post-measurement field remains Gaussian.
    """
    if type(current) is not int or current not in (-1, 0, 1):
        raise ValueError('one actual signed native clock record required')
    effective = duration*(1+current/4)
    empty_probability = math.prod(1-math.sin(pair['energy']*effective)**2*math.sin(phase)**2
                                  for pair in paired_chart()['pairs'])
    both_empty = empty_probability**2
    return {'duration': duration, 'effective_duration': effective,
            'two_endpoint_vacuum_probability': both_empty,
            'collision_probability': 1-both_empty,
            'each_coin_raw0_probability': both_empty/2,
            'each_coin_raw8_probability': (1-both_empty)/2,
            'inverse_drive_return_probability': both_empty**2+(1-both_empty)**2,
            'return_without_repair': 1,
            'formula': 'For orthogonal projectors P0,P8, return=sum_r |<psi|P_r|psi>|^2=A^2+(1-A)^2.',
            'postmeasurement_Gaussian_assumption': False}
