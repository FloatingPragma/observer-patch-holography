"""Source coin instruments and a bounded, state-independent seam selector.

The auxiliary's current/contrast and final quantum state are discarded. The
coin, selected address, rejection count and timeout remain visible. Exact
integer Kraus supports prove the instrument identities; array executions are
separate numerical checks, not physical service-accuracy certificates.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import product

import numpy as np

from model.observer import native, quantum_mixing


def _law(law):
    if type(law) is not str or law not in ('C', 'E'):
        raise ValueError('native C or E instrument required')
    return law


def _size_rounds(size, rounds):
    if (type(size) is not int or not 2 <= size <= 256
            or type(rounds) is not int or not 1 <= rounds <= 32):
        raise ValueError('2..256 addresses and 1..32 rejection rounds required')
    return size, rounds


def coin_maps(law, seam, binary=True):
    """Integer supports of K=(1/sqrt(2)) M, independently source-enumerated."""
    _law(law)
    if (type(binary) is not bool or type(seam) is not tuple or len(seam) != 2
            or any(type(v) is not int for v in seam) or seam not in quantum_mixing.seams()):
        raise ValueError('a source seam and Boolean binary-sector switch required')
    return _coin_maps(law, seam, binary)


@lru_cache(maxsize=256)
def _coin_maps(law, seam, binary):
    states = native.states(binary=binary)
    index = {state: i for i, state in enumerate(states)}
    groups = defaultdict(list)
    for coin, (i, state) in product((0, 1), enumerate(states)):
        out, (record_coin, current) = native.native(state, seam, coin)
        if record_coin != coin:
            raise ArithmeticError('source coin record changed')
        label = current if law == 'C' else state[seam[0]]-state[seam[1]]
        groups[coin, label].append((i, index[out]))
    result = tuple((coin, label, tuple(mapping)) for (coin, label), mapping in sorted(groups.items()))
    for coin in (0, 1):
        seen = [0]*len(states)
        for c, _, mapping in result:
            if c != coin:
                continue
            sources, targets = zip(*mapping, strict=True)
            if len(set(sources)) != len(sources) or len(set(targets)) != len(targets):
                raise ArithmeticError('native recorded branch is not a partial injection')
            for source in sources:
                seen[source] += 1
        # No shared output within a branch gives zero off-diagonal effects;
        # one input occurrence gives every diagonal entry exactly 1/2.
        if seen != [1]*len(states):
            raise ArithmeticError('coarse coin effect is not exactly I/2')
    return result


def effect_certificate():
    addresses = branches = 0
    for binary, law, seam in product((False, True), ('C', 'E'), quantum_mixing.seams()):
        maps = coin_maps(law, seam, binary)
        branches += len(maps)
        addresses += sum(len(mapping) for _, _, mapping in maps)
    return {'source_dimensions': [78, 66], 'laws': ['C', 'E'], 'seams': 30,
            'partial_injection_branches_checked': branches,
            'input_coin_addresses_checked': addresses, 'coarse_coin_effect': 'I/2',
            'independent_auxiliary_initial_state_required': False,
            'auxiliary_reset_required': False,
            'auxiliary_current_and_quantum_output_in_comparison': False}


def selector_certificate(size=30, rounds=7):
    """Exact reduced CP weights, retaining the finite rejection timeout."""
    size, rounds = _size_rounds(size, rounds)
    bits = (size-1).bit_length()
    denominator = 2**bits
    reject = F(denominator-size, denominator)
    rows = [{'round': r+1, 'probability_per_address': str(reject**r/denominator)}
            for r in range(rounds)]
    timeout = reject**rounds
    accepted = sum((F(row['probability_per_address']) for row in rows), F(0))
    if size*accepted+timeout != 1:
        raise ArithmeticError('finite selector lost timeout mass')
    return {'addresses': size, 'bits_per_round': bits, 'round_cap': rounds,
            'maximum_auxiliary_native_attempts': bits*rounds,
            'expected_auxiliary_native_attempts': str(bits*sum((reject**r for r in range(rounds)), F(0))),
            'round_records': rows, 'timeout_probability': str(timeout),
            'accepted_probability_per_address': str(accepted),
            'conditional_address_probability': str(F(1, size))}


def decode_coin_word(bits, size=30):
    """A fixed controller reads only native coin records; None means reject."""
    _size_rounds(size, 1)
    if (type(bits) is not tuple or len(bits) != (size-1).bit_length()
            or any(type(bit) is not int or bit not in (0, 1) for bit in bits)):
        raise ValueError('one complete binary controller word required')
    value = 0
    for bit in bits:
        value = 2*value+bit
    return value if value < size else None


def apply_coin(density, law, coin, *, seam=None, binary=True):
    """Apply the actual coarse native coin CP map, leaving minor factors alone.

    Accepts unnormalized complex operators for linear channel checks. Density
    positivity is established by the caller's preparation and the Kraus sum.
    """
    seam = quantum_mixing.seams()[0] if seam is None else seam
    if type(coin) is not int or coin not in (0, 1):
        raise ValueError('native coin zero or one required')
    maps = coin_maps(law, seam, binary)
    d = 66 if binary else 78
    if (type(density) is not np.ndarray or density.dtype.kind not in 'fc'
            or density.ndim != 2 or density.shape[0] != density.shape[1]
            or density.shape[0] % d or not d <= density.shape[0] <= 3*d
            or not np.isfinite(density).all()):
        raise ValueError('finite plain auxiliary matrix with up to three reference levels required')
    reference = density.shape[0]//d
    source = density.reshape(d, reference, d, reference)
    result = np.zeros_like(source, dtype=complex)
    with np.errstate(over='ignore', invalid='ignore'):
        for c, _, mapping in maps:
            if c == coin:
                inputs, outputs = np.array(mapping, dtype=int).T
                result[outputs[:, None], :, outputs[None, :], :] += source[inputs[:, None], :, inputs[None, :], :]/2
    if not np.isfinite(result).all():
        raise ArithmeticError('native coin output is nonfinite')
    return result.reshape(density.shape)


def auxiliary_execution(law, rounds=7):
    """All 32 words per round on the same entangled auxiliary, without reset."""
    _law(law)
    cert = selector_certificate(rounds=rounds)
    d, ref = 66, 2
    # Select opposite orientations so the discarded current genuinely steers
    # the reference if illegally retained. The coarse coin cannot do so.
    seam = quantum_mixing.seams()[0]
    states = native.states(binary=True)
    pair = next((i, j) for i, x in enumerate(states) for j, y in enumerate(states)
                if x[seam[0]] == y[seam[1]] == 1 and x[seam[1]] == y[seam[0]] == 0
                and all(x[k] == y[k] for k in range(12) if k not in seam))
    vector = np.zeros(d*ref, complex)
    vector[pair[0]*ref], vector[pair[1]*ref+1] = .6, .8j
    initial = np.outer(vector, vector.conj())
    marginal = np.trace(initial.reshape(d, ref, d, ref), axis1=0, axis2=2)
    unresolved = initial
    records, maximum = [], 0.
    for r in range(rounds):
        branches = {(): unresolved}
        for _ in range(5):
            branches = {(*word, coin): apply_coin(rho, law, coin)
                        for word, rho in branches.items() for coin in (0, 1)}
        next_unresolved = np.zeros_like(initial)
        for word, rho in branches.items():
            address = decode_coin_word(word)
            if address is None:
                next_unresolved += rho
            else:
                expected = F(cert['round_records'][r]['probability_per_address'])
                reduced = np.trace(rho.reshape(d, ref, d, ref), axis1=0, axis2=2)
                error = float(np.linalg.norm(reduced-float(expected)*marginal, 'nuc')/2)
                maximum = max(maximum, error)
                records.append({'round': r+1, 'address': address,
                                'probability': float(np.trace(rho).real)})
        unresolved = next_unresolved
    reduced = np.trace(unresolved.reshape(d, ref, d, ref), axis1=0, axis2=2)
    maximum = max(maximum, float(np.linalg.norm(reduced-float(F(cert['timeout_probability']))*marginal, 'nuc')/2))
    total = sum(row['probability'] for row in records)+float(np.trace(unresolved).real)
    if maximum > 2e-12 or abs(total-1) > 2e-12:
        raise ArithmeticError('native auxiliary/reference selector execution failed')
    return {'law': law, 'initial_auxiliary_reference_entangled': True,
            'all_accepted_round_address_records': records,
            'timeout_probability_diagnostic': float(np.trace(unresolved).real),
            'maximum_subnormalized_reference_distance_diagnostic': maximum,
            'total_probability_residual_diagnostic': abs(total-1),
            'auxiliary_replaced_between_attempts': False}
