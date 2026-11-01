"""Regional product-potential words, reusing the admitted native encoding.

No acquisition, preparation, loader, pointer or arbitrary matrix primitive is
introduced here. Source words act on the existing nine-bit address register.
The complete channel is specified by first-failure Kraus operators below.
Ideal controls preserve the code. Finite integrated-control rounding and
every resulting complement branch enter the independent complete bound.
"""
from __future__ import annotations

from fractions import Fraction as F
import math

import numpy as np

from ..local_repair_family import learned_flow, modular_clock
from ..native_program import interpreter

REGIONS = ((1, 7), (1, 2))
BASE_PORT = 7  # A single gauge for all regions and their overlap/union.
OPERATOR_CAP = F(1, 10**9)


def region(ports):
    if (type(ports) is not tuple or any(type(p) is not int for p in ports)
            or tuple(sorted(set(ports))) != ports or any(not 0 <= p < 12 for p in ports)):
        raise ValueError('sorted distinct original ports required')
    return ports


def weights(values):
    values = learned_flow._weights(values)
    return tuple(w/values[BASE_PORT] for w in values)


def diagonal(values, ports):
    """Product-potential lift; a sector-central gauge of reduced modular time."""
    w, ports = weights(values), region(ports)
    states, codes, _ = learned_flow._source_chart()
    ratios = [F(1)]*512
    for state, address in zip(states, codes, strict=True):
        ratios[address] = math.prod(w[p]**state[p] for p in ports)
    return tuple(ratios)


def word(values, ports, parameter=F(1)):
    if type(parameter) is not F or abs(parameter) > 1:
        raise ValueError('original rational parameter in [-1,1] required')
    return modular_clock.common.diagonal_word(parameter, diagonal(values, ports))


def admitted(values, ports, parameter=F(1)):
    """Literal grammar admission; this is the bound record-clock encoding."""
    return interpreter.source_word(word(values, ports, parameter), 9)


def prepare_request(values, ports, parameter=F(1)):
    """Certify before touching a source; every rejected request is identity.

    This does not condition a quantum output. Unresolved record inference is
    handled by its existing owner before these twelve relative weights exist.
    """
    from .independent import angle_certificate
    try:
        payload = admitted(values, ports, parameter)
        certificate = angle_certificate(word(values, ports, parameter), values, ports, parameter)
    except (ValueError, ArithmeticError, TypeError, OverflowError) as error:
        return {'status': 'COMPILER_STOP', 'reason': str(error), 'pulses_executed': 0,
                'source_and_untouched_reference_action': 'identity'}
    return {'status': 'ADMITTED', 'source_word': payload, 'certificate': certificate}


def executed(values, ports, parameter=F(1)):
    tape = word(values, ports, parameter)
    array = modular_clock.common.full.compiler.execute(tape, 9)
    # Binary64 numerical execution is checked separately from the analytic
    # literal-angle proof in independent.py. No polar repair is performed.
    bound = modular_clock.common.certify_diagonal(parameter, diagonal(values, ports), array)
    if bound > OPERATOR_CAP:
        raise ArithmeticError('regional source array exceeds its fixed full-address bound')
    return array


def first_failure_kraus(tape, columns=None):
    """Lazy complete work-register instrument, including every zero abort.

    At pulse t: F[t,j] = reset(row j of exp(G_t)) V[t-1] ... V[0].
    V[t] is the proper-code block, lifted on the addressed bits. The reset
    concerns only those bits; spectators and arbitrary R are identities.
    Sender reset/receiver overwrite are the existing supplied service.
    """
    compiler = modular_clock.common.full.compiler
    tape = compiler._validate(tape, 9)
    prefix = np.eye(512, dtype=complex) if columns is None else np.array(columns, dtype=complex, copy=True)
    if prefix.ndim != 2 or len(prefix) != 512 or not np.isfinite(prefix).all():
        raise ValueError('finite source-address input columns required')
    for ordinal, event in enumerate(tape):
        full = compiler.native_pulse(event['kind'], event['angle'])
        d = 4 if event['kind'] == 'J' else 2
        for j in range(d, 6):
            local = np.zeros((d, d), complex)
            local[0] = full[j, :d]
            yield ('EXPORT_ABORT', ordinal, j), compiler._apply(prefix, local, event['targets'], 9)
        prefix = compiler._apply(prefix, full[:d, :d], event['targets'], 9)
    yield ('SUCCESS',), prefix


def complete_kraus(tape):
    """Rectangular complete H66 channel; aborts keep processor work or R.

    Processor aborts have a 512-dimensional work output; address aborts have
    a one-dimensional flag output. Labels form a direct sum, not a selected
    output. Tensor any returned matrix with I_R to retain arbitrary R.
    This also permits bounded hostile-prefix tests of the complete service.
    """
    codes = learned_flow._source_chart()[1]
    injection = np.eye(512, dtype=complex)[:, codes]
    for label, branch in first_failure_kraus(tape, injection):
        if label != ('SUCCESS',):
            yield label, branch
        else:
            yield label, branch[list(codes)]
            for address in range(512):
                if address not in codes:
                    yield ('ADDRESS_ABORT', address), branch[address:address+1]


def address_export(array, density, reference_dimension=1):
    """All final address outcomes, with R retained even on address abort.

    This bounded numerical interface accepts arbitrary operators, so tests
    may compare full channels on matrix units without positivity clipping.
    Six-level first failures are supplied separately by first_failure_kraus.
    """
    if type(reference_dimension) is not int or reference_dimension < 1:
        raise ValueError('positive original reference dimension required')
    r = reference_dimension
    if (type(array) is not np.ndarray or array.shape != (512, 512)
            or not np.isfinite(array).all() or np.shape(density) != (66*r, 66*r)
            or not np.isfinite(density).all()):
        raise ValueError('full source array and matching H66/R operator required')
    codes = learned_flow._source_chart()[1]
    outside = tuple(i for i in range(512) if i not in codes)
    good, bad = array[np.ix_(codes, codes)], array[np.ix_(outside, codes)]
    lifted = np.kron(good, np.eye(r))
    abort = np.einsum('ji,irjs->rs', bad.conj().T@bad,
                      np.asarray(density).reshape(66, r, 66, r), optimize=True)
    return {('SUCCESS',): lifted@density@lifted.conj().T, ('ADDRESS_ABORT',): abort}
