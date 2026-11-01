"""Independent integer Horner replay of the consumed polynomial reference."""
from fractions import Fraction as F
import math

import numpy as np

from ..full_fock_field.field import FullFock
from ..plaquette_reference import preparation as p


def certify(data=None):
    data = p.design() if data is None else data
    h = FullFock(8).hamiltonian(p.GEOMETRY).toarray()
    if not np.array_equal(h, data.hamiltonian):
        raise ArithmeticError('prepared polynomial is not bound to the complete original action')
    frame = data.frame
    if frame.shape != (256, 256):
        raise ArithmeticError('complete Fock frame required')
    entries = []
    for row in frame:
        where = np.flatnonzero(row)
        if len(where) != 1 or row[where[0]] not in (1, -1, 1j, -1j):
            raise ArithmeticError('the exact onsite frame must be a signed Gaussian permutation')
        entries.append((int(where[0]), complex(row[where[0]])))
    if len({j for j, _ in entries}) != 256:
        raise ArithmeticError('onsite frame drops a full occupation input')
    rows = []
    maximum = F()
    for i, (a, phase_a) in enumerate(entries):
        row = []
        for j, (b, phase_b) in enumerate(entries):
            # Each operation is a sign/coordinate permutation of original
            # dyadic entries; no approximate matrix reduction is performed.
            z = phase_a*h[a, b]*phase_b.conjugate()-(1 if i == j else 0)
            if z.imag or 16*z.real != int(16*z.real):
                raise ArithmeticError('original action does not have the exact real denominator16')
            n = int(16*z.real)
            if n: row.append((j, n))
        maximum = max(maximum, F(sum(abs(n) for _, n in row), 16))
        rows.append(row)
    if maximum > 3:
        raise ArithmeticError('original action exceeds the exact centered norm premise')
    d = 1024
    denominator = d**8*math.factorial(8)
    identity = np.eye(256, dtype=object)
    value = identity.copy()
    # Independent Horner recurrence, unlike the producer's summed powers.
    for degree in reversed(range(8)):
        following = np.zeros((256, 256), dtype=object)
        for i, row in enumerate(rows):
            for j, n in row: following[i] += n*value[j]
        following += (-1)**degree*d**(8-degree)*(math.factorial(8)//math.factorial(degree))*identity
        value = following
    if denominator != data.denominator or not np.array_equal(value, np.array(data.numerator, dtype=object)):
        raise ArithmeticError('integer Horner polynomial differs from the prepared columns')
    norms = tuple(sum(int(x)**2 for x in value[:, j]) for j in range(256))
    if norms != data.column_norm_squared or sum(norms) != data.total_norm_squared or not all(norms):
        raise ArithmeticError('the preparation probabilities are not the complete polynomial column norms')
    return {'all_complete_Fock_entries_checked': 256**2,
            'independent_integer_Horner_steps': 8,
            'centered_source_row_norm_bound': maximum,
            'polynomial_denominator': denominator,
            'all_column_probabilities_exactly_checked': 256,
            'source_eigenbasis_or_floating_residual_used': False}
