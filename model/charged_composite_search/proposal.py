"""Untrusted complete eigensystem proposal, admitted only by integer checks."""
import numpy as np
from scipy import linalg

from . import source
from .certificate import VECTOR_BITS, ENERGY_BITS


def propose():
    block = source.block()
    energies, vectors = linalg.eigh(block['numerator'].astype(float)/16,
                                   block['gram'].astype(float))
    if not np.isfinite(energies).all() or not np.isfinite(vectors).all():
        raise ArithmeticError('nonfinite eigensolver proposal')
    return {'columns': block['columns'],
            'energies': np.array([round(float(v)*(1 << ENERGY_BITS)) for v in energies], dtype=np.int64),
            'vectors': np.array([[round(float(v)*(1 << VECTOR_BITS)) for v in row]
                                 for row in vectors], dtype=np.int64)}
