"""Exact complete-basis spectral certification of the compact source matrix."""
from fractions import Fraction as F

import numpy as np
from scipy import linalg as la, sparse as sp

from model.certified_winding_core import intervals as iv
from model.certified_winding_core.spectrum import certify
from .spectral import SIGMA, source
from .topology import direction


def rounded_source(radius=3, mass=64, half_separation=1, bits=26):
    if type(bits) is not int or not 8 <= bits <= 26:
        raise ValueError('dyadic precision from eight through twenty-six required')
    value = source(radius, mass, half_separation)
    n, factor = value['h'].shape[0], 1 << bits
    hr = sp.lil_matrix((n, n), dtype=np.int64)
    hi = sp.lil_matrix((n, n), dtype=np.int64)
    indices = {p: i for i, p in enumerate(value['points'])}
    error_squared = F(0)
    for p, i in indices.items():
        vector = direction(p, half_separation)
        norm2 = sum(x*x for x in vector)
        if norm2:
            denominator = iv.sqrt_interval(norm2)
            rounded = [iv.dyadic_center(iv.scale(iv.divide((F(x), F(x)), denominator), mass), bits)
                       for x in vector]
        else:
            rounded = [(0, F(0))]*3
        (x, ex), (y, ey), (z, ez) = rounded
        error_squared += 2*(ex*ex+ey*ey+ez*ez)
        hr[2*i, 2*i], hr[2*i+1, 2*i+1] = z, -z
        hr[2*i, 2*i+1] = hr[2*i+1, 2*i] = x
        hi[2*i, 2*i+1], hi[2*i+1, 2*i] = -y, y
        for axis in range(3):
            q = tuple(x+int(j == axis) for j, x in enumerate(p))
            if q not in indices:
                continue
            j = indices[q]
            for u in range(2):
                for v in range(2):
                    z = 0.5j*SIGMA[axis][u, v]*factor
                    hr[2*i+u, 2*j+v] = hr[2*j+v, 2*i+u] = int(z.real)
                    hi[2*i+u, 2*j+v], hi[2*j+v, 2*i+u] = int(z.imag), -int(z.imag)
    return hr.tocsr(), hi.tocsr(), error_squared


def complete_certificate(radius=3, mass=64, half_separation=1):
    bits = 26
    hr, hi, error = rounded_source(radius, mass, half_separation, bits)
    eigenvalues, vectors = la.eigh((hr.toarray()+1j*hi.toarray())/(1 << bits), driver='evr')
    e, qr, qi = (np.rint(x*(1 << bits)).astype(np.int64) for x in (eigenvalues, vectors.real, vectors.imag))
    result = certify(hr, hi, e, qr, qi, error, bits)
    intervals = [[F(x) for x in pair] for pair in result['ordered_energy_intervals']]
    low = [pair for pair in intervals if -2 < pair[0] <= pair[1] < 2]
    if len(low) != 4 or any(a <= 0 <= b for a, b in low):
        raise ArithmeticError('complete exact certification did not separate the four nonzero levels')
    if any(not (-2 < a <= b < 2 or b < -2 or a > 2) for a, b in intervals):
        raise ArithmeticError('complete exact certification did not resolve the core window')
    result.update(radius=radius, mass=mass, half_separation=half_separation,
                  exact_zero_mode_count=0, exact_core_band_dimension=4,
                  certified_core_abs_energy_lower=str(min(min(abs(a), abs(b)) for a, b in low)),
                  all_finite_source_eigenvalues_certified=True)
    return result
