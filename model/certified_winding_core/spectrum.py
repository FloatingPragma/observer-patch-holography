"""Finite spectral enclosures with exact integer residuals and source intervals."""
from fractions import Fraction as F
import hashlib

import numpy as np
from scipy import linalg, sparse

from ..bounded_defect_index.finite import sites
from ..native_mass_channels.algebra import SIGMA
from ..winding_sector_search.algebra import winding
from . import intervals as iv


def source(radius, degree, mass, bits=26):
    degree = winding(degree)
    if (type(radius) is not int or not 1 <= radius <= 5 or type(mass) is not int
            or not 1 <= mass <= 4 or type(bits) is not int or not 8 <= bits <= 26):
        raise ValueError('radius1..5, integer mass1..4 and precision8..26 required')
    points = sites(radius)
    addresses = {p: i for i, p in enumerate(points)}
    factor = 1 << bits
    real, imag = {}, {}
    error_squared = F(0)

    def put(i, j, re=0, im=0):
        if re: real[i, j] = int(re)
        if im: imag[i, j] = int(im)

    for point, i in addresses.items():
        texture = iv.direction(point, degree)
        rounded = [iv.dyadic_center(iv.scale(v, mass*(-1)**(point[j] % 2)), bits)
                   for j, v in enumerate(texture)]
        (x, ex), (y, ey), (z, ez) = rounded
        error_squared += 2*(ex*ex+ey*ey+ez*ez)
        put(2*i, 2*i, z)
        put(2*i+1, 2*i+1, -z)
        put(2*i, 2*i+1, x, -y)
        put(2*i+1, 2*i, x, y)
        for axis in range(3):
            target = tuple(v+int(j == axis) for j, v in enumerate(point))
            if target in addresses:
                j = addresses[target]
                for a in range(2):
                    for b in range(2):
                        entry = complex(SIGMA[axis][a, b]*1j)
                        re, im = int(entry.real)*factor//2, int(entry.imag)*factor//2
                        put(2*i+a, 2*j+b, re, im)
                        put(2*j+b, 2*i+a, re, -im)

    def matrix(data):
        keys = list(data)
        return sparse.coo_matrix(([data[p] for p in keys],
                                 ([p[0] for p in keys], [p[1] for p in keys])),
                                shape=(2*len(points),)*2, dtype=np.int64).tocsr()
    hr, hi = matrix(real), matrix(imag)
    if (hr-hr.T).nnz or (hi+hi.T).nnz:
        raise ArithmeticError('the rational source must be exactly Hermitian')
    return hr, hi, error_squared


def integer_sum_squares(*arrays):
    return sum(int(x)*int(x) for array in arrays for x in array.flat)


def exact_dot(a, b):
    """Integer dot via 13-bit limbs with exactly representable binary64 sums.

    Every partial product and sum is an integer of magnitude below2^53.
    Exact int64 recombination is bounded separately.
    """
    if a.dtype != np.int64 or b.dtype != np.int64 or a.shape[1] != b.shape[0]:
        raise ValueError('compatible integer matrices required')
    base = 1 << 13
    ah, al = a//base, a % base
    bh, bl = b//base, b % base
    maxabs = lambda x: max(abs(int(x.min())), abs(int(x.max())))
    inner = a.shape[1]
    limbs = ((ah, bh), (ah, bl), (al, bh), (al, bl))
    limits = [inner*maxabs(x)*maxabs(y) for x, y in limbs]
    if max(limits) >= 2**53 or (base*base*limits[0]+base*(limits[1]+limits[2])+limits[3]) >= 2**63:
        raise ValueError('exact limb accumulation domain exceeded')
    if np.finfo(np.float64).nmant != 52:
        raise ArithmeticError('IEEE binary64 dot arithmetic required')
    products = []
    for x, y in limbs:
        value = x.astype(np.float64)@y.astype(np.float64)
        if not np.isfinite(value).all() or not np.array_equal(value, np.rint(value)):
            raise ArithmeticError('integer-limb dot returned noninteger data')
        products.append(value.astype(np.int64))
    return base*base*products[0]+base*(products[1]+products[2])+products[3]


def certify(hr, hi, energy, qr, qi, source_error_squared, bits):
    """Accept arbitrary proposals only through exact residual/Gram inequalities."""
    n = hr.shape[0]
    if (hr.shape != (n, n) or hi.shape != (n, n) or qr.shape != (n, n)
            or qi.shape != (n, n) or energy.shape != (n,)
            or any(x.dtype != np.int64 for x in (hr, hi, qr, qi, energy))):
        raise ValueError('complete int64 matrices and ordered diagonal required')
    if not (sparse.isspmatrix_csr(hr) and sparse.isspmatrix_csr(hi)
            and hr.has_canonical_format and hi.has_canonical_format):
        raise ValueError('canonical sparse source matrices required')
    if type(bits) is not int or not 1 <= bits <= 26 or F(source_error_squared) < 0:
        raise ValueError('valid precision and nonnegative source error required')
    factor = 1 << bits
    maxabs = lambda x: max(abs(int(x.min())), abs(int(x.max()))) if x.size else 0
    maximum = max(maxabs(qr), maxabs(qi))
    diagonal_max = max(map(abs, map(int, energy)))
    source_row = (int(np.max(np.diff(hr.indptr)))*maxabs(hr.data)
                  +int(np.max(np.diff(hi.indptr)))*maxabs(hi.data))
    bound = max(2*n*maximum*maximum+factor*factor,
                (source_row+diagonal_max)*maximum, 2*maxabs(hi.data), 2*maxabs(hr.data))
    if bound >= 2**63-1:
        raise ValueError('integer intermediate overflow is not excluded')
    if (hr-hr.T).nnz or (hi+hi.T).nnz or np.any(energy[1:] < energy[:-1]):
        raise ValueError('Hermitian source and sorted energies required')
    gram_re = exact_dot(qr.T, qr)+exact_dot(qi.T, qi)-factor*factor*np.eye(n, dtype=np.int64)
    cross = exact_dot(qr.T, qi)
    gram_im = cross-cross.T
    residual_re = hr@qr-hi@qi-qr*energy
    residual_im = hr@qi+hi@qr-qi*energy
    gram_squared = F(integer_sum_squares(gram_re, gram_im), factor**4)
    residual_squared = F(integer_sum_squares(residual_re, residual_im), factor**4)
    g = iv.sqrt_interval(gram_squared, 40)[1]
    r = iv.sqrt_interval(residual_squared, 40)[1]
    source_error = iv.sqrt_interval(source_error_squared, 40)[1]
    if g >= 1:
        raise ArithmeticError('the proposed full basis is not certified invertible')
    largest = F(diagonal_max, factor)
    error = source_error+(r+2*largest*(1+g)*g)/(1-g)
    levels = [F(int(x), factor) for x in energy]
    digest = hashlib.sha256()
    for array in (energy, qr, qi):
        digest.update(np.asarray(array, dtype='<i8').tobytes())
    return {'dimension': n, 'dyadic_bits': bits,
            'source_Frobenius_error_squared_upper': str(source_error_squared),
            'integer_Gram_defect_squared': str(gram_squared),
            'integer_intertwining_residual_squared': str(residual_squared),
            'integer_intermediate_magnitude_upper': bound,
            'Gram_defect_upper': str(g), 'residual_upper': str(r),
            'ordered_eigenvalue_error_upper': str(error),
            'ordered_energy_intervals': [[str(e-error), str(e+error)] for e in levels],
            'all_source_eigenvalues_counted': True,
            'floating_residual_used_for_certification': False,
            'proposal_sha256': digest.hexdigest()}


def calculate(radius, degree, mass, bits=26):
    hr, hi, source_error = source(radius, degree, mass, bits)
    factor = 1 << bits
    # Any numerical solver may propose a basis. Its accuracy is not a premise.
    energies, columns = linalg.eigh((hr.toarray()+1j*hi.toarray())/factor, driver='evr')
    if not np.isfinite(energies).all() or not np.isfinite(columns).all():
        raise ArithmeticError('the proposal must be finite')
    if max(np.max(np.abs(energies)), np.max(np.abs(columns))) * factor >= 2**62:
        raise ArithmeticError('proposal cannot be represented as int64')
    e, qr, qi = (np.rint(x*factor).astype(np.int64)
                 for x in (energies, columns.real, columns.imag))
    result = certify(hr, hi, e, qr, qi, source_error, bits)
    result.update(radius=radius, degree=degree, mass=mass, spacing='1', velocity='1')
    bands = [[F(x) for x in pair] for pair in result['ordered_energy_intervals']]
    inside = lambda a, b, limit: -limit < a <= b < limit
    threshold = F(mass, 4)
    count = sum(inside(a, b, threshold) for a, b in bands)
    if any(not inside(a, b, threshold) and not (b < -threshold or a > threshold) for a, b in bands):
        raise ArithmeticError('the fixed search window is not resolved')
    ambiguous = sum(a <= 0 <= b for a, b in bands)
    exact_two = degree % 2 != 0 and ambiguous == 2
    gap = min((min(abs(a), abs(b)) for a, b in bands if not a <= 0 <= b), default=F(0))
    result.update(window_absolute_energy=str(threshold), window_rank_exact=count,
                  intervals_containing_zero=ambiguous,
                  exact_kernel_dimension=2 if exact_two else (0 if ambiguous == 0 else None),
                  inversion_index=2 if degree % 2 else None,
                  nonzero_gap_lower=str(gap) if ambiguous == 0 or exact_two else None,
                  infinite_volume_kernel_or_gap_certified=False)
    return result
