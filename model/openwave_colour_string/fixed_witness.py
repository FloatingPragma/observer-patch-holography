"""Recheck a complete integer proof frame against the actual exact source.

This is the six-quark complete-spectrum residual proof with supplied proof
coordinates instead of a fresh LAPACK frame. No stored error or interval is
accepted. The original six-quark producer and comparison rule are unchanged.
"""
from fractions import Fraction as F
import numpy as np
from mpmath import iv
from model.openwave_six_quark.certificate import frobenius_ceil, scalar


def integers(values, length, bound):
    """Canonical integer strings, with finite dimensions and explicit bounds."""
    assert isinstance(values, list) and len(values) == length
    out = []
    for value in values:
        assert isinstance(value, str)
        number = int(value)
        assert str(number) == value and abs(number) <= bound
        out.append(number)
    return out


def spectral_enclosure(exact, witness):
    n = exact.rows
    assert n > 0 and exact.cols == n
    assert exact == exact.T
    D = witness['dyadic_denominator']
    assert type(D) is int and D == 2**40
    frame = witness['dyadic_eigenframe']
    assert isinstance(frame, list) and len(frame) == n
    O = np.array([integers(row, n, D) for row in frame], dtype=object)
    E = np.array(integers(witness['dyadic_eigenvalues'], n, 2**62), dtype=object)
    assert all(E[i] <= E[i+1] for i in range(n-1))
    h = np.asarray(exact.evalf(18), float)
    assert np.all(np.isfinite(h)) and np.max(abs(h))*D < 2**62
    H = np.rint(h*D).astype(np.int64)
    iv.dps = 55
    scalar.cache_clear()
    entry = F(0)
    for i in range(n):
        for j in range(n):
            if exact[i, j] == 0:
                continue
            x = scalar(exact[i, j])
            lo = F.from_float(float(np.nextafter(float(x.a), -np.inf)))
            hi = F.from_float(float(np.nextafter(float(x.b), np.inf)))
            mid = F(int(H[i, j]), D)
            entry = max(entry, abs(lo-mid), abs(hi-mid))
    delta = n*entry
    eta = F(frobenius_ceil(O.T@O-D*D*np.eye(n, dtype=object)), D*D)
    # This certifies invertibility of the entire frame, not selected columns.
    assert eta < F(1, 2)
    residual = np.asarray(H, dtype=object)@O-O*E
    rho = F(frobenius_ceil(residual), D*D)+delta*(1+eta)
    K = (F(max(sum(abs(int(x)) for x in row) for row in H), D)+delta).__ceil__()+1
    error = 2*rho+4*K*eta
    assert max(abs(int(x)) for x in E) < K*D
    return {
        'complete_reduced_dimension': n, 'dyadic_denominator': D,
        'source_operator_rounding_upper': str(delta),
        'exact_Gram_error_upper': str(eta), 'exact_residual_upper': str(rho),
        'all_eigenvalue_error_upper': str(error),
        'all_eigenvalue_intervals': [[str(F(int(x), D)-error), str(F(int(x), D)+error)] for x in E],
        'dyadic_eigenframe': [[str(int(x)) for x in row] for row in O],
        'dyadic_eigenvalues': [str(int(x)) for x in E],
    }
