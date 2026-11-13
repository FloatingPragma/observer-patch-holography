from fractions import Fraction as F
import numpy as np
from scipy.linalg import eigh_tridiagonal
from .spectrum import full_eigen_interval, jacobi, full_inertia_bounds, record

def test_exact_infinite_tail_encloses_independent_large_cap():
    for g in (F(1, 2), F(8)):
        d, b = jacobi(128, g)
        e = eigh_tridiagonal(np.array(d, float), np.full(128, float(b)), eigvals_only=True)
        for k in (0, 1):
            lo, hi, n = full_eigen_interval(k, g)
            assert float(lo)-1e-12 < e[k] < float(hi)+1e-12
            assert hi-lo < F(1, 2**40)
            assert full_inertia_bounds(lo, n, g) == (k, k)

def test_curvature_has_full_tail_error_and_noncollapse():
    row = record(F(1))
    assert F(row['curvature_variance_interval_exact'][0]) > 0
    assert len(row['all_65_Galerkin_centered_curvature_weights']) == 65
    assert all(r['infinite_source_absolute_error_upper'] < 1e-7 for r in row['full_source_Euclidean_correlations'])
