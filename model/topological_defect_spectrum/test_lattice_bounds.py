from copy import deepcopy
from fractions import Fraction as F

import numpy as np
import pytest
from scipy.linalg import eigh_tridiagonal
import sympy as s

from tools import strict_json
from . import lattice, lattice_bounds as bounds


def test_independent_source_matrix_and_exact_sturm_counts():
    for radius in (2, 3, 8):
        diagonal, off = bounds.tridiagonal(radius)
        original = lattice.jacobi(radius)
        assert diagonal == [F(original[j, j]) for j in range(original.rows)]
        assert off == [F(original[j, j+1]) for j in range(original.rows-1)]
        # Sturm's determinant recurrence is independent of division-based LDL.
        for point in (F(-7, 3), F(-1, 3), F(2, 3), F(11, 3)):
            determinants = [F(1), diagonal[0]-point]
            for j in range(1, len(diagonal)):
                determinants.append((diagonal[j]-point)*determinants[-1]-off[j-1]**2*determinants[-2])
            assert all(determinants)
            changes = sum(x*y < 0 for x, y in zip(determinants, determinants[1:], strict=False))
            assert bounds.sturm_count(diagonal, off, point) == changes
    with pytest.raises(ValueError, match='zero LDL pivot'):
        bounds.sturm_count([F(0), F(1)], [F(1)], F(0))


def test_infinite_residual_retains_the_omitted_boundary_rows():
    radius = 2
    matrix = lattice.jacobi(radius)
    kernel = matrix.nullspace()[0]
    vector = [F(x) for x in kernel]
    diagonal, off = bounds.tridiagonal(radius)
    # This exact finite-box eigenvector has zero internal residual, but a
    # strictly positive infinite-source residual at both omitted sites.
    assert matrix*kernel == s.zeros(2*radius+1, 1)
    squared = (vector[0]**2+vector[-1]**2)/(4*sum(x*x for x in vector))
    error = bounds.residual_bound(diagonal, off, F(0), vector)
    assert squared > 0
    assert error*error >= squared
    assert error == bounds.sqrt_upper(squared)
    with pytest.raises(ValueError, match='complete finite witness'):
        bounds.residual_bound(diagonal, off, F(0), vector[:-1])
    with pytest.raises(ValueError, match='nonzero witness'):
        bounds.residual_bound(diagonal, off, F(0), [F(0)]*len(vector))


def test_all_six_infinite_levels_against_independent_larger_cutoff_diagnostic():
    result = bounds.certificate()
    # This numerical check is independent evidence only; acceptance above
    # uses exact residuals and tails, and does not rely on this eigensolver.
    positions = np.arange(-64, 65)
    parity = np.where(positions % 2, -1., 1.)
    numerical, vectors = eigh_tridiagonal(positions*parity, parity[:-1]/2)
    positive = np.flatnonzero(numerical > .5)[:6]
    diagonal, off = bounds.tridiagonal(64)
    previous = F(0)
    for row, index in zip(result['positive_levels'], positive, strict=True):
        lo, hi = map(F, row['energy_interval'])
        center = F(float(numerical[index]))
        vector = [F(float(x)) for x in vectors[:, index]]
        # A floating eigenvalue can lie just outside the exact enclosure.
        # Bound its error from its actual vector with exact arithmetic,
        # including the two omitted infinite-source boundary rows.
        residual = bounds.residual_bound(diagonal, off, center, vector)
        assert residual < F(1, 10**11)
        assert previous < lo < hi
        assert center-residual <= hi and lo <= center+residual
        assert hi-lo <= F(2, 10**12)
        assert F(row['infinite_source_residual_upper']) < F(1, 10**11)
        previous = hi
    assert result['infinite_tail_retained_by_resolvent'] is True
    assert result['three_dimensional_zero_multiplicity'] == 2


def test_corrupted_witnesses_and_unprotected_endpoints_are_rejected():
    witness = strict_json.load(bounds.WITNESS)
    mutations = []
    for key, value in (('schema', 'wrong'), ('radius', 31), ('velocity', '2')):
        candidate = deepcopy(witness)
        candidate[key] = value
        mutations.append(candidate)
    candidate = deepcopy(witness)
    candidate['candidates'].pop()
    mutations.append(candidate)
    candidate = deepcopy(witness)
    candidate['candidates'][0]['positive_level'] = 2
    mutations.append(candidate)
    candidate = deepcopy(witness)
    candidate['candidates'][0]['center'] = candidate['candidates'][1]['center']
    mutations.append(candidate)
    candidate = deepcopy(witness)
    candidate['candidates'][0]['vector'][32] = '0'
    mutations.append(candidate)
    candidate = deepcopy(witness)
    candidate['candidates'][0]['vector'][0] = '1'
    mutations.append(candidate)
    for candidate in mutations:
        with pytest.raises(ValueError):
            bounds.certify(candidate)
    diagonal, off = bounds.tridiagonal(32)
    with pytest.raises(ValueError, match='tail-resolvent gap'):
        bounds.protected_endpoint(diagonal, off, 32, F(32))
    with pytest.raises(ValueError, match='spectral separation'):
        bounds.protected_endpoint(diagonal, off, 32, F(0))
    with pytest.raises(ValueError):
        bounds.sqrt_upper(F(-1))


def test_outward_rounding_encloses_both_positive_and_negative_rationals():
    for value in (F(1, 3), F(-1, 3), F(2), F(-2), F(0)):
        lower = bounds.outward_decimal(value, 12)
        upper = bounds.outward_decimal(value, 12, upper=True)
        assert lower <= value <= upper
        assert upper-lower <= F(1, 10**12)
