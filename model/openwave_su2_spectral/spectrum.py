"""Infinite Jacobi tails are bounded before any finite eigenvalue is accepted."""
from fractions import Fraction as F
from functools import lru_cache
from math import factorial, isqrt, lcm, nextafter, inf
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh_tridiagonal
from model.openwave_pure_su2_source.source import coefficients


@lru_cache(None)
def witnesses():
    return json.loads(Path(__file__).with_name('numerical_witnesses.json').read_text())


def sqrt_upper(x, bits=70):
    x = F(x)
    if x < 0:
        raise ValueError('negative radicand')
    k = isqrt((x.numerator << (2*bits))//x.denominator)
    return F(k+1, 1 << bits)


def exp_bounds(x):
    x = F(x)
    if x < 0:
        lo, hi = exp_bounds(-x)
        return 1/hi, 1/lo
    n = max(80, int(x)+80)
    term = total = F(1)
    for j in range(1, n+1):
        term *= x/j; total += term
    nxt = term*x/(n+1)
    return total, total+nxt/(1-x/F(n+2))


def jacobi(n, g2=F(8), spacing=F(1)):
    alpha, beta = coefficients(g2, spacing)
    diagonal = [alpha*k*(k+2)+beta for k in range(n+1)]
    return diagonal, -beta/2


def sturm(diagonal, off, x, last_shift=F(0)):
    prev, current = F(0), F(1)
    signs = [1]
    for k, d in enumerate(diagonal):
        value = (d-x-(last_shift if k == len(diagonal)-1 else 0))*current-off*off*prev
        if value:
            signs.append(1 if value > 0 else -1)
        prev, current = current, value
    return sum(a != b for a, b in zip(signs[:-1], signs[1:], strict=True))


def full_inertia_bounds(x, n, g2, spacing=F(1)):
    alpha, beta = coefficients(g2, spacing)
    floor = alpha*(n+1)*(n+3)
    if x >= floor:
        raise ValueError('tail positivity is required before the Schur reduction')
    d, b = jacobi(n, g2, spacing)
    lo = sturm(d, b, x)
    hi = sturm(d, b, x, b*b/(floor-x))
    assert lo <= hi
    return lo, hi


def full_eigen_interval(k, g2, spacing=F(1), bits=40):
    d, b = jacobi(64, g2, spacing)
    proposals = eigh_tridiagonal(np.array(d, float), np.full(64, float(b)), eigvals_only=True)
    center = F(round(proposals[k]*2**24), 2**24)
    lo, hi = center-F(1, 100), center+F(1, 100)
    n = 16
    while full_inertia_bounds(lo, n, g2, spacing) != (k, k) or full_inertia_bounds(hi, n, g2, spacing) != (k+1, k+1):
        n += 8
        assert n <= 96
    while hi-lo > F(1, 2**bits):
        mid = (lo+hi)/2
        count = full_inertia_bounds(mid, n, g2, spacing)
        while count[0] != count[1]:
            n += 8
            assert n <= 128
            count = full_inertia_bounds(mid, n, g2, spacing)
        if count[0] <= k:
            lo = mid
        else:
            hi = mid
    assert full_inertia_bounds(lo, n, g2, spacing) == (k, k)
    assert full_inertia_bounds(hi, n, g2, spacing) == (k+1, k+1)
    return lo, hi, n


def rational_trial(g2, n=64):
    d, b = jacobi(n, g2)
    if n != 64:
        raise ValueError('the accepted complete witness has65 character coordinates')
    w = witnesses()[str(g2)]
    assert w['ground_denominator'] == 2**48 and w['basis_denominator'] == 2**52 and w['energy_denominator'] == 2**44
    e = np.array(w['energy_integer_coordinates'], float)/2**44
    v = np.array(w['basis_integer_coordinates'], float)/2**52
    u = [F(z, 2**48) for z in w['ground_integer_coordinates']]
    assert len(e) == len(u) == 65 and v.shape == (65, 65)
    norm = sum(z*z for z in u)
    Hu = [d[k]*u[k]+(b*u[k-1] if k else 0)+(b*u[k+1] if k < n else 0) for k in range(n+1)]
    mean = sum(x*y for x, y in zip(u, Hu, strict=True))/norm
    residual2 = (sum((h-mean*z)**2 for h, z in zip(Hu, u, strict=True))+(b*u[-1])**2)/norm
    W = [(u[k-1] if k else 0)/2+(u[k+1] if k < n else 0)/2 for k in range(n+1)]
    wm = sum(x*y for x, y in zip(u, W, strict=True))/norm
    centered = [w-wm*z for w, z in zip(W, u, strict=True)]
    variance = (sum(z*z for z in centered)+(u[-1]/2)**2)/norm
    return d, b, e, v, u, norm, mean, residual2, wm, centered, variance


def exact_basis_error(d, b, values, vectors):
    """Polar-factor enclosure from exact integer Gram and residual matrices."""
    n = len(d)
    vi = np.array([[round(float(x)*2**52) for x in row] for row in vectors], dtype=object)
    ei = np.array([round(float(x)*2**44) for x in values], dtype=object)
    gram = vi.T@vi-(2**104)*np.eye(n, dtype=object)
    eta = sqrt_upper(F(sum(int(x)**2 for x in gram.flat), 2**208))
    den = lcm(*(x.denominator for x in d), b.denominator)
    hvi = np.array([int(z*den)*vi[k] for k, z in enumerate(d)], dtype=object)
    hvi[1:] += int(b*den)*vi[:-1]
    hvi[:-1] += int(b*den)*vi[1:]
    residual = hvi*2**44-den*vi*ei[None, :]
    r = sqrt_upper(F(sum(int(x)**2 for x in residual.flat), den*den*2**192))
    assert eta < F(1, 2)
    energy = max(abs(F(int(x), 2**44)) for x in ei)
    # Q=V(V*V)^(-1/2); both exact source and Q D Q* are positive.
    generator_error = 2*r+4*energy*eta
    return eta, generator_error, vi, ei


def record(g2):
    intervals = [full_eigen_interval(k, g2) for k in range(4)]
    d, b, ev, V, u, norm, mean, residual2, wm, centered, variance = rational_trial(g2)
    separation = intervals[1][0]-mean
    assert separation > 0
    distance = sqrt_upper(residual2)/separation
    assert distance < F(1, 10**9)
    varlo, varhi = max(F(0), variance-6*distance), min(F(1), variance+6*distance)
    assert varlo > 0
    v = np.array(centered, float)/np.sqrt(float(norm))
    weights = abs(V.T@v)**2
    orthogonality, generator_error, vi, ei = exact_basis_error(d, b, ev, V)
    exact_weights = [sum(F(int(vi[k, j]), 2**52)*centered[k] for k in range(65))**2/norm
                     for j in range(65)]
    # Full finite Galerkin spectral measure, including the residual vacuum overlap.
    assert abs(weights.sum()-float(sum(z*z for z in centered)/norm)) < 1e-12
    rows = []
    alpha, beta = coefficients(g2)
    m = 16
    tail = sqrt_upper(sum(z*z for z in centered[m+1:])/norm)
    initial_outside = abs(u[-1])/2*sqrt_upper(1/norm)
    for time in (F(0), F(1, 4), F(1, 2), F(1), F(2)):
        # Contractive diagonal Dyson expansion must make n-m hops to the cut.
        reach = exp_bounds(beta*time)[1]*(beta*time)**(64-m)/factorial(64-m)
        derror = 2*abs(b)*time*(tail+2*reach)+4*initial_outside
        energy_uncertainty = max(abs(mean-intervals[0][0]), abs(mean-intervals[0][1]))
        amplification = exp_bounds(time*max(mean, intervals[0][1]))[1]
        finite_evolution_error = time*generator_error+3*orthogonality
        error = 12*distance+amplification*(derror+time*energy_uncertainty+finite_evolution_error)
        scalar_lo = scalar_hi = F(0)
        for j, weight in enumerate(exact_weights):
            exponent = time*(F(int(ei[j]), 2**44)-mean)
            bounds = (F(0), F(1, 2**92)) if exponent > 64 else exp_bounds(-exponent)
            scalar_lo += weight*bounds[0]; scalar_hi += weight*bounds[1]
        value = float((scalar_lo+scalar_hi)/2)
        assert abs(value-float(np.sum(weights*np.exp(-float(time)*(ev-float(mean)))))) < 1e-10
        full_lo, full_hi = max(F(0), scalar_lo-error), scalar_hi+error
        rows.append({'time_exact': str(time), 'finite_full_spectral_value': value,
                     'infinite_source_absolute_error_upper': float(error+(scalar_hi-scalar_lo)/2),
                     'infinite_source_interval': [max(0., nextafter(float(full_lo), -inf)), nextafter(float(full_hi), inf)],
                     'exponential_evaluation': 'exact rational Taylor bounds; for exponent>64, 0<=exp(-x)<=2^-92'})
    gaplo = intervals[1][0]-intervals[0][1]
    gaphi = intervals[1][1]-intervals[0][0]
    assert gaplo > 0
    return {'g_squared_exact': str(g2), 'spacing_exact': '1',
            'full_Haar_character_eigenvalue_intervals_exact': [[str(a), str(b)] for a, b, _ in intervals],
            'Schur_tail_cutoffs': [n for _, _, n in intervals],
            'full_gap_interval_exact': [str(gaplo), str(gaphi)],
            'full_gap_interval_decimal': [float(gaplo), float(gaphi)],
            'ground_trial_integer_coordinates': [int(z*2**48) for z in u],
            'ground_trial_denominator': 2**48,
            'full_residual_squared_exact': str(residual2),
            'finite_spectral_basis_orthogonality_error_upper_exact': str(orthogonality),
            'finite_spectral_generator_error_upper_exact': str(generator_error),
            'true_ground_projector_distance_upper_exact': str(distance),
            'curvature_variance_interval_exact': [str(varlo), str(varhi)],
            'curvature_centered_energy_upper_exact': str(alpha),
            'all_65_Galerkin_energy_approximations': ev.tolist(),
            'all_65_Galerkin_centered_curvature_weights': weights.tolist(),
            'full_source_Euclidean_correlations': rows,
            'complete_positive_transfer_spectral_measure': True,
            'Galerkin_lines_identified_with_all_infinite_lines': False}


def run():
    rows = [record(g) for g in (F(1, 2), F(1), F(4), F(8))]
    return {'schema': 'oph.untruncated_SU2_plaquette_spectral.v1',
            'verdict': 'PASS_FULL_SU2_CHARACTER_GAPS_AND_CURVATURE_CORRELATIONS',
            'source_domain': 'one four-link plaquette, exact Gauss-invariant full Haar class-function space',
            'same_coefficients_as_pure_cube': True,
            'character_basis': 'chi_n(theta)=sin((n+1)theta)/sin(theta), every n>=0',
            'character_matrix': 'H_nn=alpha*n*(n+2)+beta; H_n,n+1=-beta/2',
            'executions': rows,
            'spacing_identity': 'at fixed g squared, H(a)=H(1)/a and every nonzero gap scales as 1/a',
            'scope': {'finite_representation_tail_certified': True,
                      'selected_probe_alone_proves_dense_decay': False,
                      'infinite_volume_or_four_dimensional_continuum_gap': False,
                      'physical_GeV_glueball_mass_claimed': False,
                      'Clay_problem_solved': False}}
