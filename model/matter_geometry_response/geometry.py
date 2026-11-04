"""Exact spectral derivatives and rational finite-frame enclosures."""
from fractions import Fraction as F
from math import isqrt

from model.common_source.slow_geometry import source

SIZE = 12
PRECISION = 80
ITERATIONS = 9


def matmul(a, b):
    return [[sum((a[i][k]*b[k][j] for k in range(len(b))), 0)
             for j in range(len(b[0]))] for i in range(len(a))]


def laplacian(theta=F(0)):
    if type(theta) is not F or abs(theta) > F(1, 8):
        raise ValueError('an exact parameter in [-1/8,1/8] is required')
    matrix = [[F(0) for _ in range(SIZE)] for _ in range(SIZE)]
    for a, b in source().seams():
        weight = 1+theta if (a, b) == (0, 1) else 1-theta if (a, b) == (0, 5) else F(1)
        matrix[a][a] += weight
        matrix[b][b] += weight
        matrix[a][b] -= weight
        matrix[b][a] -= weight
    return matrix


def exact_derivative():
    """Differentiate the complete isolated band over Q(sqrt(5))."""
    c = source()
    p = c.slow_band_projector_exact()
    bands = [c.band_projector_exact(i) for i in range(4)]
    eigenvalues = [c.q5(0), c.q5(5, -1), c.q5(6), c.q5(5, 1)]
    d = [[c.q5(0) for _ in range(SIZE)] for _ in range(SIZE)]
    for (a, b), sign in (((0, 1), 1), ((0, 5), -1)):
        for i, vi in ((a, 1), (b, -1)):
            for j, vj in ((a, 1), (b, -1)):
                d[i][j] = c.q5_add(d[i][j], c.q5(sign*vi*vj))
    dp = [[c.q5(0) for _ in range(SIZE)] for _ in range(SIZE)]
    for k in (0, 2, 3):
        coefficient = c.q5_inv(c.q5_sub(eigenvalues[1], eigenvalues[k]))
        term = c.q5_matmul(c.q5_matmul(bands[k], d), p)
        for i in range(SIZE):
            for j in range(SIZE):
                dp[i][j] = c.q5_add(dp[i][j], c.q5_mul(
                    coefficient, c.q5_add(term[i][j], term[j][i])))
    gram = c.intrinsic_gram_exact()
    dg = [[c.q5_sub(c.q5_mul(c.q5(4), dp[i][j]),
            c.q5_mul(c.q5(2), c.q5_mul(gram[i][j],
                c.q5_add(dp[i][i], dp[j][j]))))
           for j in range(SIZE)] for i in range(SIZE)]
    if dg[0][1] != c.q5(F(-33, 100), F(29, 100)):
        raise ArithmeticError('native normalized Gram derivative changed')
    if c.q5_sign(c.q5_sub(dg[0][1], c.q5(F(3, 10)))) <= 0:
        raise ArithmeticError('strict Gram response bound failed')
    return dp, dg


def _nearest(numerator, denominator):
    return (2*numerator+denominator)//(2*denominator)


def _sqrt_bounds(value, digits=100):
    if value < 0:
        raise ValueError('nonnegative rational square root required')
    scale = 1 << digits
    floor = isqrt(value.numerator*scale*scale//value.denominator)
    lower = F(floor, scale)
    upper = lower if lower*lower == value else F(floor+1, scale)
    return lower, upper


def enclosed_projector(theta, *, precision=PRECISION, iterations=ITERATIONS):
    """Integer-rounded sign iteration with an analytic operator error bound."""
    if type(precision) is not int or not 40 <= precision <= 120:
        raise ValueError('fixed-point precision must be an integer from 40 to 120')
    if type(iterations) is not int or not 1 <= iterations <= 10:
        raise ValueError('one to ten sign iterations are required')
    matrix = laplacian(theta)
    scale = 1 << precision
    exact = [[(2*matrix[i][j]-9*int(i == j))/9 for j in range(SIZE)]
             for i in range(SIZE)]
    x = [[_nearest(value.numerator*scale, value.denominator) for value in row]
         for row in exact]
    for _ in range(iterations):
        cube = matmul(matmul(x, x), x)
        x = [[_nearest(3*scale*scale*x[i][j]-cube[i][j], 2*scale*scale)
              for j in range(SIZE)] for i in range(SIZE)]
    rounding = F(2*(4**(iterations+1)-1), scale)
    if rounding > F(1, 4):
        raise ArithmeticError('sign iteration rounding left its Lipschitz domain')
    error = (F(15, 16)**(2**iterations)+rounding)/2
    projector = [[(F(int(i == j))-F(x[i][j], scale))/2-F(1, SIZE)
                  for j in range(SIZE)] for i in range(SIZE)]
    return projector, error


def enclosed_gram(theta):
    p, error = enclosed_projector(theta)
    result = []
    for i in range(SIZE):
        row = []
        for j in range(SIZE):
            if i == j:
                row.append((F(1), F(1)))
                continue
            a, b = p[i][i], p[j][j]
            if min(a, b)-error <= 0:
                raise ArithmeticError('positive coordinate norm was not certified')
            low = _sqrt_bounds((a-error)*(b-error))[0]
            high = _sqrt_bounds((a+error)*(b+error))[1]
            quotients = [(p[i][j]+s*error)/denominator
                         for s in (-1, 1) for denominator in (low, high)]
            row.append((min(quotients), max(quotients)))
        result.append(row)
    return result, error


def certificate():
    c = source()
    dp, dg = exact_derivative()
    endpoints = {}
    intervals = {}
    for energy, theta in ((0, F(3, 25)), (1, F(-3, 25))):
        gram, error = enclosed_gram(theta)
        intervals[energy] = gram
        endpoints[str(energy)] = {
            'theta': str(theta), 'projector_operator_error_upper': str(error),
            'normalized_Gram_intervals': [[[str(lo), str(hi)] for lo, hi in row]
                                          for row in gram],
        }
    contrast = intervals[0][0][1][0]-intervals[1][0][1][1]
    if contrast <= F(3, 40):
        raise ArithmeticError('finite native Gram contrast was not certified')
    return {
        'source_laplacian_spectrum': ['0 (1)', '5-sqrt(5) (3)', '6 (5)', '5+sqrt(5) (3)'],
        'physical_theta_interval': ['-3/25', '3/25'],
        'theorem_theta_interval': ['-1/8', '1/8'],
        'band_rank': 3, 'low_to_next_band_gap_lower': 'sqrt(5)+1/2',
        'projector_change_operator_norm_upper': '1/5',
        'projector_diagonal_lower': '1/20',
        'normalized_Gram_diagonal': 1,
        'exact_projector_derivative': [[c.q5_str(x) for x in row] for row in dp],
        'exact_normalized_Gram_derivative': [[c.q5_str(x) for x in row] for row in dg],
        'derivative_G01_at_zero': ['-33/100', '29/100'],
        'derivative_G01_strict_lower': '3/10',
        'sign_iteration': {'scale_bits': PRECISION, 'iterations': ITERATIONS,
                           'squared_sign_defect_initial_upper': '15/16'},
        'energy_endpoint_frames': endpoints,
        'G01_energy0_minus_energy1_strict_lower': '3/40',
        'finite_contrast_certified_lower': str(contrast),
        'numerical_eigenvectors_used_in_certificate': False,
    }
