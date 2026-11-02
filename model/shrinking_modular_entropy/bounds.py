"""Exact constants and symbolic (never underflowed) refinement moduli."""
from fractions import Fraction as F

GAMMA_PI3 = F(4249388, 150088295)
BALL_ERROR = F(3, 14) * 195689447424 * 64 / GAMMA_PI3
ENTROPY_ERROR = F(165, 14) / GAMMA_PI3
INVERSE_NORM = F(68, 5)
RADIUS_RATIO = F(1, 2**32)


def rational(value, name):
    if type(value) not in (int, F):
        raise ValueError(name + ' must be an exact rational, not bool/float')
    return F(value)


def integer(value, name, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError(name + ' must be an integer >= ' + str(minimum))
    return value


def regularity(depth):
    """Represent K_depth and 1/(6*12**(depth+1)) without large powers."""
    j = integer(depth, 'depth')
    return {
        'K_prefactor': 2048,
        'K_8192_exponent': j,
        'K_12_exponent': j * (j + 1) // 2,
        'holder_denominator_factor': 6,
        'holder_denominator_base': 12,
        'holder_denominator_exponent': j + 1,
    }


def finite_K(depth):
    """Small-depth exact control only; production moduli stay symbolic."""
    j = integer(depth, 'depth')
    if j > 32:
        raise ValueError('materializing large regularity constants is forbidden')
    return 2048 * 8192**j * 12**(j * (j + 1) // 2)


def log_tail_bound(delta, terms):
    d = rational(delta, 'delta')
    n = integer(terms, 'terms', 1)
    if not 0 < d <= F(1, 4):
        raise ValueError('0 < delta <= 1/4 required')
    # The coarser tail avoids materializing the exponential power.
    return 1 / (n * d)


def modulus(epsilon, sobolev_bound=6):
    """Certify equation (12), preserving the positive mesh as a formula.

    The mesh upper bound is min(1/(6B), (base/K_depth)**power).
    K_depth and power are both represented by exact integer exponents.
    No floating arithmetic is used and no zero mesh is substituted.
    """
    e = rational(epsilon, 'epsilon')
    b = rational(sobolev_bound, 'Sobolev bound')
    if not 0 < e <= 1 or b < 1:
        raise ValueError('0 < epsilon <= 1 and B >= 1 required')
    d = (e / 32)**2
    required = 16 / (e * d)
    n = (required.numerator + required.denominator - 1) // required.denominator
    depth = 2*n
    endpoint = e/4
    tail = 4*log_tail_bound(d, n)
    spatial = e/4
    if endpoint + tail + spatial >= e:
        raise ArithmeticError('strict entropy-error budget failed')
    return {
        'epsilon': str(e), 'B': str(b), 'delta': str(d), 'series_terms': n,
        'regularity': regularity(depth),
        'mesh_normalization_bound': str(1/(6*b)),
        'mesh_power_base_numerator_over_K': str(e*d/(1024*b*b)),
        'mesh_power_factor': 6, 'mesh_power_base': 12,
        'mesh_power_exponent': depth+1,
        'mesh_representation': 'min(normalization_bound,(base_over_K/K_depth)^(6*12^(depth+1)))',
        'mesh_is_strictly_positive': True,
        'mesh_float_value_supplied': False,
        'endpoint_error_upper': str(endpoint),
        'series_error_upper': str(tail),
        'spatial_error_upper': str(spatial),
        'total_error_upper': str(endpoint + tail + spatial),
    }


def normalized_error(radius_ratio):
    rho = rational(radius_ratio, 'radius ratio')
    if not 0 < rho <= F(1, 8):
        raise ValueError('0 < R/L <= 1/8 required')
    return ENTROPY_ERROR*rho + BALL_ERROR*rho*rho


def neumann_budget(radius_ratio):
    return 9*INVERSE_NORM*normalized_error(radius_ratio)


def certificate():
    if not neumann_budget(RADIUS_RATIO) < F(1, 4):
        raise ArithmeticError('finite entropy matrix is not certified invertible')
    # Explicit scalar checks used by the fractional boundary estimate.
    if not 960*F(22, 7) < 3072:
        raise ArithmeticError('fractional integral bound')
    # Multiply 1+3145728/s <= 2048^2/s^2 by s^2; s<=1.
    if not 1+3145728 <= 2048**2:
        raise ArithmeticError('fractional Sobolev bound')
    return {
        'projector_Hs_constant': 3,
        'sharp_ball_Hs_to_Hs_over_12_constant': 2048,
        'regional_defect_constant': 12,
        'natural_shell_defect_constant': 8,
        'natural_shell_log_tail_constant': 2,
        'entropy_scalar_error': '8sqrt(delta)+4tau_N+256B^2*K_(2N)*eta^(1/(6*12^(2N+1)))/delta',
        'modulus_example': modulus(F(1, 2), 6),
        'modulus_formula_valid_for': '0<epsilon<=1, B>=1',
        'chosen_epsilon_and_theta': '(R/L)^5',
        'normalized_entropy_error_coefficient': str(ENTROPY_ERROR),
        'normalized_ball_error_coefficient': str(BALL_ERROR),
        'radius_ratio': str(RADIUS_RATIO),
        'matrix_Neumann_budget': str(neumann_budget(RADIUS_RATIO)),
        'matrix_inverse_infinity_norm_upper': str(F(272, 15)),
        'positive_mesh_materialized_as_float': False,
    }
