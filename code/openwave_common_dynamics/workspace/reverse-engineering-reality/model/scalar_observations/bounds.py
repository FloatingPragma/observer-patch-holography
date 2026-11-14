"""Rational bounds for one fixed seven-mode interacting observation family.

The finite-code and unrestricted spatial comparisons have different domains.
Neither estimates finite occupation error.  Time is Hamiltonian model time;
these bounds concern ideal operations, not floating or physical services.
"""
from __future__ import annotations

from fractions import Fraction as F
from math import factorial, isqrt


MODES = 7
OMEGA = F(2)


def _rational(value, name, *, positive=False):
    if type(value) not in (int, F) or value < 0 or (positive and value == 0):
        raise ValueError(f"{name} must be an exact nonnegative rational, positive where required")
    return F(value)


def _dimensions(dimensions):
    if type(dimensions) is int:
        dimensions = (dimensions,) * MODES
    else:
        dimensions = tuple(dimensions)
    if len(dimensions) != MODES or any(type(d) is not int or d < 2 for d in dimensions):
        raise ValueError("seven integer occupation dimensions of at least two required")
    return dimensions


def sqrt_upper(value, digits=36):
    """Directed rational square root, checked without floating arithmetic."""
    value = _rational(value, "radicand")
    if type(digits) is not int or digits < 1:
        raise ValueError("positive integer root precision required")
    scale = 10**digits
    numerator = isqrt(value.numerator * scale * scale // value.denominator)
    result = F(numerator, scale)
    if result * result < value:
        result += F(1, scale)
    if result * result < value:
        raise ArithmeticError("directed square-root enclosure failed")
    return result


def _coordinate_power(d, power):
    """Bargmann matrix of the full Q**power before compression, Omega=2."""
    matrix = [[F(0) for _ in range(d)] for _ in range(d)]
    for column in range(d):
        terms = {column: F(1)}
        for _ in range(power):
            following = {}
            for occupation, coefficient in terms.items():
                following[occupation + 1] = following.get(occupation + 1, F(0)) + coefficient / 2
                if occupation:
                    following[occupation - 1] = following.get(occupation - 1, F(0)) + coefficient * occupation / 2
            terms = following
        for row, coefficient in terms.items():
            if row < d:
                matrix[row][column] = coefficient
    return matrix


def _row_norm_upper(matrix):
    # The input monomial |n>=(a*)**n|0> has norm sqrt(n!).  Transform to
    # normalized Fock entries before applying the row-sum bound.  All uses
    # are Hermitian or anti-Hermitian, so absolute entries are symmetric.
    d = len(matrix)
    return max(sum((abs(matrix[i][j]) * sqrt_upper(F(factorial(i), factorial(j)))
                    for j in range(d)), F(0)) for i in range(d))


def coordinate_square_norm_upper(d):
    """Bound ||Pi_d Q**2 Pi_d||; this is not (Pi_d Q Pi_d)**2."""
    if type(d) is not int or d < 2:
        raise ValueError("integer occupation dimension at least two required")
    return _row_norm_upper(_coordinate_power(d, 2))


def _centered_onsite(d, frequency_squared, quartic_coefficient):
    q2, q4 = _coordinate_power(d, 2), _coordinate_power(d, 4)
    out = [[frequency_squared * q2[i][j] / 2 + quartic_coefficient * q4[i][j]
            for j in range(d)] for i in range(d)]
    # P**2/2 at Omega=2 has diagonal n+1/2, and raising/lowering
    # monomial coefficients -1/2 and -n(n-1)/2 respectively.
    for n in range(d):
        out[n][n] += n + F(1, 2)
        if n + 2 < d:
            out[n + 2][n] -= F(1, 2)
        if n >= 2:
            out[n - 2][n] -= F(n * (n - 1), 2)
    scalar = F(1, 2) + frequency_squared / 8 + 3 * quartic_coefficient / 16
    for n in range(d):
        out[n][n] -= scalar
    return out


def onsite_norm_upper(d, mode, coupling=1, mass=1):
    """Uniform centered onsite norm for every 0<=kappa<=1, Omega=2."""
    if type(d) is not int or d < 2 or type(mode) is not int or not 0 <= mode < MODES:
        raise ValueError("valid local dimension and seven-mode index required")
    coupling = _rational(coupling, "coupling")
    mass = _rational(mass, "mass", positive=True)
    coefficient = coupling / (24 if mode == 0 else 16)
    endpoints = (mass * mass,) if mode == 0 else (mass * mass, mass * mass + 1)
    # The centered matrix is affine in kappa.  Convexity of the norm
    # bounds its whole interval by the maximum of the two endpoint bounds.
    return max(_row_norm_upper(_centered_onsite(d, w2, coefficient)) for w2 in endpoints)


def pair_coefficient(i, j, coupling=1):
    if any(type(k) is not int for k in (i, j)) or not 0 <= i < j < MODES:
        raise ValueError("two ordered distinct seven-mode indices required")
    coupling = _rational(coupling, "coupling")
    same_axis = i > 0 and (i - 1) // 2 == (j - 1) // 2
    return coupling / (8 if same_axis else 4)


def generator_sum_bound(dimensions, coupling=1, mass=1):
    """Bound ||A-centered_scalar||+||B|| uniformly in spatial grid."""
    dims = _dimensions(dimensions)
    coupling = _rational(coupling, "coupling")
    mass = _rational(mass, "mass", positive=True)
    q2 = [coordinate_square_norm_upper(d) for d in dims]
    onsites = sum((onsite_norm_upper(d, i, coupling, mass) for i, d in enumerate(dims)), F(0))
    pairs = sum((pair_coefficient(i, j, coupling) * q2[i] * q2[j]
                 for i in range(MODES) for j in range(i + 1, MODES)), F(0))
    return onsites + pairs


def _matmul(left, right):
    d = len(left)
    return [[sum((left[i][k] * right[k][j] for k in range(d)), F(0))
             for j in range(d)] for i in range(d)]


def _commutator(left, right):
    lr, rl = _matmul(left, right), _matmul(right, left)
    return [[lr[i][j] - rl[i][j] for j in range(len(left))] for i in range(len(left))]


def commutator_bounds(dimensions, coupling=1, mass=1):
    """Exact local row bounds for the two global nested commutators.

    Put X_i=Pi Q_i**2 Pi, C_i=[a_i,X_i], D_i=[a_i,C_i] and
    E_i=[X_i,C_i].  C_i and E_i are independent of kappa because the
    kappa-dependent onsite term is proportional to X_i.  D_i is affine,
    so the two endpoints certify the entire spatial family.
    """
    dims = _dimensions(dimensions)
    coupling = _rational(coupling, "coupling")
    mass = _rational(mass, "mass", positive=True)
    qnorms, cnorms, dnorms, enorms = [], [], [], []
    for mode, d in enumerate(dims):
        q2 = _coordinate_power(d, 2)
        quartic = coupling / (24 if mode == 0 else 16)
        onsite = _centered_onsite(d, mass * mass, quartic)
        c = _commutator(onsite, q2)
        endpoints = (mass * mass,) if mode == 0 else (mass * mass, mass * mass + 1)
        qnorms.append(_row_norm_upper(q2))
        cnorms.append(_row_norm_upper(c))
        dnorms.append(max(_row_norm_upper(_commutator(_centered_onsite(d, w2, quartic), c))
                          for w2 in endpoints))
        enorms.append(_row_norm_upper(_commutator(q2, c)))
    coefficients = [[F(0) if i == j else pair_coefficient(min(i, j), max(i, j), coupling)
                     for j in range(MODES)] for i in range(MODES)]
    rnorms = [sum((coefficients[i][j] * qnorms[j] for j in range(MODES)), F(0))
              for i in range(MODES)]
    aab = (sum((dnorms[i] * rnorms[i] for i in range(MODES)), F(0))
           + sum((coefficients[i][j] * cnorms[i] * cnorms[j]
                  for i in range(MODES) for j in range(MODES)), F(0)))
    bab = sum((enorms[i] * rnorms[i]**2 for i in range(MODES)), F(0))
    return {"A_A_B_norm_upper": aab, "B_A_B_norm_upper": bab}


def product_error(dimensions, time, steps, coupling=1, mass=1):
    """Operator error of one symmetric-product evolution interval."""
    time = _rational(time, "model time")
    if type(steps) is not int or steps < 1:
        raise ValueError("positive integer symmetric step count required")
    bound = generator_sum_bound(dimensions, coupling, mass)
    nested = commutator_bounds(dimensions, coupling, mass)
    commutator_error = time**3 * (nested["A_A_B_norm_upper"] + 2 * nested["B_A_B_norm_upper"]) / (24 * steps**2)
    return min(F(2), bound**3 * time**3 / (3 * steps**2), commutator_error)


def splitting_history_bound(dimensions, times, steps, coupling=1, mass=1):
    """Complete recorded-output trace/TV error; every branch is retained."""
    times = tuple(times)
    if not times:
        raise ValueError("at least one evolution interval required")
    return min(F(1), sum((product_error(dimensions, time, steps, coupling, mass)
                          for time in times), F(0)))


def kappa_error_upper(grid):
    """1-kappa_n <= pi**2/(3n**2) < 484/(147n**2); no float input."""
    if grid is None:
        return F(0)
    if type(grid) is not int or grid < 5:
        raise ValueError("integer spatial grid at least five, or None for the fixed-band target")
    return F(484, 147 * grid * grid)


def spatial_history_bound(dimensions, grid, total_time):
    """Finite-occupation complete-history comparison with kappa=1.

    All preparation, instruments, feedback and times must agree.  The
    comparison target retains the same finite occupation code.
    """
    dims = _dimensions(dimensions)
    total_time = _rational(total_time, "total model time")
    generator_error = kappa_error_upper(grid) * sum(
        (coordinate_square_norm_upper(d) for d in dims[1:]), F(0)) / 2
    return min(F(1), total_time * generator_error)


def blank_energy(mass=1, coupling=1, initial_kick=0):
    """Exact kappa=1 energy of the Omega=2 product blank and one mode kick."""
    mass = _rational(mass, "mass", positive=True)
    coupling = _rational(coupling, "coupling")
    initial_kick = _rational(initial_kick, "absolute initial kick")
    return F(7, 2) + (7 * mass * mass + 6) / 8 + 49 * coupling / 128 + initial_kick**2 / 2


def unrestricted_spatial_history(grid, times, eta=2, feedback_kick=1,
                                 initial_kick=1, coupling=1, mass=1):
    """Finite-energy interacting spatial bound for two unrestricted reads.

    This compares the same seven unrestricted oscillators at kappa_n and
    kappa=1.  It does not compare them to any finite occupation execution.
    A signed initial momentum kick acts on one normalized mode; the first
    complete Ramsey read has strength eta on one mode; its committed bit
    may trigger one momentum kick before the second evolution/read.
    """
    times = tuple(_rational(t, "model time") for t in times)
    if len(times) != 2:
        raise ValueError("exactly two evolution intervals required")
    coupling = _rational(coupling, "coupling", positive=True)
    eta = _rational(eta, "detector strength")
    feedback_kick = _rational(feedback_kick, "absolute feedback kick")
    initial = blank_energy(mass, coupling, initial_kick)
    measured = initial + eta**2 / 8
    continued = (sqrt_upper(measured) + feedback_kick * sqrt_upper(F(1, 2)))**2
    delta = kappa_error_upper(grid)
    errors = [delta * times[0] * sqrt_upper(24 * initial / coupling) / 2,
              delta * times[1] * sqrt_upper(24 * continued / coupling) / 2]
    return {
        "comparison": "same seven unrestricted interacting modes, spatial gradient kappa_n versus one",
        "grid": grid,
        "model_times": list(times),
        "kappa_error_upper": delta,
        "initial_energy_exact": initial,
        "after_first_read_energy_exact": measured,
        "after_conditional_feedback_energy_upper": continued,
        "segment_vector_errors_upper": errors,
        "complete_history_trace_distance_upper": min(F(1), sum(errors, F(0))),
        "complete_history_tv_upper": min(F(1), sum(errors, F(0))),
        "positive_interaction_required": True,
        "every_record_and_disturbed_state_retained": True,
        "finite_occupation_error_included": False,
        "splitting_error_included": False,
        "mode_band_removal_included": False,
        "floating_or_physical_service_errors_included": False,
    }


def serializable(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {key: serializable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [serializable(item) for item in value]
    return value
