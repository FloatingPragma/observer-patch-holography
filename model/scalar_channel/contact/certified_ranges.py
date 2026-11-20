"""Exact enclosures of finite conditional radar means on the shared source.

The integer chart is A = 2 S^-1 H S, with S = diag(I, sqrt(2) I).
Its station basis vectors are unchanged by S. Integer moments therefore
give the station matrix elements of exp(-i n A/4). The physical H is
Hermitian with norm at most eight. The integral Taylor remainder after
degree 240 is bounded by 48^241/241! for all twelve half-unit waits.
This bound uses the Hermitian physical source, not the Euclidean norm of A.

Every polynomial term, renewal product, squared amplitude and conditional
ratio is enclosed by directed integer arithmetic on a fixed dyadic grid.
The exported endpoints are integers divided by 10^8. No floating-point
estimate is an assumption of these intervals. They describe conditional
simulation-generated records with their success and timeout probabilities.
"""

from fractions import Fraction
from math import factorial, isqrt
from numbers import Integral

from ..carrier import require

BITS = 80
DEGREE = 240
GRID = 10 ** 8
NORM_BOUND = 8


def _ceil_div(numerator, denominator):
    require(denominator > 0, 'interval division requires a positive denominator')
    return -(-numerator // denominator)


def _fraction(value):
    require(isinstance(value, (Integral, Fraction)) and not isinstance(value, bool),
            'a certificate parameter must be an integer or a Fraction')
    return Fraction(value)


def _add(left, right):
    return left[0] + right[0], left[1] + right[1]


def _sub(left, right):
    return left[0] - right[1], left[1] - right[0]


def _mul(left, right, scale):
    corners = [a * b for a in left for b in right]
    return min(corners) // scale, _ceil_div(max(corners), scale)


def _square(value, scale):
    low, high = value
    lower = 0 if low <= 0 <= high else min(low * low, high * high)
    return lower // scale, _ceil_div(max(low * low, high * high), scale)


def _complex_sub(left, right):
    return _sub(left[0], right[0]), _sub(left[1], right[1])


def _complex_mul(left, right, scale):
    real = _sub(_mul(left[0], right[0], scale), _mul(left[1], right[1], scale))
    imag = _add(_mul(left[0], right[1], scale), _mul(left[1], right[0], scale))
    return real, imag


def _probability(value, scale):
    low, high = _add(_square(value[0], scale), _square(value[1], scale))
    # A first-detection branch is a product of unitaries and projectors.
    return max(0, low), min(scale, high)


def _ratio(numerator, denominator, scale):
    require(0 <= numerator[0] <= numerator[1] and 0 < denominator[0] <= denominator[1],
            'conditional interval requires a positive success lower bound')
    return (numerator[0] * scale // denominator[1],
            _ceil_div(numerator[1] * scale, denominator[0]))


def _rational_scale(interval, coefficient):
    coefficient = _fraction(coefficient)
    require(coefficient >= 0, 'interval scaling requires a nonnegative coefficient')
    return (interval[0] * coefficient.numerator // coefficient.denominator,
            _ceil_div(interval[1] * coefficient.numerator, coefficient.denominator))


def _root_ratio_upper(numerator, denominator):
    """A rational upper bound for sqrt(numerator/denominator), proved by squaring."""
    scale = 1 << 24
    root = isqrt(numerator * scale * scale // denominator)
    if root * root * denominator < numerator * scale * scale:
        root += 1
    require(root * root * denominator >= numerator * scale * scale,
            'the rational square-root bound is not outward')
    return Fraction(root, scale)


def certify_chart(chart, weights, stations, norm_bound):
    """Check exact Hermiticity and an operator-norm bound in the physical metric."""
    rows = [list(row) for row in chart]
    size = len(rows)
    require(size > 0 and all(len(row) == size for row in rows), 'chart must be square')
    require(all(isinstance(x, Integral) and not isinstance(x, bool) for row in rows for x in row),
            'the source chart must have integer entries')
    rows = [[int(x) for x in row] for row in rows]
    weights = list(weights)
    require(len(weights) == size and all(isinstance(x, Integral) and not isinstance(x, bool)
                                        and x > 0 for x in weights), 'invalid integer similarity metric')
    weights = [int(x) for x in weights]
    stations = list(stations)
    require(stations and all(isinstance(x, Integral) and not isinstance(x, bool)
                            and 0 <= x < size for x in stations)
            and len(set(stations)) == len(stations), 'invalid station basis')
    stations = [int(x) for x in stations]
    require(all(weights[x] == 1 for x in stations), 'station vectors must be unchanged by the similarity')
    require(all(weights[i] * rows[i][j] == weights[j] * rows[j][i]
                for i in range(size) for j in range(size)), 'chart is not Hermitian in its declared metric')
    norm_bound = _fraction(norm_bound)
    require(norm_bound >= 0, 'negative source norm bound')
    root_bounds = {(a, b): _root_ratio_upper(a, b) for a in set(weights) for b in set(weights)}
    row_bounds = [sum((Fraction(abs(value), 2) * root_bounds[weights[i], weights[j]]
                       for j, value in enumerate(row)), Fraction(0)) for i, row in enumerate(rows)]
    require(max(row_bounds) <= norm_bound, 'physical row sums exceed the proposed operator-norm bound')
    return rows, weights, stations, max(row_bounds)


def integer_moments(chart, stations, degree):
    """Return (A^k)[target, source] for every station pair, using integer products."""
    size = len(chart)
    sparse_rows = [[(j, value) for j, value in enumerate(row) if value] for row in chart]
    current = [[int(i == station) for station in stations] for i in range(size)]
    moments = []
    for order in range(degree + 1):
        moments.append([current[station][:] for station in stations])
        if order != degree:
            current = [[sum(value * current[j][source] for j, value in row)
                        for source in range(len(stations))] for row in sparse_rows]
    return moments


def enclose_propagators(chart, weights, stations, ticks=12, dt=Fraction(1, 2),
                        degree=DEGREE, bits=BITS, norm_bound=NORM_BOUND):
    """Enclose every unmonitored station propagator in directed dyadic rectangles.

    The amplitude axes are tick, target, source, real/imaginary, lower/upper.
    Each endpoint is an integer divided by ``dyadic_scale``.
    """
    require(type(ticks) is int and ticks > 0, 'invalid number of detector ticks')
    require(type(degree) is int and degree >= 0, 'invalid Taylor degree')
    require(type(bits) is int and bits >= 16, 'invalid dyadic precision')
    dt = _fraction(dt)
    require(dt > 0, 'detector interval must be positive')
    norm_bound = _fraction(norm_bound)
    chart, weights, stations, row_upper = certify_chart(chart, weights, stations, norm_bound)
    moments = integer_moments(chart, stations, degree)
    size, count, scale = len(chart), len(stations), 1 << bits
    # Cayley-Hamilton makes a zero prefix of length size an exact all-time zero.
    identically_zero = [[degree >= size - 1 and all(moments[k][b][a] == 0 for k in range(size))
                         for a in range(count)] for b in range(count)]
    amplitudes = [[[( (scale * int(a == b), scale * int(a == b)), (0, 0))
                    for a in range(count)] for b in range(count)]]
    remainder_ceilings = [0]
    maximum_remainder = Fraction(0)
    for tick in range(1, ticks + 1):
        time = tick * dt
        remainder = (time * norm_bound) ** (degree + 1) / factorial(degree + 1)
        maximum_remainder = max(maximum_remainder, remainder)
        radius = _ceil_div(remainder.numerator * scale, remainder.denominator)
        remainder_ceilings.append(radius)
        rectangles = [[[0, 0, 0, 0] for _ in stations] for _ in stations]
        coefficient = Fraction(1)
        for order in range(degree + 1):
            component = 0 if order % 2 == 0 else 2
            sign = (1, -1, -1, 1)[order % 4]
            numerator = sign * coefficient.numerator * scale
            denominator = coefficient.denominator
            for target in range(count):
                for source in range(count):
                    term = numerator * moments[order][target][source]
                    rectangles[target][source][component] += term // denominator
                    rectangles[target][source][component + 1] += _ceil_div(term, denominator)
            coefficient *= time / (2 * (order + 1))
        enclosed = []
        for target in range(count):
            row = []
            for source in range(count):
                if identically_zero[target][source]:
                    row.append(((0, 0), (0, 0)))
                else:
                    rl, rh, il, ih = rectangles[target][source]
                    # Unitarity bounds each component by one on these unit vectors.
                    row.append(((max(-scale, rl - radius), min(scale, rh + radius)),
                                (max(-scale, il - radius), min(scale, ih + radius))))
            enclosed.append(row)
        amplitudes.append(enclosed)
    return dict(dyadic_scale=scale, bits=bits, degree=degree, ticks=ticks,
                dt_numerator=dt.numerator, dt_denominator=dt.denominator,
                station_indices=stations, amplitudes=amplitudes,
                physical_norm_upper=[norm_bound.numerator, norm_bound.denominator],
                physical_row_sum_upper=[row_upper.numerator, row_upper.denominator],
                maximum_taylor_remainder=[maximum_remainder.numerator, maximum_remainder.denominator],
                component_remainder_dyadic_ceilings=remainder_ceilings)


def renewal_enclosures(propagators):
    """Directed renewal products with all hit, timeout and undefined-mean outcomes."""
    kernels = propagators['amplitudes']
    scale, ticks = propagators['dyadic_scale'], propagators['ticks']
    dt = Fraction(propagators['dt_numerator'], propagators['dt_denominator'])
    count = len(propagators['station_indices'])
    zero = (0, 0)
    first = [[[zero for _ in range(ticks)] for _ in range(count)] for _ in range(count)]
    success = [[zero for _ in range(count)] for _ in range(count)]
    timeout = [[zero for _ in range(count)] for _ in range(count)]
    means = [[zero for _ in range(count)] for _ in range(count)]
    defined = [[False for _ in range(count)] for _ in range(count)]
    for source in range(count):
        for target in range(count):
            if source == target:
                continue
            paths = [(zero, zero)]
            for tick in range(1, ticks + 1):
                path = kernels[tick][target][source]
                for earlier in range(1, tick):
                    path = _complex_sub(path, _complex_mul(kernels[tick - earlier][target][target],
                                                          paths[earlier], scale))
                paths.append(path)
                first[source][target][tick - 1] = _probability(path, scale)
            probabilities = first[source][target]
            total = (sum(p[0] for p in probabilities), min(scale, sum(p[1] for p in probabilities)))
            require(0 <= total[0] <= total[1] <= scale, 'empty probability enclosure')
            success[source][target] = total
            timeout[source][target] = scale - total[1], scale - total[0]
            if total[0] > 0:
                numerator = (sum((n + 1) * p[0] for n, p in enumerate(probabilities)),
                             sum((n + 1) * p[1] for n, p in enumerate(probabilities)))
                low, high = _ratio(numerator, total, scale)
                means[source][target] = max(scale, low), min(ticks * scale, high)
                defined[source][target] = True
            else:
                # This is a conditional range if success occurs, not an asserted mean.
                means[source][target] = scale, ticks * scale
    distances = [[zero for _ in range(count)] for _ in range(count)]
    round_success = [[zero for _ in range(count)] for _ in range(count)]
    pair_defined = [[False for _ in range(count)] for _ in range(count)]
    for left in range(count):
        for right in range(left + 1, count):
            interval = _rational_scale(_add(means[left][right], means[right][left]), dt / 2)
            distances[left][right] = distances[right][left] = interval
            probability = _mul(success[left][right], success[right][left], scale)
            round_success[left][right] = round_success[right][left] = probability
            pair_defined[left][right] = pair_defined[right][left] = defined[left][right] and defined[right][left]
    return dict(dyadic_scale=scale, first_hit_probability=first, leg_success_probability=success,
                leg_timeout_probability=timeout, leg_conditional_mean_ticks=means,
                leg_positive_success_proved=defined, distance=distances,
                round_trip_success_probability=round_success, distance_defined=pair_defined)


def export_grid(enclosures, grid=GRID):
    """Return integer lower and upper arrays on a declared rational output grid."""
    require(type(grid) is int and grid > 0, 'invalid certificate output grid')
    scale = enclosures['dyadic_scale']
    def split(tree):
        if isinstance(tree, tuple) and len(tree) == 2:
            low, high = tree
            return low * grid // scale, _ceil_div(high * grid, scale)
        parts = [split(child) for child in tree]
        return [part[0] for part in parts], [part[1] for part in parts]
    result = dict(scale=grid, leg_positive_success_proved=enclosures['leg_positive_success_proved'],
                  distance_defined=enclosures['distance_defined'])
    for name in ('first_hit_probability', 'leg_success_probability', 'leg_timeout_probability',
                 'leg_conditional_mean_ticks', 'distance', 'round_trip_success_probability'):
        result[name + '_lower'], result[name + '_upper'] = split(enclosures[name])
    return result


def run(spec):
    """Certify both frozen source preparations; no sampled records enter this function."""
    from .contact_response import exact_chart
    from .joint_source import ENDPOINTS, MARK, validate_source

    validate_source(spec)
    require(tuple(spec['endpoints']) == ENDPOINTS and spec['mark'] == MARK,
            'certified source stations differ from the declared experiment')
    require(spec['dt'] == .5 and spec['ticks'] == 12, 'certified source requires twelve half-unit ticks')
    chart, weights = exact_chart(spec)
    settings = []
    for marked in (0, 1):
        current = chart.copy()
        current[spec['mark'], spec['mark']] += 2 * marked
        propagators = enclose_propagators(current, weights, spec['endpoints'])
        result = export_grid(renewal_enclosures(propagators))
        result.update(marked_contact=marked,
                      physical_row_sum_upper=propagators['physical_row_sum_upper'],
                      maximum_taylor_remainder=propagators['maximum_taylor_remainder'],
                      component_remainder_dyadic_ceilings=propagators['component_remainder_dyadic_ceilings'])
        settings.append(result)
    return dict(scale=GRID, dyadic_precision_bits=BITS, taylor_degree=DEGREE,
                physical_source_norm_bound=NORM_BOUND, source_dimensions=len(chart),
                station_indices=list(ENDPOINTS), marked_seam=MARK,
                dt_numerator=1, dt_denominator=2, ticks=12,
                arithmetic='integer moments and outward dyadic operations; exact rational output endpoints',
                remainder_rule='Hermitian integral Taylor bound (time*norm)^(degree+1)/(degree+1)!',
                similarity='A=2*S^-1*H*S; S=diag(I,sqrt(2)*I); station vectors have metric weight one',
                first_detection_rule='f_ba(n)=K_ba(n)-sum_(j=1)^(n-1) K_bb(n-j)*f_ba(j)',
                distance_rule='dt*(conditional_tick_mean_ab+conditional_tick_mean_ba)/2',
                diagonal_menu='unused zero sentinel',
                conditional_domain='a positive success lower bound is proved for each defined directed leg',
                settings=settings)
