"""Private arithmetic contexts and bounded interval evaluation of observations.

The acceptance check uses directed mpmath interval arithmetic, rather than
agreement between successive point approximations. It is a numerical enclosure
calculation, not a formal proof of the underlying transcendental library.
"""

from fractions import Fraction

from mpmath.ctx_iv import MPIntervalContext
from mpmath.ctx_mp import MPContext
from mpmath.libmp import finf, fnan, fninf


OUTPUT_DIGITS = 23
PRECISION_STEPS = (90, 180, 360, 720, 1440, 2880, 5760)
OBSERVATION_TIMES = tuple(range(1, 33))


def _positive_rational(value, name):
    if type(value) is not int and not isinstance(value, Fraction):
        raise ValueError(f"{name} must be a positive int or Fraction")
    value = Fraction(value)
    if value <= 0:
        raise ValueError(f"{name} must be positive")
    return value


def validate_inputs(spectrum, tick, order):
    """Retain exact supplied scalars before checking the certified domain."""
    if type(order) is not int or order not in (2, 4):
        raise ValueError("order must be the integer 2 or 4")
    tick = _positive_rational(tick, "tick")
    try:
        entries = list(spectrum)
    except TypeError as exc:
        raise ValueError("spectrum must be a nonempty sequence of mode pairs") from exc
    if not entries:
        raise ValueError("spectrum must be nonempty")
    result = []
    for entry in entries:
        try:
            eigenvalue, multiplicity = entry
        except (TypeError, ValueError) as exc:
            raise ValueError("each spectrum entry must be an eigenvalue/count pair") from exc
        eigenvalue = _positive_rational(eigenvalue, "eigenvalue")
        if type(multiplicity) is not int or multiplicity <= 0:
            raise ValueError("multiplicity must be a positive integer")
        if tick*tick*eigenvalue > Fraction(1, 100):
            raise ValueError("outside certified domain: tick squared times eigenvalue exceeds 1/100")
        result.append((eigenvalue, multiplicity))
    return tuple(result), tick, order


def rational(ctx, value):
    """Enclose an original rational with directed rounding on both endpoints."""
    value = Fraction(value)
    return ctx.mpf(value.numerator)/ctx.mpf(value.denominator)


def _finite(value):
    return all(endpoint not in (finf, fninf, fnan) for endpoint in value._mpi_)


def log1p_positive(ctx, value):
    """Enclose log(1+x) for nonnegative x without forming a rounded 1+x.

    For y=x/(2+x), the positive atanh series has remainder at most
    2*y**(2*k+1)/((2*k+1)*(1-y*y)) after k terms. Every operation and the
    addition of this remainder use outward interval rounding.
    """
    if not _finite(value) or value.a < 0:
        raise ValueError("log1p requires a finite nonnegative interval")
    if value.b == 0:
        return ctx.mpf(0)
    if value.b > ctx.mpf(1)/8:
        # This branch has no small-signal cancellation; interval addition
        # still encloses 1+x when the input interval is broad.
        return ctx.log(1+value)
    y = value/(2+value)
    y_squared = y*y
    power = y
    partial = 2*y
    count = 1
    while True:
        power *= y_squared
        remainder = 2*power/((2*count+1)*(1-y_squared))
        target = partial.b*ctx.mpf(2)**(-ctx.prec)
        if remainder.b <= target.a:
            return partial+ctx.mpf([0, remainder.b])
        partial += 2*power/(2*count+1)
        count += 1


def _formatted_interval(value, formatter):
    if not _finite(value) or not value.a > 0:
        return None
    # make_mpf retains each binary endpoint exactly, without a precision-
    # changing conversion through float or the caller's global context.
    endpoints = [formatter.nstr(formatter.make_mpf(endpoint), OUTPUT_DIGITS,
                                strip_zeros=False) for endpoint in value._mpi_]
    return endpoints[0] if endpoints[0] == endpoints[1] else None


def quantum_observations(spectrum, tick, order, mode_numbers):
    """Evaluate a mode-number callback afresh until decimal output is resolved.

    The callback receives (private_interval_context, spectrum, tick, order)
    and yields (omega_interval, multiplicity, 32 particle-number intervals).
    Original exact inputs are reused for every precision attempt.
    """
    spectrum, tick, order = validate_inputs(spectrum, tick, order)
    formatter = MPContext()
    for precision in PRECISION_STEPS:
        ctx = MPIntervalContext()
        ctx.dps = precision
        rows = [[ctx.mpf(0) for _ in range(3)] for _ in OBSERVATION_TIMES]
        unresolved = False
        for omega, multiplicity, numbers in mode_numbers(ctx, spectrum, tick, order):
            numbers = tuple(numbers)
            if len(numbers) != len(OBSERVATION_TIMES):
                raise ValueError("mode callback must retain all 32 observation times")
            if not _finite(omega) or not omega.a > 0:
                unresolved = True
                break
            for row, number in zip(rows, numbers):
                # A positive-domain mode cannot be declared stationary
                # merely because its computed enclosure includes zero.
                if not _finite(number) or not number.a > 0:
                    unresolved = True
                    break
                row[0] += multiplicity*number
                row[1] += multiplicity*omega*number
                row[2] += multiplicity*log1p_positive(ctx, number)/2
            if unresolved:
                break
        if unresolved:
            continue
        means = [sum((row[i] for row in rows), ctx.mpf(0))/len(rows) for i in range(3)]
        encoded_rows = [[_formatted_interval(value, formatter) for value in row] for row in rows]
        encoded_means = [_formatted_interval(value, formatter) for value in means]
        if all(value is not None for row in encoded_rows for value in row) and all(
                value is not None for value in encoded_means):
            return dict(times=list(OBSERVATION_TIMES), observations=encoded_rows, mean=encoded_means)
    raise ValueError("vacuum observations unresolved within the interval precision budget")
