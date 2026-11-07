"""Exact finite controls and bounds independent of global state enumeration.

The repeated weight word is inherited from the declared weighted extension.
Mean occupation two or four per cell is an explicit input. Bounds concern
diagonal canonical states and complete local experiments, not preparation,
physical spatial volume, or a source-selected field theory.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
from math import comb, isqrt

from ..weighted_source_family import family as owning


WORD = (25, 25)+(8,)*10


def integer(value, lower=0):
    if type(value) is not int or value < lower:
        raise ValueError('an exact integer in the stated domain is required')
    return value


def native_word():
    word = owning.local_family().weights
    if word != WORD:
        raise ArithmeticError('the retained native weight word differs')
    return word


def declared_mean(value):
    integer(value, 1)
    if value not in (2, 4):
        raise ValueError('the proved mean occupation is exactly two or four')
    return value


@lru_cache(maxsize=8, typed=True)
def parameters(bits=192, *, mean=4):
    """Rational enclosures for the unique positive root and occupations."""
    integer(bits, 8)
    declared_mean(mean)
    discriminant, shift, denominator = (6401, 1, 1600) if mean == 4 else (1256, -16, 1000)
    scale = 1 << bits
    root_floor = isqrt(discriminant*scale*scale)
    root = F(root_floor, scale), F(root_floor+1, scale)
    if not root[0]**2 < discriminant < root[1]**2:
        raise ArithmeticError('the algebraic root enclosure failed')
    z = tuple((shift+x)/denominator for x in root)
    high = tuple(25*x/(1+25*x) for x in z)
    low = tuple(8*x/(1+8*x) for x in z)
    hmin, hmax, lmin, lmax = ((F(1,2), F(3,5), F(1,4), F(1,3)) if mean == 4
                             else (F(3,10), F(1,3), F(1,8), F(1,6)))
    if not hmin < high[0] < high[1] < hmax:
        raise ArithmeticError('the high-weight occupation range failed')
    if not lmin < low[0] < low[1] < lmax:
        raise ArithmeticError('the low-weight occupation range failed')
    polynomial = (lambda t:800*t*t-t-2) if mean == 4 else (lambda t:1000*t*t+32*t-1)
    if not polynomial(z[0]) < 0 < polynomial(z[1]):
        raise ArithmeticError('the fugacity polynomial does not straddle zero')
    return {'root':root, 'fugacity':z, 'high':high, 'low':low}


def sqrt_upper(value, bits=256):
    """Outward rational square-root bound; no binary floating arithmetic."""
    if type(value) is not F or value < 0:
        raise ValueError('a nonnegative exact rational is required')
    integer(bits, 1)
    if not value:
        return F()
    shift = max(0, (value.denominator.bit_length()-value.numerator.bit_length()+1)//2)
    scale = 1 << (bits+shift)
    floor = isqrt(value.numerator*scale*scale//value.denominator)
    if floor*floor*value.denominator == value.numerator*scale*scale:
        return F(floor, scale)
    return F(floor+1, scale)


def log_upper(value):
    """ln(value) < (7/10)*bit_length(value), for positive integers."""
    integer(value, 1)
    if value == 1:
        return F()
    return F(7*value.bit_length(), 10)


def local_bound(cells, window_cells=1, *, mean=4):
    """All-domain entropy and half-trace bounds on any retained cell set."""
    integer(cells, 1); integer(window_cells)
    declared_mean(mean)
    if window_cells > cells:
        raise ValueError('the retained cell count exceeds the source')
    log_bound = log_upper(240*cells)
    entropy_bound = F(window_cells, cells)*log_bound
    squared = min(F(1), entropy_bound/2)
    return {'cells':cells, 'ports':12*cells, 'mean_per_cell':mean,
            'total_occupation':mean*cells,
            'window_cells':window_cells, 'window_ports_upper':12*window_cells,
            'conditioning_probability_lower':F(1, 240*cells),
            'global_relative_entropy_upper':log_bound,
            'local_relative_entropy_upper':entropy_bound,
            'half_trace_distance_squared_upper':squared,
            'half_trace_distance_upper':sqrt_upper(squared)}


def level_bound(level, window_cells=1, *, mean=4):
    """Use the retained cell-count grammar without constructing its graph."""
    integer(level)
    result = local_bound(20*(4**level), window_cells, mean=mean)
    return {'level':level, **result}


def two_size_bound(first_cells, second_cells, window_cells=1, *, mean=4):
    """A common local labeling compares both references to one product state."""
    first = local_bound(first_cells, window_cells, mean=mean)
    second = local_bound(second_cells, window_cells, mean=mean)
    return min(F(1), first['half_trace_distance_upper']+
               second['half_trace_distance_upper'])


def instrument_error(input_error, *, gain=F(1), implementation_error=F()):
    """Propagate through a certified CP trace gain, without normalization."""
    if (any(type(value) is not F for value in (input_error,gain,implementation_error))
            or not 0 <= input_error <= 1 or gain < 0 or implementation_error < 0):
        raise ValueError('exact nonnegative bounds and an input state distance at most one are required')
    return gain*input_error+implementation_error


@lru_cache(maxsize=256, typed=True)
def coefficient(high_ports, low_ports, occupied):
    """Finite type-count coefficient, used only for bounded examples."""
    integer(high_ports); integer(low_ports)
    if type(occupied) is not int:
        raise ValueError('an integer coefficient degree is required')
    if not 0 <= occupied <= high_ports+low_ports:
        return 0
    return sum(comb(high_ports, a)*25**a*comb(low_ports, occupied-a)*8**(occupied-a)
               for a in range(max(0, occupied-low_ports), min(high_ports, occupied)+1))


def canonical_counts(cells, window_cells, *, mean=4):
    """Exact full local state grouped by its two exchangeable weight types.

    Each count pair represents all equally weighted configurations of that
    pair. Zero sectors are retained. This finite-control evaluator is not
    used by local_bound or level_bound at large capacity.
    """
    local_bound(cells, window_cells, mean=mean)
    denominator = coefficient(2*cells, 10*cells, mean*cells)
    outside_high, outside_low = 2*(cells-window_cells), 10*(cells-window_cells)
    return {(a, b):F(comb(2*window_cells, a)*comb(10*window_cells, b)*25**a*8**b*
                      coefficient(outside_high, outside_low, mean*cells-a-b), denominator)
            for a in range(2*window_cells+1) for b in range(10*window_cells+1)}


def bernoulli_count_interval(ports, occupied, probability):
    """Outward interval for the complete binomial count probability."""
    integer(ports); integer(occupied)
    if (occupied > ports or type(probability) is not tuple or len(probability) != 2
            or any(type(x) is not F for x in probability)
            or not 0 < probability[0] <= probability[1] < 1):
        raise ValueError('a valid count and faithful rational probability interval are required')
    lo, hi = probability
    count = comb(ports, occupied)
    return (count*lo**occupied*(1-hi)**(ports-occupied),
            count*hi**occupied*(1-lo)**(ports-occupied))


def product_counts(window_cells, bits=192, *, mean=4):
    integer(window_cells)
    values = parameters(bits, mean=mean)
    high = [bernoulli_count_interval(2*window_cells, a, values['high'])
            for a in range(2*window_cells+1)]
    low = [bernoulli_count_interval(10*window_cells, b, values['low'])
           for b in range(10*window_cells+1)]
    return {(a, b):(h[0]*l[0], h[1]*l[1])
            for a, h in enumerate(high) for b, l in enumerate(low)}


def distance_interval(cells, window_cells=1, bits=192, *, mean=4):
    """Half trace distance for the entire local diagonal state, with rounding bounds."""
    canonical = canonical_counts(cells, window_cells, mean=mean)
    product = product_counts(window_cells, bits, mean=mean)
    lower = upper = F()
    for key, mass in canonical.items():
        lo, hi = product[key]
        lower += max(lo-mass, mass-hi, F())
        upper += max(abs(mass-lo), abs(mass-hi))
    return lower/2, min(F(1), upper/2)


def conditioning_interval(cells, bits=192, *, mean=4, occupation=None):
    """Exact-algebraic enclosure for a bounded finite conditioning example."""
    integer(cells, 1)
    declared_mean(mean)
    occupied = mean*cells if occupation is None else integer(occupation)
    if occupied > 12*cells:
        raise ValueError('occupation exceeds the full carrier')
    zlo, zhi = parameters(bits, mean=mean)['fugacity']
    numerator = coefficient(2*cells, 10*cells, occupied)
    return (numerator*zlo**occupied/((1+25*zhi)**(2*cells)*(1+8*zhi)**(10*cells)),
            numerator*zhi**occupied/((1+25*zlo)**(2*cells)*(1+8*zlo)**(10*cells)))


def adjacent_count_control(cells, bits=192, *, mean=4):
    """Finite exact check of the two counts used in the all-domain proof."""
    integer(cells, 1)
    values = parameters(bits, mean=mean)
    lo, hi = (2*cells*p for p in values['high'])
    a, second = lo.numerator//lo.denominator, hi.numerator//hi.denominator
    if a != second:
        raise ValueError('the algebraic interval does not isolate the required count')
    b = mean*cells-a
    high = bernoulli_count_interval(2*cells, a, values['high'])
    low = bernoulli_count_interval(10*cells, b, values['low'])
    return {'high_count':a, 'low_count':b,
            'event_probability_lower':high[0]*low[0],
            'theorem_conditioning_lower':F(1, 240*cells)}
