"""Exact cut and complete-record coherence bounds, with no fitted rate."""
from fractions import Fraction as F
from math import comb


def graph(vertices, edges, weights=None):
    if type(vertices) is not int or vertices < 2:
        raise ValueError('an exact vertex count of at least two is required')
    if (type(edges) is not tuple or not edges or len(set(edges)) != len(edges)
            or any(type(e) is not tuple or len(e) != 2 or any(type(x) is not int for x in e)
                   or not 0 <= e[0] < e[1] < vertices for e in edges)):
        raise ValueError('distinct sorted graph edges are required')
    if weights is None:
        weights = (F(1, len(edges)),)*len(edges)
    if (type(weights) is not tuple or len(weights) != len(edges)
            or any(type(w) not in (int, F) or w <= 0 for w in weights)
            or sum(weights) != 1):
        raise ValueError('strictly positive exact selector weights must sum to one')
    reached = {0}
    while True:
        new = reached | {b for a, b in edges if a in reached} | {a for a, b in edges if b in reached}
        if new == reached:
            break
        reached = new
    if len(reached) != vertices:
        raise ValueError('connected support is required for a positive cut bound')
    return tuple(map(F, weights))


def cut_table(vertices, edges, weights=None):
    weights = graph(vertices, edges, weights)
    if vertices > 20:
        raise ValueError('exhaustive cut enumeration is limited to twenty vertices')
    return tuple(sum((w for (a, b), w in zip(edges, weights, strict=True)
                      if ((mask >> a) ^ (mask >> b)) & 1), F())
                 for mask in range(1 << vertices))


def survival(first, second, vertices, edges, weights=None, *, law='C'):
    weights = graph(vertices, edges, weights)
    if (type(first) is not int or type(second) is not int
            or not 0 <= first < 1 << vertices or not 0 <= second < 1 << vertices
            or type(law) is not str or law not in ('C', 'E')):
        raise ValueError('two full occupation addresses and the C or E law are required')
    difference = [((first >> x) & 1)-((second >> x) & 1) for x in range(vertices)]
    killed = sum((w*(F(abs(difference[a]-difference[b]), 2) if law == 'C'
                     else int(difference[a] != difference[b]))
                  for (a, b), w in zip(edges, weights, strict=True)), F())
    return 1-killed


def coherence_upper(initial_l1, cuts, ticks, *, same_number=True, repair_probability=F(1)):
    if (type(initial_l1) not in (int, F) or initial_l1 < 0
            or type(ticks) is not int or ticks < 0 or type(same_number) is not bool
            or type(repair_probability) not in (int, F) or not 0 <= repair_probability <= 1
            or type(cuts) is not tuple or len(cuts) < 4 or len(cuts) & (len(cuts)-1)
            or any(type(v) not in (int, F) for v in cuts) or cuts[0] or cuts[-1]
            or any(not 0 < v <= 1 for v in cuts[1:-1])):
        raise ValueError('exact nonnegative norm, cut certificate and tick domain required')
    rate = min(cuts[1:-1])*F(repair_probability)/(1 if same_number else 2)
    return F(initial_l1)*(1-rate)**ticks


def mixing_upper(vertices, number, ticks, gap_lower, cut_lower):
    """Rational trace-norm bound; d replaces sqrt(d-1) conservatively."""
    if (any(type(v) is not int for v in (vertices, number, ticks))
            or vertices < 2 or not 0 < number < vertices or ticks < 0
            or type(gap_lower) not in (int, F) or not 0 < gap_lower <= 1
            or type(cut_lower) not in (int, F) or not 0 < cut_lower <= 1):
        raise ValueError('a nontrivial sector and proved positive rate bounds are required')
    dimension = comb(vertices, number)
    return min(F(2), dimension*(1-F(gap_lower))**ticks
               +(dimension-1)*(1-F(cut_lower))**ticks)


def polynomial_clock_limit(vertices, edges, ticks_per_time, weights=None):
    """The relevant rate is the actual cut per tick times ticks per unit time."""
    cuts = cut_table(vertices, edges, weights)
    if type(ticks_per_time) not in (int, F) or ticks_per_time <= 0:
        raise ValueError('a supplied positive exact clock conversion is required')
    return min(cuts[1:-1])*ticks_per_time


def diffusion_coherence_exponent(diffusion, spacing, duration):
    """Cubic one-particle coherence <= exp(-12 D t/a^2) under Poisson updates."""
    if (any(type(x) not in (int, F) for x in (diffusion, spacing, duration))
            or diffusion <= 0 or spacing <= 0 or duration < 0):
        raise ValueError('positive exact diffusion/spacing and nonnegative duration required')
    return 12*F(diffusion)*F(duration)/F(spacing)**2
