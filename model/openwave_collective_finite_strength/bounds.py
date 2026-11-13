"""Outward rational elementary functions and complete spin-sector inertia."""
from fractions import Fraction as F
from math import factorial, isqrt

GRID = 10**100


def down(x):
    x = F(x)
    return F((x*GRID).__floor__(), GRID)


def up(x):
    return -down(-F(x))


def root_bounds(x):
    x = F(x)
    if x < 0:
        raise ValueError('nonnegative radicand required')
    k = isqrt(x.numerator*GRID**2//x.denominator)
    lo = F(k, GRID)
    return lo, lo if lo*lo == x else F(k+1, GRID)


def exp_bounds(x):
    x = F(x)
    if abs(x) > 1:
        raise ValueError('the exponential certificate requires |x|<=1')
    term = total = F(1)
    for k in range(1, 81):
        term *= x/k
        total += term
    tail = abs(x)**81/(factorial(81)*(1-abs(x)/82))
    return down(total-tail), up(total+tail)


def cos_bounds(x):
    x = F(x)
    if abs(x) > 8:
        raise ValueError('the cosine certificate requires |x|<=8')
    total = sum(((-1)**k*x**(2*k)/factorial(2*k) for k in range(81)), F(0))
    tail = abs(x)**162/(factorial(162)*(1-x*x/(163*164)))
    return max(F(-1), down(total-tail)), min(F(1), up(total+tail))


def moments(n, q, m2):
    if type(n) is not int or n < 1 or not 0 <= m2 <= 1:
        raise ValueError('complete positive register count and squared spin mean required')
    q, m2 = F(q), F(m2)
    mu2 = q*q*(F(1,n)+(1-F(1,n))*m2)
    mu4 = q**4/F(n**4)*(n*(3*n-2)+n*(n-1)*(6*n-8)*m2
                          +n*(n-1)*(n-2)*(n-3)*m2*m2)
    return mu2, mu4


def jacobi_diagonal(n, twice_spin, terms, q=F(1,8), omega=F(1)):
    if type(n) is not int or n < 1 or not 0 <= twice_spin <= n or (n-twice_spin)%2:
        raise ValueError('an actual complete spin sector is required')
    if not 0 < q <= F(1,2) or omega <= 0:
        raise ValueError('positive frequency and q<=1/2 required')
    lo, hi = [], []
    for k in range(twice_spin+1):
        twice_m = 2*k-twice_spin
        a = b = omega*n/2
        for f, energy in terms:
            elo, ehi = map(F, energy)
            if not 0 <= f <= 1 or not 0 <= elo <= ehi:
                raise ValueError('positive native energies and profiles in [0,1] required')
            e0, e1 = exp_bounds(F(f)*q*twice_m/n)
            a += elo*e0; b += ehi*e1
        lo.append(down(a)); hi.append(up(b))
    squares = [omega*omega*F((k+1)*(twice_spin-k),4) for k in range(twice_spin)]
    return lo, hi, squares


def inertia_below(diagonal, squares, x):
    """Sturm sign variations; exact zero entries are removed, including the endpoint."""
    x = F(x); previous, current = F(1), diagonal[0]-x
    signs = [1]
    if current: signs.append(1 if current > 0 else -1)
    for d, b2 in zip(diagonal[1:], squares, strict=True):
        previous, current = current, (d-x)*current-b2*previous
        if current: signs.append(1 if current > 0 else -1)
    return sum(a != b for a,b in zip(signs, signs[1:]))


def bracket_eigenvalue(diagonal, squares, index, width=F(1,2**38)):
    if not 0 <= index < len(diagonal):
        raise ValueError('eigenvalue index outside complete block')
    # Two times the largest off-diagonal bound encloses every Gershgorin row.
    radius = 2*max((root_bounds(x)[1] for x in squares), default=F(0))+1
    lo, hi = min(diagonal)-radius, max(diagonal)+radius
    while hi-lo > width:
        middle = (lo+hi)/2
        if inertia_below(diagonal, squares, middle) <= index:
            lo = middle
        else:
            hi = middle
    return lo, hi


def block_eigenvalue(n, twice_spin, terms, index, **kwargs):
    lo, hi, squares = jacobi_diagonal(n, twice_spin, terms, **kwargs)
    return bracket_eigenvalue(lo,squares,index)[0], bracket_eigenvalue(hi,squares,index)[1]


def complete_gap(n, terms, **kwargs):
    ground = block_eigenvalue(n,n,terms,0,**kwargs)
    candidates = [block_eigenvalue(n,n,terms,1,**kwargs)]
    if n > 1: candidates.append(block_eigenvalue(n,n-2,terms,0,**kwargs))
    excited = (min(x[0] for x in candidates),min(x[1] for x in candidates))
    return {'ground':ground,'first_excited':excited,
            'gap':(excited[0]-ground[1],excited[1]-ground[0]),
            'candidates':candidates}
