"""Candidate evidence from Fourier modes and executed signed gate words."""

from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import itertools
import json
import math

import mpmath as mp
import sympy as sp

from . import numerics


CONFIGS = (("axis3", 3, "axis"), ("axis4", 4, "axis"),
           ("ball6", 6, "ball"), ("two_scale6", 6, "two_scale"))
TICKS = (F(1, 256), F(1, 512))
TIMES = tuple(range(1, 33))


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                    allow_nan=False).encode()).hexdigest()


def stencil(kind):
    if kind == "ball":
        return sorted(v for v in itertools.product(range(-2, 3), repeat=3)
                      if 0 < sum(x*x for x in v) <= 4)
    radii = (1,) if kind == "axis" else (1, 2)
    return sorted(tuple(s*k if j == i else 0 for j in range(3))
                  for i in range(3) for k in radii for s in (-1, 1))


def spectrum(q, vectors):
    tables = {3: (F(1), F(-1, 2), F(-1, 2)),
              4: (F(1), F(0), F(-1), F(0)),
              6: (F(1), F(1, 2), F(-1, 2), F(-1), F(-1, 2), F(1, 2))}
    alpha = F(6*q*q, sum(sum(t*t for t in v) for v in vectors))
    counts = Counter()
    for k in itertools.product(range(q), repeat=3):
        value = 1+alpha*sum(1-tables[q][sum(a*b for a, b in zip(k, v)) % q]
                            for v in vectors)
        counts[value] += 1
    return sorted(counts.items())


@lru_cache(None)
def algebra():
    z = sp.Symbol("z")
    r = 2**sp.Rational(1, 3)
    w = (4+2*r+r*r)/6
    v = -r*w
    def step(t):
        return sp.Matrix([[1-t*t/2, t], [-t+t**3/4, 1-t*t/2]])
    result = {}
    for order, matrix in ((2, step(z)), (4, step(w*z)*step(v*z)*step(w*z))):
        rows = {}
        for label, i, j in (("a", 0, 0), ("b", 0, 1), ("c", 1, 0), ("d", 1, 1)):
            coeffs = []
            for n in range(8):
                expr = sp.to_number_field(sp.expand(matrix[i, j]).coeff(z, n), r)
                raw = list(expr.native_coeffs())
                raw = [sp.Rational(0)]*(3-len(raw))+raw
                coeffs.append([str(sp.Rational(x)) for x in reversed(raw)])
            rows[label] = coeffs
        result[str(order)] = rows
    return result


def number(x):
    return mp.nstr(x, 23, strip_zeros=False)


def mp_fraction(x):
    return mp.mpf(x.numerator)/x.denominator


def mode_step(z, order):
    def s(t):
        return mp.matrix([[1-t*t/2, t], [-t+t**3/4, 1-t*t/2]])
    if order == 2:
        return s(z)
    w = 1/(2-mp.root(2, 3))
    return s(w*z)*s(-mp.root(2, 3)*w*z)*s(w*z)


def quantum_case(spec, tick, order):
    if tuple(TIMES) != numerics.OBSERVATION_TIMES:
        raise ValueError("all 32 declared observation times are required")
    return numerics.quantum_observations(spec, tick, order, _mode_numbers)


def _mode_numbers(ctx, spec, tick, order):
    """Factored symplectic defect times the Chebyshev power polynomial.

    For S=[[a,b],[c,a]], det(S)=1 gives
    n_j=(b+c)^2 U_(j-1)(a)^2/4. The off-diagonal defect is factored
    before numerical evaluation, retaining arbitrarily small rational ticks.
    """
    if order == 4:
        r = ctx.mpf(2)**(ctx.mpf(1)/3)
        a6 = ctx.mpf(1)/48+5*r/288+r*r/72
        defect5 = ctx.mpf(1)/36+r/48+r*r/72
        defect7 = ctx.mpf(25)/1728+5*r/432+r*r/108
    for eigenvalue, multiplicity in spec:
        omega = ctx.sqrt(numerics.rational(ctx, eigenvalue))
        squared = numerics.rational(ctx, tick*tick*eigenvalue)
        if order == 2:
            a = 1-squared/2
            one_number = squared**3/64
        else:
            a = 1-squared/2+squared**2/24+a6*squared**3
            one_number = squared**5*(defect5+defect7*squared)**2/4
        previous, current = ctx.mpf(0), ctx.mpf(1)
        populations = []
        for _ in TIMES:
            populations.append(one_number*current**2)
            previous, current = current, 2*a*current-previous
        yield omega, multiplicity, populations


def gate_word(order, steps):
    r = mp.root(2, 3)
    weights = [mp.mpf(1)] if order == 2 else [1/(2-r), -r/(2-r), 1/(2-r)]
    word = []
    for _ in range(steps):
        for weight in weights:
            for kind, coefficient in (("K", weight/2), ("D", weight), ("K", weight/2)):
                if word and word[-1][0] == kind:
                    word[-1] = (kind, word[-1][1]+coefficient)
                else:
                    word.append((kind, coefficient))
    return word


def program(order, steps=3):
    """Retain each scalar writer and every read of a 27-site reference program."""
    with mp.workdps(65):
        sites = list(itertools.product(range(3), repeat=3))
        index = {v: i for i, v in enumerate(sites)}
        vectors = stencil("axis")
        neighbors = [[index[tuple((x[k]+v[k]) % 3 for k in range(3))]
                      for v in vectors] for x in sites]
        alpha = mp.mpf(9)
        q = [mp.mpf(1+x+2*y-3*z)/7 for x, y, z in sites]
        p = [mp.mpf(x-2*y+z)/11 for x, y, z in sites]
        q_version = p_version = 0
        records = []
        word = gate_word(order, steps)
        for stage, (kind, coefficient) in enumerate(word, 1):
            if kind == "K":
                new = [p[i]-coefficient/64*((1+6*alpha)*q[i]
                       -alpha*sum(q[k] for k in neighbors[i])) for i in range(27)]
                for i in range(27):
                    reads = [["P", p_version, i], ["Q", q_version, i]]
                    reads += [["Q", q_version, k] for k in neighbors[i]]
                    records.append([stage, kind, i, reads, number(new[i])])
                p, p_version = new, stage
            else:
                new = [q[i]+coefficient/64*p[i] for i in range(27)]
                for i in range(27):
                    records.append([stage, kind, i,
                                    [["Q", q_version, i], ["P", p_version, i]], number(new[i])])
                q, q_version = new, stage
        kicks = sum(kind == "K" for kind, _ in word)
        final_q, final_p = [number(x) for x in q], [number(x) for x in p]
        for kind, coefficient in reversed(word):
            if kind == "K":
                p = [p[i]+coefficient/64*((1+6*alpha)*q[i]
                     -alpha*sum(q[k] for k in neighbors[i])) for i in range(27)]
            else:
                q = [q[i]-coefficient/64*p[i] for i in range(27)]
        expected_q = [mp.mpf(1+x+2*y-3*z)/7 for x, y, z in sites]
        expected_p = [mp.mpf(x-2*y+z)/11 for x, y, z in sites]
        if max(abs(a-b) for a, b in zip(q+p, expected_q+expected_p)) >= mp.mpf("1e-55"):
            raise ValueError("inverse execution exceeded the declared residual bound")
        return dict(steps=steps, sites=27, kicks=kicks, drifts=len(word)-kicks,
                    signed_word=[[k, number(c)] for k, c in word], writes=len(records),
                    nonlocal_reads=kicks*27*6, local_reads=len(records)*2,
                    two_mode_phase_gates=kicks*27*3, onsite_phase_gates=kicks*27,
                    onsite_kinetic_gates=(len(word)-kicks)*27, inverse_writes=len(records),
                    inverse_max_error_bound="1e-55",
                    trace_sha256=digest(records),
                    final_q=final_q, final_p=final_p)


def ball_moments(radius):
    count = moment = 0
    for x in range(-radius, radius+1):
        for y in range(-radius, radius+1):
            remaining = radius*radius-x*x-y*y
            if remaining < 0:
                continue
            z = math.isqrt(remaining)
            count += 2*z+1
            moment += (2*z+1)*(x*x+y*y)+z*(z+1)*(2*z+1)//3
    return count-1, moment


def sparse_first_level():
    q, m, k, r = 2**16, 2**8, 2**6, 2**7
    dc, mc = ball_moments(k)
    df, mf = ball_moments(r)
    degree, moment = dc+df, m*m*mc+mf
    weight = F(6*q*q*degree, moment)
    tick = F(1, 64*2**10)
    return dict(q=q, m=m, K=k, r=r, coarse_degree=dc, fine_degree=df,
                degree=degree, second_moment=moment, outgoing_weight=str(weight),
                repaired_tick=str(tick), stability_square_bound=str(tick*tick*(1+2*weight)))


def exponents():
    result = {}
    for family, degree, sigma in (("sparse", F(21, 16), F(7, 32)),
                                   ("dense", F(3, 2), F(1, 2))):
        result[family] = {}
        for order in (2, 4):
            vacuum = F(3, 2*order)+sigma
            energy = F(3, 2*order)+(1+F(1, 2*order))*sigma
            result[family][str(order)] = dict(vacuum_tick_power=str(vacuum),
                energy_tick_power=str(energy), fixed_tolerance_read_power=str(3+degree+energy))
    result["explicit_repair"] = dict(tick_power="5/8", read_power="79/16",
        vacuum_error_power=str(3+8*F(7, 32)-8*F(5, 8)),
        energy_error_power=str(3+9*F(7, 32)-8*F(5, 8)),
        second_order_energy_power=str(3+5*F(7, 32)-4*F(5, 8)))
    return result


def extreme_signals():
    """Execute impulses on a non-wrapping ray control, retaining all sites."""
    result = {}
    with mp.workdps(90):
        for order in (2, 4):
            for steps in (1, 2, 3):
                q = [mp.mpf(int(i == 0)) for i in range(32)]
                p = [mp.mpf(0) for _ in range(32)]
                for kind, coefficient in gate_word(order, steps):
                    if kind == "K":
                        p = [p[i]-coefficient/64*(3*q[i]-q[(i-1) % 32]-q[(i+1) % 32])
                             for i in range(32)]
                    else:
                        q = [q[i]+coefficient/64*p[i] for i in range(32)]
                q_radius = max(min(i, 32-i) for i, x in enumerate(q) if x != 0)
                p_radius = max(min(i, 32-i) for i, x in enumerate(p) if x != 0)
                result[f"{order}:{steps}"] = dict(sites=32, tick="1/64", max_q_distance=q_radius,
                    max_p_distance=p_radius, q=[number(x) for x in q], p=[number(x) for x in p],
                    outermost_signal=number(p[p_radius]))
    return result


def candidate():
    graphs = {}
    for name, q, kind in CONFIGS:
        vectors = stencil(kind)
        spec = spectrum(q, vectors)
        graphs[name] = dict(q=q, kind=kind, vectors=[list(v) for v in vectors],
            sites=q**3, degree=len(vectors), spectrum=[[str(v), n] for v, n in spec],
            moments={str(j): str(sum(n*v**j for v, n in spec)) for j in range(1, 6)},
            quantum={f"{p}:{dt}": quantum_case(spec, dt, p) for p in (2, 4) for dt in TICKS})
    return dict(algebra=algebra(), graphs=graphs,
                programs={str(p): program(p) for p in (2, 4)},
                first_sparse_level=sparse_first_level(), exponents=exponents(),
                extreme_signals=extreme_signals())
