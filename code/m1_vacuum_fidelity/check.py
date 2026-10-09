"""Independent exact algebra, spatial matrices and analytic Gaussian replay.

This module never imports the candidate producer.
"""

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import itertools
import json
import math

import mpmath as mp
import numpy as np

from . import numerics


@dataclass(frozen=True)
class Cubic:
    """Q[r]/(r^3-2), with rational arithmetic throughout."""
    coefficients: tuple = (F(0), F(0), F(0))

    @staticmethod
    def lift(value):
        return value if isinstance(value, Cubic) else Cubic((F(value), F(0), F(0)))

    def __add__(self, other):
        other = self.lift(other)
        return Cubic(tuple(a+b for a, b in zip(self.coefficients, other.coefficients)))

    __radd__ = __add__

    def __neg__(self):
        return Cubic(tuple(-a for a in self.coefficients))

    def __sub__(self, other):
        return self+-self.lift(other)

    def __rsub__(self, other):
        return self.lift(other)+-self

    def __mul__(self, other):
        other = self.lift(other)
        result = [F(0)]*5
        for i, a in enumerate(self.coefficients):
            for j, b in enumerate(other.coefficients):
                result[i+j] += a*b
        for n in (4, 3):
            result[n-3] += 2*result[n]
        return Cubic(tuple(result[:3]))

    __rmul__ = __mul__

    def __truediv__(self, denominator):
        return self*F(1, denominator)

    def __pow__(self, exponent):
        if type(exponent) is not int or exponent < 0:
            raise ValueError("nonnegative integer field exponent required")
        result = self.lift(1)
        for _ in range(exponent):
            result *= self
        return result

    def real(self):
        r = mp.root(2, 3)
        return sum(mp.mpf(c.numerator)/c.denominator*r**i
                   for i, c in enumerate(self.coefficients))

    def interval(self):
        lo, hi = F(1259, 1000), F(126, 100)
        lower = upper = F(0)
        for i, c in enumerate(self.coefficients):
            lower += c*(lo**i if c >= 0 else hi**i)
            upper += c*(hi**i if c >= 0 else lo**i)
        return lower, upper


ZERO, ONE, ROOT = Cubic(), Cubic.lift(1), Cubic((F(0), F(1), F(0)))
W = (4+2*ROOT+ROOT*ROOT)/6
V = -ROOT*W


def need(condition, message):
    if not condition:
        raise ValueError(message)


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def same(left, right, message):
    need(canonical(left) == canonical(right), message)


def digest(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def decimal(value):
    return mp.nstr(value, 23, strip_zeros=False)


def padd(a, b):
    return [(a[i] if i < len(a) else ZERO)+(b[i] if i < len(b) else ZERO)
            for i in range(max(len(a), len(b)))]


def pmul(a, b):
    out = [ZERO]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out


def mmul(a, b):
    return [[padd(pmul(a[i][0], b[0][j]), pmul(a[i][1], b[1][j]))
             for j in range(2)] for i in range(2)]


def polynomial_step(weight):
    return [[[ONE, ZERO, -weight**2/2], [ZERO, weight]],
            [[ZERO, -weight, ZERO, weight**3/4], [ONE, ZERO, -weight**2/2]]]


@lru_cache(None)
def exact_algebra():
    need(F(1259, 1000)**3 < 2 < F(126, 100)**3, "invalid cubic root enclosure")
    need(2*W+V == ONE and 2*W**3+V**3 == ZERO, "order cancellation failed")
    need(V.interval()[1] < 0, "inverse stage must be retained")
    answer = {}
    for order in (2, 4):
        matrix = polynomial_step(ONE)
        if order == 4:
            matrix = mmul(mmul(polynomial_step(W), polynomial_step(V)), polynomial_step(W))
        determinant = padd(pmul(matrix[0][0], matrix[1][1]),
                           [-x for x in pmul(matrix[0][1], matrix[1][0])])
        need(determinant[0] == ONE and all(x == ZERO for x in determinant[1:]),
             "symplectic identity failed")
        rows = {}
        for key, i, j in (("a", 0, 0), ("b", 0, 1), ("c", 1, 0), ("d", 1, 1)):
            need(all(x == ZERO for x in matrix[i][j][8:]), "unreported high coefficient")
            coeffs = (matrix[i][j]+[ZERO]*8)[:8]
            rows[key] = [[str(x) for x in c.coefficients] for c in coeffs]
        same(rows["a"], rows["d"], "reversibility diagonal mismatch")
        answer[str(order)] = rows
    certify_stability(answer)
    return answer


def certify_stability(algebra):
    """Bind rational interval bounds to the actual multiplied gate polynomials."""
    # Rational bounds prove a single domain, not just sampled eigenvalue stability.
    a6 = ONE/48+5*ROOT/288+ROOT**2/72
    b5 = ONE/36+ROOT/36+ROOT**2/48
    c5 = (ROOT+ROOT**2)/144
    d7 = 25*ONE/1728+5*ROOT/432+ROOT**2/108
    defect = ONE/36+ROOT/48+ROOT**2/72
    expected = {}
    for order, a, b, c in (
        ("2", [ONE, ZERO, -ONE/2], [ZERO, ONE], [ZERO, -ONE, ZERO, ONE/4]),
        ("4", [ONE, ZERO, -ONE/2, ZERO, ONE/24, ZERO, a6],
         [ZERO, ONE, ZERO, -ONE/6, ZERO, -b5],
         [ZERO, -ONE, ZERO, ONE/6, ZERO, c5, ZERO, -d7]),
    ):
        expected[order] = {
            name: [[str(x) for x in coefficient.coefficients]
                   for coefficient in (row+[ZERO]*8)[:8]]
            for name, row in (("a", a), ("b", b), ("c", c), ("d", a))
        }
    same(algebra, expected, "stability coefficients differ from exact gate multiplication")
    need(b5-c5 == defect, "stability defect is not the actual off-diagonal sum")
    need(0 < a6.interval()[0] < a6.interval()[1] < F(1, 10), "a6 enclosure")
    s = F(1, 100)
    need(1-s/2 > 0 and 1-s/4 > F(9, 10), "second-order stability domain")
    need(F(1, 24)*s+a6.interval()[1]*s*s < F(1, 2), "a<1 bound")
    b_lo = 1-s/6-b5.interval()[1]*s*s
    c_lo = 1-s/6-c5.interval()[1]*s*s
    c_hi = 1+d7.interval()[1]*s**3
    need(b_lo > F(99, 100) and c_lo > F(99, 100) and c_hi < F(101, 100),
         "off-diagonal sign bound")
    need(b_lo*c_lo > F(9, 10) and c_hi < F(11, 10), "phase bound")
    need(defect.interval()[0] > F(1, 16)
         and defect.interval()[1]+d7.interval()[1]*s < F(1, 8), "nonzero defect bound")


def vectors_for(kind):
    out = []
    for x, y, z in itertools.product(range(-2, 3), repeat=3):
        s = x*x+y*y+z*z
        if ((kind == "axis" and s == 1)
            or (kind == "ball" and 0 < s <= 4)
            or (kind == "two_scale" and s in (1, 4) and sum(v != 0 for v in (x, y, z)) == 1)):
            out.append((x, y, z))
    return out


def spatial(q, kind):
    vectors = vectors_for(kind)
    sites = list(itertools.product(range(q), repeat=3))
    links = [[sum(((x[i]+v[i]) % q)*q**(2-i) for i in range(3)) for v in vectors]
             for x in sites]
    need(all(len(set(row)) == len(vectors) and i not in row for i, row in enumerate(links)),
         "wrapped self-read or duplicate edge")
    alpha = F(6*q*q, sum(sum(x*x for x in v) for v in vectors))
    numerator, denominator = alpha.numerator, alpha.denominator
    diagonal = denominator+numerator*len(vectors)
    b = np.eye(q**3)*float(F(diagonal, denominator))
    for i, row in enumerate(links):
        b[i, row] -= float(alpha)
    # Closed-walk moments are exact integer matrix applications, not Fourier sums.
    vector = [0]*(q**3)
    vector[0] = 1
    moments = {}
    for power in range(1, 6):
        nxt = [diagonal*x for x in vector]
        for i, row in enumerate(links):
            for j in row:
                nxt[j] -= numerator*vector[i]
        vector = nxt
        moments[str(power)] = str(F(q**3*vector[0], denominator**power))
    twice_cos = {3: (2, -1, -1), 4: (2, 0, -2, 0), 6: (2, 1, -1, -2, -1, 1)}[q]
    histogram = Counter()
    for k in sites:
        twice_symbol = sum(2-twice_cos[sum(v[i]*k[i] for i in range(3)) % q] for v in vectors)
        histogram[1+alpha*F(twice_symbol, 2)] += 1
    spectrum = sorted(histogram.items())
    for power in range(1, 6):
        need(sum(n*l**power for l, n in spectrum) == F(moments[str(power)]),
             "spatial closed-walk and Fourier moment disagreement")
    eigenvalues, eigenvectors = np.linalg.eigh(b)
    expanded = np.array([float(l) for l, n in spectrum for _ in range(n)])
    need(np.allclose(eigenvalues, expanded, rtol=1e-12, atol=1e-11), "full spatial spectrum differs")
    return vectors, spectrum, moments, b, eigenvalues, eigenvectors


def _matrix_mode_numbers(ctx, spectrum, tick, order):
    """Enclose actual matrix powers from the independently multiplied gates.

    Cancellation in the off-diagonal sum is retained inside directed intervals;
    the reporting driver increases precision when that sum is unresolved. No
    producer defect factorization or Chebyshev recurrence is used here.
    """
    root = ctx.mpf(2)**(ctx.mpf(1)/3)
    powers = (ctx.mpf(1), root, root*root)
    coefficients = {
        key: [sum((numerics.rational(ctx, F(value))*power
                   for value, power in zip(triple, powers)), ctx.mpf(0))
              for triple in row]
        for key, row in exact_algebra()[str(order)].items()
    }

    def polynomial(row, z):
        value = ctx.mpf(0)
        for coefficient in reversed(row):
            value = value*z+coefficient
        return value

    def multiply(left, right):
        return [[left[i][0]*right[0][j]+left[i][1]*right[1][j]
                 for j in range(2)] for i in range(2)]

    for value, multiplicity in spectrum:
        omega = ctx.sqrt(numerics.rational(ctx, value))
        z = ctx.sqrt(numerics.rational(ctx, tick*tick*value))
        a, b, c, d = (polynomial(coefficients[key], z) for key in ("a", "b", "c", "d"))
        step = [[a, b], [c, d]]
        state = [[ctx.mpf(1), ctx.mpf(0)], [ctx.mpf(0), ctx.mpf(1)]]
        numbers = []
        for _ in range(32):
            state = multiply(step, state)
            numbers.append(((state[0][0]-state[1][1])**2
                            +(state[0][1]+state[1][0])**2)/4)
        yield omega, multiplicity, numbers


def quantum(spectrum, tick, order):
    return numerics.quantum_observations(spectrum, tick, order, _matrix_mode_numbers)


def word(order, steps):
    # Construct the already merged macro word directly, independently of S2 expansion.
    if order == 2:
        macro = [("K", ONE/2), ("D", ONE), ("K", ONE/2)]
    else:
        macro = [("K", W/2), ("D", W), ("K", (W+V)/2), ("D", V),
                 ("K", (W+V)/2), ("D", W), ("K", W/2)]
    result = []
    for i in range(steps):
        if i:
            result[-1] = ("K", result[-1][1]+macro[0][1])
            result.extend(macro[1:])
        else:
            result.extend(macro)
    return result


def spatial_quantum(b, eigenvalues, basis, order, reference):
    """Apply every kick/drift to the full site-space canonical matrix."""
    n = len(b)
    q_rows = np.hstack((np.eye(n), np.zeros((n, n))))
    p_rows = np.hstack((np.zeros((n, n)), np.eye(n)))
    omega = np.sqrt(eigenvalues)
    root = np.sqrt(omega)
    for j in range(1, 5):
        for kind, coefficient in word(order, 1):
            s = float(coefficient.real())/256
            if kind == "K":
                p_rows -= s*(b@q_rows)
            else:
                q_rows += s*p_rows
        if j not in (1, 2, 4):
            continue
        a = (basis.T@q_rows[:, :n]@basis)*(root[:, None]/root[None, :])
        d = (basis.T@p_rows[:, n:]@basis)*(root[None, :]/root[:, None])
        bb = (basis.T@q_rows[:, n:]@basis)*(root[:, None]*root[None, :])
        c = (basis.T@p_rows[:, :n]@basis)/(root[:, None]*root[None, :])
        beta = ((a-d)+1j*(bb+c))/2
        population = np.sum(np.abs(beta)**2, axis=1)
        singular_squared = np.linalg.eigvalsh(beta@beta.conj().T)
        result = [population.sum(), np.dot(omega, population),
                  np.log1p(np.maximum(singular_squared, 0)).sum()/2]
        target = [float(x) for x in reference["observations"][j-1]]
        need(np.allclose(result, target, rtol=2e-3, atol=1e-24),
             f"full site covariance disagrees for order {order}, step {j}: {result} != {target}")


def program(order):
    with mp.workdps(75):
        sites = list(itertools.product(range(3), repeat=3))
        links = [[sum(((x[i]+v[i]) % 3)*3**(2-i) for i in range(3))
                  for v in vectors_for("axis")] for x in sites]
        q = [Cubic.lift(F(1+x+2*y-3*z, 7)) for x, y, z in sites]
        p = [Cubic.lift(F(x-2*y+z, 11)) for x, y, z in sites]
        q_version = p_version = 0
        records = []
        instructions = word(order, 3)
        for stage, (kind, coefficient) in enumerate(instructions, 1):
            if kind == "K":
                out = []
                for i in range(27):
                    out.append(p[i]-coefficient/64*(55*q[i]-9*sum((q[k] for k in links[i]), ZERO)))
                    reads = [["P", p_version, i], ["Q", q_version, i]]
                    reads += [["Q", q_version, k] for k in links[i]]
                    records.append([stage, kind, i, reads, decimal(out[-1].real())])
                p, p_version = out, stage
            else:
                out = []
                for i in range(27):
                    out.append(q[i]+coefficient/64*p[i])
                    records.append([stage, kind, i, [["Q", q_version, i], ["P", p_version, i]],
                                    decimal(out[-1].real())])
                q, q_version = out, stage
        kick_count = 4 if order == 2 else 10
        drift_count = 3 if order == 2 else 9
        need(len(instructions) == kick_count+drift_count, "gate count mismatch")
        final_q, final_p = [decimal(x.real()) for x in q], [decimal(x.real()) for x in p]
        for kind, coefficient in reversed(instructions):
            if kind == "K":
                p = [p[i]+coefficient/64*(55*q[i]-9*sum((q[k] for k in links[i]), ZERO))
                     for i in range(27)]
            else:
                q = [q[i]-coefficient/64*p[i] for i in range(27)]
        need(q == [Cubic.lift(F(1+x+2*y-3*z, 7)) for x, y, z in sites], "inverse positions fail")
        need(p == [Cubic.lift(F(x-2*y+z, 11)) for x, y, z in sites], "inverse momenta fail")
        return dict(steps=3, sites=27, kicks=kick_count, drifts=drift_count,
            signed_word=[[k, decimal(c.real())] for k, c in instructions],
            writes=27*(kick_count+drift_count), nonlocal_reads=kick_count*27*6,
            local_reads=54*(kick_count+drift_count), trace_sha256=digest(records),
            two_mode_phase_gates=kick_count*27*3, onsite_phase_gates=kick_count*27,
            onsite_kinetic_gates=drift_count*27, inverse_writes=27*(kick_count+drift_count),
            inverse_max_error_bound="1e-55",
            final_q=final_q, final_p=final_p)


def positive_octant_ball(radius):
    # Different enumeration: only nonnegative coordinates, with signed multiplicities.
    degree = moment = 0
    for x in range(radius+1):
        for y in range(math.isqrt(radius*radius-x*x)+1):
            last = math.isqrt(radius*radius-x*x-y*y)
            multiplier = (1 if x == 0 else 2)*(1 if y == 0 else 2)
            degree += multiplier*(2*last+1)
            moment += multiplier*((x*x+y*y)*(2*last+1)+last*(last+1)*(2*last+1)//3)
    return degree-1, moment


def first_level():
    coarse, mc = positive_octant_ball(64)
    fine, mf = positive_octant_ball(128)
    moment = 65536*mc+mf
    weight = F(6*65536**2*(coarse+fine), moment)
    tick = F(1, 65536)
    need(tick*tick*(1+2*weight) < F(1, 100), "explicit repair stability failed")
    return dict(q=65536, m=256, K=64, r=128, coarse_degree=coarse, fine_degree=fine,
        degree=coarse+fine, second_moment=moment, outgoing_weight=str(weight),
        repaired_tick=str(tick), stability_square_bound=str(tick*tick*(1+2*weight)))


def exponents():
    # Independent rational solutions of the number/energy balance equations.
    out = {}
    for family, degree, variance in (("sparse", F(21, 16), F(7, 16)),
                                      ("dense", F(3, 2), F(1))):
        out[family] = {}
        for p in (2, 4):
            vacuum = (3+p*variance)/(2*p)
            energy = (3+F(2*p+1, 2)*variance)/(2*p)
            out[family][str(p)] = dict(vacuum_tick_power=str(vacuum), energy_tick_power=str(energy),
                                      fixed_tolerance_read_power=str(3+degree+energy))
    out["explicit_repair"] = dict(tick_power="5/8", read_power=str(F(69, 16)+F(5, 8)),
        vacuum_error_power=str(F(96+56-160, 32)), energy_error_power=str(F(96+63-160, 32)),
        second_order_energy_power=str(F(96+35-80, 32)))
    return out


def extreme_signals():
    out = {}
    a6 = ONE/48+5*ROOT/288+ROOT**2/72
    d7 = 25*ONE/1728+5*ROOT/432+ROOT**2/108
    with mp.workdps(100):
        for order in (2, 4):
            for steps in (1, 2, 3):
                q = [ONE if i == 0 else ZERO for i in range(32)]
                p = [ZERO]*32
                for kind, coefficient in word(order, steps):
                    if kind == "K":
                        p = [p[i]-coefficient/64*(3*q[i]-q[(i-1) % 32]-q[(i+1) % 32])
                             for i in range(32)]
                    else:
                        q = [q[i]+coefficient/64*p[i] for i in range(32)]
                q_radius = steps if order == 2 else 3*steps
                p_radius = q_radius+1
                predicted = (Cubic.lift(F(1, 4*64**(2*steps+1))) if order == 2 else
                             -d7*(2*a6)**(steps-1)*F((-1)**(3*steps+1), 64**(6*steps+1)))
                need(p[p_radius] == predicted and predicted != ZERO, "extremal signal coefficient")
                need(q[q_radius] != ZERO, "position support must be attained")
                need(all(x == ZERO for i, x in enumerate(q) if min(i, 32-i) > q_radius),
                     "position influence outside exact support")
                need(all(x == ZERO for i, x in enumerate(p) if min(i, 32-i) > p_radius),
                     "momentum influence outside exact support")
                out[f"{order}:{steps}"] = dict(sites=32, tick="1/64", max_q_distance=q_radius,
                    max_p_distance=p_radius, q=[decimal(x.real()) for x in q],
                    p=[decimal(x.real()) for x in p], outermost_signal=decimal(predicted.real()))
    return out


@lru_cache(None)
def reconstruct():
    graphs = {}
    for name, q, kind in (("axis3", 3, "axis"), ("axis4", 4, "axis"),
                          ("ball6", 6, "ball"), ("two_scale6", 6, "two_scale")):
        vectors, spectrum, moments, b, eigenvalues, basis = spatial(q, kind)
        cases = {f"{p}:{dt}": quantum(spectrum, dt, p) for p in (2, 4)
                 for dt in (F(1, 256), F(1, 512))}
        for order in (2, 4):
            spatial_quantum(b, eigenvalues, basis, order, cases[f"{order}:1/256"])
        graphs[name] = dict(q=q, kind=kind, vectors=[list(v) for v in vectors], sites=q**3,
            degree=len(vectors), spectrum=[[str(v), n] for v, n in spectrum], moments=moments, quantum=cases)
    return dict(algebra=exact_algebra(), graphs=graphs, programs={str(p): program(p) for p in (2, 4)},
                first_sparse_level=first_level(), exponents=exponents(), extreme_signals=extreme_signals())


def verify_evidence(candidate):
    same(candidate, reconstruct(), "evidence differs from complete independent reconstruction")
