"""Exact arithmetic in Q(sqrt 5): quaternions, the binary icosahedral group, the carrier icosahedron."""
from fractions import Fraction as F
import itertools


class Q5:
    """a + b sqrt(5) with rational a, b."""
    __slots__ = ('a', 'b')

    def __init__(self, a=0, b=0):
        self.a, self.b = F(a), F(b)

    def __add__(self, o):
        o = _q5(o)
        return Q5(self.a + o.a, self.b + o.b)

    __radd__ = __add__

    def __neg__(self):
        return Q5(-self.a, -self.b)

    def __sub__(self, o):
        return self + (-_q5(o))

    def __rsub__(self, o):
        return _q5(o) - self

    def __mul__(self, o):
        o = _q5(o)
        return Q5(self.a * o.a + 5 * self.b * o.b, self.a * o.b + self.b * o.a)

    __rmul__ = __mul__

    def __eq__(self, o):
        o = _q5(o)
        return self.a == o.a and self.b == o.b

    def __hash__(self):
        return hash((self.a, self.b))

    def galois(self):
        """sqrt 5 -> -sqrt 5, the nontrivial automorphism of Q(sqrt 5)."""
        return Q5(self.a, -self.b)

    def sign(self):
        """Exact sign of a + b sqrt 5."""
        a, b = self.a, self.b
        if b == 0:
            return (a > 0) - (a < 0)
        if a == 0 or (a > 0) == (b > 0):
            return 1 if (a > 0 or (a == 0 and b > 0)) else -1
        # opposite signs: compare a^2 with 5 b^2
        dominant = a if a * a > 5 * b * b else b
        return 1 if dominant > 0 else -1

    def __repr__(self):
        return f'{self.a}+{self.b}r5'


def _q5(x):
    return x if isinstance(x, Q5) else Q5(x)


ZERO, ONE = Q5(0), Q5(1)
PHI = Q5(F(1, 2), F(1, 2))           # golden ratio
PHI_INV = Q5(F(-1, 2), F(1, 2))      # 1 / phi = phi - 1
HALF = Q5(F(1, 2))
assert PHI * PHI_INV == ONE


# ------------------------------------------------------------------ quaternions (w, x, y, z)
def qmul(p, q):
    a1, b1, c1, d1 = p
    a2, b2, c2, d2 = q
    return (a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2,
            a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2,
            a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2,
            a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2)


def qconj(q):
    return (q[0], -q[1], -q[2], -q[3])


def qnorm(q):
    return sum((c * c for c in q), ZERO)


def qgalois(q):
    return tuple(c.galois() for c in q)


QONE = (ONE, ZERO, ZERO, ZERO)
QMINUS = (-ONE, ZERO, ZERO, ZERO)
# generators of order 6 and 10: (1 + i + j + k) / 2 and (phi + phi^-1 i + j) / 2
GENERATORS = ((HALF, HALF, HALF, HALF), (PHI * HALF, PHI_INV * HALF, HALF, ZERO))


def binary_icosahedral():
    """Closure of the two generators: the 120 unit quaternions of 2I."""
    group = {QONE}
    frontier = [QONE]
    while frontier:
        new = []
        for g in frontier:
            for h in GENERATORS:
                gh = qmul(g, h)
                if gh not in group:
                    group.add(gh)
                    new.append(gh)
        frontier = new
    return sorted(group, key=lambda q: tuple((c.a, c.b) for c in q))


def rotation(q):
    """3 x 3 matrix of v -> q v q^-1 on the imaginary quaternions."""
    qc = qconj(q)
    cols = []
    for e in ((ZERO, ONE, ZERO, ZERO), (ZERO, ZERO, ONE, ZERO), (ZERO, ZERO, ZERO, ONE)):
        cols.append(qmul(qmul(q, e), qc)[1:])
    return tuple(tuple(cols[c][r] for c in range(3)) for r in range(3))


def det3(m):
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))


def apply3(m, v):
    return tuple(sum((m[r][c] * v[c] for c in range(3)), ZERO) for r in range(3))


def is_identity3(m):
    return all(m[r][c] == (ONE if r == c else ZERO) for r in range(3) for c in range(3))


# ------------------------------------------------------------------ the carrier: an icosahedron
def icosahedron():
    """Ports (vertices), seams (edges) and faces with outward orientation."""
    verts = []
    for cyc in range(3):
        for s1, s2 in itertools.product((1, -1), repeat=2):
            v = [ZERO, Q5(s1), PHI * s2]
            verts.append(tuple(v[(i - cyc) % 3] for i in range(3)))

    def d2(u, v):
        return sum(((x - y) * (x - y) for x, y in zip(u, v)), ZERO)

    edges = [(i, j) for i, j in itertools.combinations(range(12), 2) if d2(verts[i], verts[j]) == Q5(4)]
    adj = {i: set() for i in range(12)}
    for i, j in edges:
        adj[i].add(j)
        adj[j].add(i)
    faces = []
    for a, b, c in itertools.combinations(range(12), 3):
        if b in adj[a] and c in adj[a] and c in adj[b]:
            # outward orientation: det[v_a, v_b, v_c] > 0
            if det3((verts[a], verts[b], verts[c])).sign() < 0:
                b, c = c, b
            faces.append((a, b, c))
    return verts, edges, faces
