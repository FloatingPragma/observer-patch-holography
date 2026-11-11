"""The class of a face defect selects the pair on one port (section 7).

A class-k flux through the hexagon of the face f is the link phase w^k, w = exp(i pi/3), on the hexagon seam
(2, 3): the completion at both addresses of that seam is conjugated by diag(w^(k n_3)), n_3 the load on port 3.
Consensus states are the fixed vectors of the discriminant D = A^dag Q A, the uniform-address states fixed by the
tick. Arithmetic is exact in Q(w), with w^2 = w - 1.
"""
from fractions import Fraction as F
from itertools import combinations_with_replacement

from . import geometry, hexagon as hx, walk as wk


class Eis:
    """a + b w with rational a, b and w = exp(i pi/3)."""
    __slots__ = ('a', 'b')

    def __init__(self, a, b=0):
        self.a, self.b = F(a), F(b)

    @staticmethod
    def of(x):
        return x if isinstance(x, Eis) else Eis(x)

    def __add__(self, o):
        o = Eis.of(o)
        return Eis(self.a + o.a, self.b + o.b)

    __radd__ = __add__

    def __neg__(self):
        return Eis(-self.a, -self.b)

    def __sub__(self, o):
        return self + (-Eis.of(o))

    def __rsub__(self, o):
        return Eis.of(o) - self

    def __mul__(self, o):
        o = Eis.of(o)
        return Eis(self.a * o.a - self.b * o.b, self.a * o.b + self.b * o.a + self.b * o.b)

    __rmul__ = __mul__

    def conjugate(self):
        return Eis(self.a + self.b, -self.b)

    def norm(self):
        """|a + b w|^2 = a^2 + a b + b^2."""
        return self.a * self.a + self.a * self.b + self.b * self.b

    def __truediv__(self, o):
        o = Eis.of(o)
        c = self * o.conjugate()
        n = o.norm()
        return Eis(c.a / n, c.b / n)

    def __bool__(self):
        return bool(self.a) or bool(self.b)

    def __eq__(self, o):
        o = Eis.of(o)
        return self.a == o.a and self.b == o.b

    def __hash__(self):
        return hash((self.a, self.b))

    def __repr__(self):
        return f'{self.a}+{self.b}w'


def w_power(k):
    out = Eis(1)
    for _ in range(k % 6):
        out = out * Eis(0, 1)
    return out


FLUX_SEAM, FLUX_PORT = (2, 3), 3
CONFIGS = tuple(sorted(wk.config(a, b) for a, b in combinations_with_replacement(geometry.hexagon(), 2)))
STATES = tuple((c, r) for c in CONFIGS for r in range(4))
RING = frozenset(frozenset(x) for x in hx.RING)
E_TYPE = frozenset(frozenset(x) for x in geometry.summary()['e_type_seams'])


def repair_with_flux(state, k, pairs=hx.ADDRESSES):
    """The plain law on the hexagon with the class-k link phase on the seam (2, 3)."""
    phase = w_power(k)
    out = {}
    for (cfg, rec, kk), amp in state.items():
        a, b = pairs[kk]
        on_flux = {a, b} == set(FLUX_SEAM)
        g_in = phase.conjugate() if (on_flux and cfg[FLUX_PORT] == 1) else Eis(1)
        if on_flux and cfg[FLUX_PORT] == 2:
            g_in = (phase * phase).conjugate()
        for c2, r2, u in wk.images(cfg, rec, a, b):
            g_out = Eis(1)
            if on_flux:
                g_out = w_power(k * c2[FLUX_PORT])
            wk.add(out, (c2, r2, kk), Eis.of(amp) * g_in * g_out * u)
    return out


def discriminant(k):
    """D[(c2, r2), (c, r)] = (1/12) sum over addresses of the phased completion."""
    D = {}
    for cfg, rec in STATES:
        col = {}
        for kk in range(len(hx.ADDRESSES)):
            for (c2, r2, _), v in repair_with_flux({(cfg, rec, kk): Eis(1)}, k).items():
                col[(c2, r2)] = col.get((c2, r2), Eis(0)) + v * F(1, len(hx.ADDRESSES))
        D[(cfg, rec)] = {key: v for key, v in col.items() if v}
    return D


def sectors():
    separated = [(c, r) for c, r in STATES if max(c) == 1 and r in (0, 1)]
    pair = [(c, r) for c, r in STATES if (max(c) == 2 and r == 0)
            or (max(c) == 1 and r in (2, 3) and frozenset(wk.ports_of(c)) in RING)]
    return separated, pair


def kernel_dimension_and_vector(D, sector):
    """Dimension of the fixed space of D on the sector, and one fixed vector (exact Gaussian elimination)."""
    index = {s: i for i, s in enumerate(sector)}
    n = len(sector)
    rows = [[Eis(0)] * n for _ in range(n)]
    for j, s_ in enumerate(sector):
        for key, v in D[s_].items():
            if key in index:
                rows[index[key]][j] = rows[index[key]][j] + v
    for i in range(n):
        rows[i][i] = rows[i][i] - 1
    pivots, r = [], 0
    for c in range(n):
        p = next((i for i in range(r, n) if rows[i][c]), None)
        if p is None:
            continue
        rows[r], rows[p] = rows[p], rows[r]
        inv = Eis(1) / rows[r][c]
        rows[r] = [x * inv for x in rows[r]]
        for i in range(n):
            if i != r and rows[i][c]:
                f = rows[i][c]
                rows[i] = [x - f * y for x, y in zip(rows[i], rows[r])]
        pivots.append(c)
        r += 1
    free = [c for c in range(n) if c not in pivots]
    vector = None
    if free:
        v = [Eis(0)] * n
        v[free[0]] = Eis(1)
        for i, c in enumerate(pivots):
            v[c] = -rows[i][free[0]]
        vector = {sector[i]: x for i, x in enumerate(v) if x}
    return len(free), vector


def closed_and_hermitian(D, sector):
    inside = set(sector)
    closed = all(key in inside for s_ in sector for key in D[s_])
    hermitian = all(D.get(key, {}).get(s_, Eis(0)) == v.conjugate() for s_ in D for key, v in D[s_].items())
    return closed, hermitian


def emission_class(key):
    """Flags one tick after blank doubles: in-dimer across an E-type seam; forward or backward across a C-type seam
    by whether the emitting port lies in the third-turn orbit {2, 9, 0}."""
    cfg, rec, kk = key
    if rec == 0 and max(cfg) == 2:
        return 'double'
    if rec not in (2, 3):
        return 'other'
    a, b = hx.ADDRESSES[kk]
    seam = frozenset((a, b))
    if seam not in RING or set(wk.ports_of(cfg)) != set(seam):
        return 'other'
    emitter = a if rec == 2 else b
    if seam in E_TYPE:
        return 'in-dimer'
    return 'forward' if emitter in hx.ORBIT else 'backward'


def result():
    separated, pair = sectors()
    table, vectors, structure = {}, {}, True
    for k in range(6):
        D = discriminant(k)
        for sector in (separated, pair):
            closed, hermitian = closed_and_hermitian(D, sector)
            structure &= closed and hermitian
        n_sep, _ = kernel_dimension_and_vector(D, separated)
        n_pair, v = kernel_dimension_and_vector(D, pair)
        table[k] = (n_sep, n_pair)
        vectors[k] = v
    v3 = vectors[3]
    norms = {x.norm() for x in v3.values()}
    blank = {key: x for key, x in v3.items() if key[1] == 0}
    blank_weight = sum((x.norm() for x in blank.values()), F(0)) / sum((x.norm() for x in v3.values()), F(0))
    uniform_state = {(cfg, rec, kk): x for (cfg, rec), x in v3.items() for kk in range(len(hx.ADDRESSES))}
    walk = wk.Walk(hx.ADDRESSES)
    ticked = repair_with_flux(walk.reflect(uniform_state), 3)
    fixed_by_tick = ticked == uniform_state
    u6 = {(wk.config(p, p), 0, kk): Eis(1) for p in geometry.hexagon() for kk in range(len(hx.ADDRESSES))}
    t6 = repair_with_flux(walk.reflect(u6), 3)
    n2 = F(len(u6))
    weights = {}
    for key, amp in t6.items():
        weights[emission_class(key)] = weights.get(emission_class(key), F(0)) + amp.norm() / n2
    back = walk.reflect(repair_with_flux(t6, 3)) == u6
    p, q = weights.get('forward', F(0)), weights.get('in-dimer', F(0))
    expected = {0: (2, 1), 1: (0, 0), 2: (0, 0), 3: (0, 1), 4: (0, 0), 5: (0, 0)}
    checks = dict(
        sectors_closed_and_discriminant_hermitian_for_every_class=structure,
        consensus_states_by_class_0_2_1_class_3_0_1_others_none=table == expected,
        class_3_consensus_state_is_uniform_over_18_states=len(v3) == 18 and len(norms) == 1,
        its_blank_record_part_is_the_uniform_double_with_weight_1_3=(
            blank_weight == F(1, 3) and len(blank) == 6 and len(set(blank.values())) == 1
            and all(max(cfg) == 2 for cfg, _ in blank)),
        class_3_consensus_state_is_fixed_by_the_tick=fixed_by_tick,
        uniform_double_tick_weights_2_3_1_6_1_12_1_12_count_2_9=(
            weights == {'double': F(2, 3), 'in-dimer': F(1, 6), 'forward': F(1, 12), 'backward': F(1, 12)}
            and p + q - 2 * p * q == F(2, 9) and back))
    data = dict(
        flux_seam=list(FLUX_SEAM), states={'separated loads, records 0 and 1': len(separated),
                                           'doubles and their emission states': len(pair)},
        consensus_states={str(k): {'separated loads': a, 'pair on one port': b} for k, (a, b) in table.items()},
        class_3_blank_record_weight=str(blank_weight),
        uniform_double_tick_weights={k: str(v) for k, v in sorted(weights.items())},
        uniform_double_count=str(p + q - 2 * p * q))
    return data, checks
