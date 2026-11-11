"""Law C, the Z6 register and its integer lift on a small closed federation: exact reference implementation.

All state is Python integers. numpy draws the preparation and the schedule with the generator calls of the public
simulator, so a run of this module equals the public engine attempt for attempt.

Objects
    carrier     icosahedron with 12 ports, 30 seams and 20 oriented faces
    federation  N = 20 * 4^level carriers, three ports of each glued to a neighbour by a glue seam; every glue
                seam bounds two screen plaquettes
    law C       attempt on seam (a, b) with loads x_a, x_b and coin c: t = x_a + x_b, x_a becomes ceil(t/2) if
                c = 1 else floor(t/2), x_b = t - x_a; current j = old x_a - new x_a
    register    a repair with current j proposes delta = j mod 6; the face that traverses the seam a to b changes
                its class by +delta, the face that traverses it b to a by -delta; applied when
                cost[new] - cost[old] <= kappa * j * (d - j) with d = x_a - x_b
    lift        integer link n_s per seam, zero at the start; every applied twist adds j. Face value
                m_f = sum of n_s over the seams the face traverses a to b minus the sum over those it traverses
                b to a
"""
from __future__ import annotations

import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np

COST = (0, 5, 8, 10, 8, 5)
KAPPA = 5
LOAD_SEED_BASE = 20260909
SCHEDULE_SEED_BASE = 909000
GLUING_INDEX = 1
MINIMAL = {0: (0,), 1: (1,), 2: (2,), 3: (3, -3), 4: (-2,), 5: (-1,)}
DATA = json.loads(Path(__file__).with_name('federation.json').read_text())
TEMPLATE_SEAMS = [tuple(s) for s in DATA['template_seams']]
TEMPLATE_FACES = [tuple(f) for f in DATA['template_faces']]


def hypercharge_integers():
    """6Y = 3 n2 - 2 n3 on the 32 occupation patterns of three colour levels and two weak levels."""
    return sorted({3 * n2 - 2 * n3 for n3 in range(4) for n2 in range(3)})


def template_face_sides():
    """For each template seam (a, b), a < b: the face that traverses it a to b and the face that traverses it b to a."""
    key = {s: k for k, s in enumerate(TEMPLATE_SEAMS)}
    plus, minus = [None] * 30, [None] * 30
    for f, (p, q, r) in enumerate(TEMPLATE_FACES):
        for u, v in ((p, q), (q, r), (r, p)):
            k = key[(min(u, v), max(u, v))]
            side = plus if u < v else minus
            assert side[k] is None
            side[k] = f
    return plus, minus


class Federation:
    """Seam and face tables of the federation at one level (levels 0, 1 and 2 are stored)."""

    def __init__(self, level: int):
        fed = DATA[f'level_{level}']
        plus, minus = template_face_sides()
        self.level = level
        self.carriers = n = fed['carriers']
        self.plaquettes = fed['plaquettes']
        self.ports = 12 * n
        self.intra_seams = 30 * n
        self.intra_faces = 20 * n
        self.end_a, self.end_b, self.face_plus, self.face_minus = [], [], [], []
        for c in range(n):
            for k, (a, b) in enumerate(TEMPLATE_SEAMS):
                self.end_a.append(12 * c + a)
                self.end_b.append(12 * c + b)
                self.face_plus.append(20 * c + plus[k])
                self.face_minus.append(20 * c + minus[k])
        for a, b, p, m in zip(fed['glue_a'], fed['glue_b'], fed['glue_plus'], fed['glue_minus']):
            self.end_a.append(a)
            self.end_b.append(b)
            self.face_plus.append(20 * n + p)
            self.face_minus.append(20 * n + m)
        self.seams = len(self.end_a)
        self.faces = 20 * n + self.plaquettes

    def structure(self) -> dict:
        """Closed-surface facts used by the gates."""
        directed = Counter()
        for p, q, r in TEMPLATE_FACES:
            directed.update(((p, q), (q, r), (r, p)))
        degree = Counter(v for s in TEMPLATE_SEAMS for v in s)
        glued = Counter(self.end_a[self.intra_seams:] + self.end_b[self.intra_seams:])
        per_carrier = Counter(p // 12 for p in glued)
        around = Counter(self.face_plus[self.intra_seams:] + self.face_minus[self.intra_seams:])
        return dict(
            ports=12, seams=len(TEMPLATE_SEAMS), faces=len(TEMPLATE_FACES), euler_characteristic=12 - 30 + 20,
            every_directed_seam_once=len(directed) == 60 and set(directed.values()) == {1}
            and all((b, a) in directed for a, b in directed),
            port_degrees=sorted(set(degree.values())),
            every_seam_has_two_distinct_faces=all(p != m for p, m in zip(self.face_plus, self.face_minus)),
            internal_seam_faces_on_one_carrier=all(
                self.face_plus[s] // 20 == self.face_minus[s] // 20 == s // 30 for s in range(self.intra_seams)),
            glue_seam_faces_are_plaquettes=all(
                self.face_plus[s] >= self.intra_faces and self.face_minus[s] >= self.intra_faces
                for s in range(self.intra_seams, self.seams)),
            glue_seams=self.seams - self.intra_seams, ports_glued_once=set(glued.values()) == {1},
            glue_seams_per_carrier=sorted(set(per_carrier.values())),
            glue_seams_join_two_carriers=all(
                self.end_a[s] // 12 != self.end_b[s] // 12 for s in range(self.intra_seams, self.seams)),
            seams_around_a_plaquette=sorted(set(around.values())),
            screen_euler_characteristic=self.plaquettes - (self.seams - self.intra_seams) + self.carriers)


def attempt(xa: int, xb: int, coin: int) -> tuple[int, int, int]:
    """Law C on one seam: new loads of both ends and the current j from the first end to the second."""
    tot = xa + xb
    lo = tot // 2
    na = tot - lo if coin else lo
    return na, tot - na, xa - na


def accepts(k1: int, k2: int, j: int, d: int) -> bool:
    """Register rule: the twist delta = j mod 6 on faces of classes k1 (plus side) and k2 (minus side) is applied
    when its cost change is at most kappa * j * (d - j), with d the contrast before the repair."""
    delta = j % 6
    n1, n2 = (k1 + delta) % 6, (k2 - delta) % 6
    return COST[n1] + COST[n2] - COST[k1] - COST[k2] <= KAPPA * j * (d - j)


def prepare(fed: Federation) -> list[int]:
    """Preparation exact06 of the public simulator."""
    rng = np.random.default_rng(LOAD_SEED_BASE + fed.level)
    return rng.integers(0, 6, size=fed.ports).tolist()


def schedule(fed: Federation, sweeps: int, index: int = 1):
    """Seam and coin draws of every sweep, exactly as the public simulator draws them."""
    rng = np.random.default_rng(SCHEDULE_SEED_BASE + 1000 * fed.level + 100 * GLUING_INDEX + index)
    for _ in range(sweeps):
        seq = rng.integers(0, fed.seams, size=fed.seams, dtype=np.int64)
        coin = rng.integers(0, 2, size=fed.seams, dtype=np.int64)
        yield seq.tolist(), coin.tolist()


def run(fed: Federation, sweeps: int, index: int = 1, max_depth: int | None = None, after_sweep=None) -> dict:
    """Execute the schedule. With ``max_depth`` only the attempts of causal depth at most that value act: the result
    is the state on that consistent cut. Depth of an attempt = 1 + the largest depth of the earlier attempts on its
    resources (its two ports, and for a glue seam its two plaquettes)."""
    A, B, FP, FM = fed.end_a, fed.end_b, fed.face_plus, fed.face_minus
    SI, FI, P = fed.intra_seams, fed.intra_faces, fed.ports
    x = prepare(fed)
    x0 = list(x)
    reg = [0] * fed.seams            # the Z6 link of the run
    cls = [0] * fed.faces            # the Z6 class of the run
    link = [0] * fed.seams           # integer link of the lift
    val = [0] * fed.faces            # integer face value of the lift
    depth = [0] * (P + fed.plaquettes)
    twists = []                      # (sweep, position, seam, j, depth, value of plus face before, of minus face before)
    max_abs_j = 0
    for s, (seq, coin) in enumerate(schedule(fed, sweeps, index)):
        for i, sm in enumerate(seq):
            a, b = A[sm], B[sm]
            d = depth[a] if depth[a] > depth[b] else depth[b]
            if sm >= SI:
                ra, rb = P + FP[sm] - FI, P + FM[sm] - FI
                if depth[ra] > d:
                    d = depth[ra]
                if depth[rb] > d:
                    d = depth[rb]
                d += 1
                depth[ra] = depth[rb] = d
            else:
                d += 1
            depth[a] = depth[b] = d
            if max_depth is not None and d > max_depth:
                continue
            xi, xj = x[a], x[b]
            tot = xi + xj
            lo = tot // 2
            ni = tot - lo if coin[i] else lo
            j = xi - ni
            if not j:
                continue
            x[a], x[b] = ni, tot - ni
            delta = j % 6
            fp, fm = FP[sm], FM[sm]
            k1, k2 = cls[fp], cls[fm]
            n1, n2 = (k1 + delta) % 6, (k2 - delta) % 6
            if COST[n1] + COST[n2] - COST[k1] - COST[k2] <= KAPPA * j * (xi - xj - j):
                reg[sm] = (reg[sm] + delta) % 6
                cls[fp], cls[fm] = n1, n2
                twists.append((s, i, sm, j, d, val[fp], val[fm]))
                link[sm] += j
                val[fp] += j
                val[fm] -= j
                if abs(j) > max_abs_j:
                    max_abs_j = abs(j)
        if after_sweep is not None:
            after_sweep(s, x, reg, cls, link, val)
    return dict(loads=x, loads_initial=x0, register=reg, classes=cls, link=link, values=val, twists=twists,
                max_abs_j=max_abs_j, depth_max=max(depth[:P]), depth_min=min(depth[:P]))


def curl(fed: Federation, link: list[int]) -> list[int]:
    """Face values of an integer link field: the signed sum over the seams of each face."""
    out = [0] * fed.faces
    for s, n in enumerate(link):
        out[fed.face_plus[s]] += n
        out[fed.face_minus[s]] -= n
    return out


def gates(fed: Federation, state: dict) -> dict:
    """The exact gates of the lift on one state."""
    val, cls, link, reg = state['values'], state['classes'], state['link'], state['register']
    fi = fed.intra_faces
    return dict(
        value_mod_6_is_the_class=all(v % 6 == k for v, k in zip(val, cls)),
        link_mod_6_is_the_register=all(n % 6 == r for n, r in zip(link, reg)),
        value_is_the_curl_of_the_link=curl(fed, link) == val,
        every_carrier_total_is_zero=all(sum(val[20 * c:20 * c + 20]) == 0 for c in range(fed.carriers)),
        screen_total_is_zero=sum(val[fi:]) == 0,
        class_total_zero_mod_6_per_carrier=all(sum(cls[20 * c:20 * c + 20]) % 6 == 0 for c in range(fed.carriers)))


def fmt(v: int) -> str:
    return f'{v:+d}' if v else '0'


def word(values: list[int]) -> tuple[str, str]:
    """Composition and integer word of one carrier: defect values by class, then after a bar the hidden values."""
    defects = sorted((v % 6, -v) for v in values if v % 6)
    hidden = sorted((v for v in values if v and v % 6 == 0), reverse=True)
    comp = '{' + ','.join(str(k) for k, _ in defects) + '}'
    text = '(' + ','.join(fmt(-nv) for _, nv in defects) + ('|' + ','.join(fmt(h) for h in hidden) if hidden else '') + ')'
    return comp, text


def census(fed: Federation, state: dict) -> dict:
    """Histogram of the face values per class and tier, and the integer words of the carriers."""
    val = state['values']
    fi = fed.intra_faces
    hist = {'internal': Counter(), 'screen': Counter()}
    for f, v in enumerate(val):
        if v:
            hist['internal' if f < fi else 'screen'][v] += 1
    words = Counter()
    for c in range(fed.carriers):
        row = val[20 * c:20 * c + 20]
        if any(row):
            words[word(row)] += 1
    allowed = set(hypercharge_integers())

    def table(h):
        out = {str(k): {} for k in range(6)}
        for v in sorted(h):
            out[str(v % 6)][fmt(v)] = h[v]
        return out
    defects = [v for v in val if v % 6]
    return dict(
        values_by_class={t: table(h) for t, h in hist.items()},
        defect_faces=len(defects), hidden_faces=sum(1 for v in val if v and v % 6 == 0),
        defect_faces_at_a_hypercharge_integer=sum(1 for v in defects if v in allowed),
        defect_faces_minimal=sum(1 for v in defects if v in MINIMAL[v % 6]),
        words={f'{c} {w}': n for (c, w), n in sorted(words.items(), key=lambda kv: (-kv[1], kv[0]))})


def digest(seq) -> str:
    """sha256 of a list of small integers as signed bytes (the int8 arrays of the public engine)."""
    return hashlib.sha256(bytes(v & 0xFF for v in seq)).hexdigest()


def twist_digest(twists) -> str:
    """sha256 over (sweep, position, seam, j, depth) of every applied twist, as decimal text."""
    h = hashlib.sha256()
    for s, i, sm, j, d, _, _ in twists:
        h.update(f'{s},{i},{sm},{j},{d};'.encode())
    return h.hexdigest()
