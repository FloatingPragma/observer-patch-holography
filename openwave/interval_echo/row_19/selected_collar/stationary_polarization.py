"""Read the stationary ordinary material left after a serial record is erased."""

from fractions import Fraction as F
from itertools import combinations, product
from math import gcd, lcm

import networkx as nx
import numpy as np
import sympy as sp

from . import material_polarization as motion

exp = motion.exp
local = exp.local
FLAGS = list(product(range(2), repeat=2))
SEED = 28122914
STEP = F(1, 32)


def symmetric_instrument(board):
    routes = exp.flight.routes
    steps, _, ports = routes.layout(board)
    match = routes.receiver_codes(board, steps)
    opened = set(local.survey.routes.cell().seams) - set(
        local.survey.routes.codebook()[0][0]
    )
    graph = nx.Graph(list(ports))
    hops = nx.single_source_shortest_path_length(graph, 0)

    def part(c, g):
        net = nx.Graph(sorted(opened))
        nx.set_node_attributes(
            net, {v: (v == g, c != steps[c, v]) for v in net}, "role"
        )
        return net

    for (u, v), a in sorted(ports.items()):
        if (hops[u], hops[v]) != (2, 3):
            continue
        b = ports[v, u]
        maps = nx.isomorphism.GraphMatcher(
            part(u, a % 12), part(v, b % 12), node_match=lambda x, y: x == y
        )
        if not maps.is_isomorphic():
            continue
        flip = min(tuple(f[j] for j in range(12)) for f in maps.isomorphisms_iter())
        match[v] = [tuple(sorted((flip[x], flip[y]))) for x, y in match[u]]

        def partner(c, g):
            return 12 * c + next(
                x + y - g % 12 for x, y in match[c] if g % 12 in (x, y)
            )

        meta = {
            "held_edge": [u, v],
            "probe_ports": [a, b],
            "ordinary_partners": [partner(u, a), partner(v, b)],
        }
        return exp.flight.native_instrument(board, meta, match), list(flip)
    raise ValueError("No exchange-symmetric native crossing in the enrolled shells")


def parity(board, edge):
    _, hist, _ = exp.flight.routes.layout(board)
    p = sp.Matrix(hist) / 12
    alpha = sp.Rational(31, 32)
    n = board.cells
    hats = {}
    for i, j in FLAGS:
        signs = sp.ones(n, 1)
        signs[edge[0]], signs[edge[1]] = (-1) ** i, (-1) ** j
        hats[i, j] = (sp.eye(n) - alpha * sp.diag(*signs) * p).inv() * (
            (1 - alpha) * signs
        )
    rows = [
        {
            (i, j): F(
                sum((-1) ** (i * a + j * b) * h[c] for (a, b), h in hats.items()) / 4
            )
            for i, j in FLAGS
        }
        for c in range(n)
    ]
    mirror = next(
        c for c in range(n) if all(rows[c][i, j] == rows[0][j, i] for i, j in FLAGS)
    )
    assert min(rows[0].values()) > 0 and sum(rows[0].values()) == 1
    return rows[0], rows[mirror], mirror


class EnergyAverage:
    def __init__(self, h):
        self.h = h
        self.eye = np.eye(len(h), dtype=object)
        t = sp.Symbol("t")
        minimal = (
            sp.Poly(sp.Matrix(h.tolist()).charpoly(t).as_expr(), t).sqf_part().monic()
        )
        self.p = [int(z) for z in minimal.all_coeffs()]
        inverse = sp.invert(minimal.diff(), minimal)
        rational = [F(c) for c in inverse.all_coeffs()]
        self.den = lcm(*(x.denominator for x in rational))
        self.b = np.zeros_like(h)
        for c in rational:
            self.b = h @ self.b + int(c * self.den) * self.eye
        intervals = minimal.intervals(eps=sp.Rational(1, 2**80))
        assert len(intervals) == minimal.degree() and all(
            mult == 1 for _, mult in intervals
        )
        self.gap = min(
            F(y[0][0] - x[0][1]) / 2
            for x, y in zip(intervals, intervals[1:], strict=False)
        )
        assert self.gap > 0
        # All intervals below the alias frequency for the native step.
        assert F(intervals[-1][0][1] - intervals[0][0][0]) * STEP / 2 < 3
        self.rank = len(intervals)
        pi_lo = F("3.14159265358979323846264338327950288")
        pi_hi = pi_lo + F(1, 10**35)
        waits = []
        errors = []
        for ((a, b), _), ((c, d), _) in combinations(intervals, 2):
            low, high = F(c - b) / 2, F(d - a) / 2
            duration = F(round((pi_lo + pi_hi) / (low + high) * 2**43), 2**43)
            waits.append(duration)
            errors.append(
                max(abs(low * duration - pi_hi), abs(high * duration - pi_lo)) / 2
            )
        self.waits = waits
        self.timing_error = max(errors) * (self.rank - 1)
        self.energies, self.vectors = np.linalg.eigh(np.array(h, float) / 2)
        self.clusters = np.split(
            np.arange(len(h)),
            np.flatnonzero(np.diff(self.energies) > float(self.gap) / 2) + 1,
        )
        assert len(self.clusters) == self.rank

    def pinch(self, z):
        # The divided difference of p(x) has value p'(x) on equal
        # eigenvalues and vanishes on distinct roots. Its diagonal inverse
        # from Bezout normalizes every surviving block, including degeneracies.
        source, den = z
        accum = np.zeros_like(source)
        poly = self.eye * self.p[0]
        for c in self.p[1:]:
            accum = self.h @ accum + source @ poly
            poly = self.h @ poly + c * self.eye
        assert not np.any(poly)
        num = self.b @ accum
        outden = den * self.den
        np.testing.assert_array_equal(self.h @ num, num @ self.h)
        np.testing.assert_array_equal(num.T, num)
        assert np.trace(num) * den == np.trace(source) * outden
        common = gcd(outden, *num.flat)
        return num // common, outden // common

    def finite_wait(self, rho):
        v = self.vectors
        z = v.T @ rho @ v
        # Random waits from independent fair bits. Propagate complete blocks
        # on degenerate energies together. The separately proved trace bound
        # encloses the exact finite instrument uniformly in subsequent delay.
        energies = self.energies.copy()
        for c in self.clusters:
            energies[c] = energies[c].mean()
        de = energies[:, None] - energies[None, :]
        for wait in self.waits:
            z *= (1 + np.exp(-1j * float(wait) * de)) / 2
        return v @ z @ v.T


def native_energies(deg, edges, second):
    operators = {}
    for flags in FLAGS:
        h = np.diag(16 * (6 - deg) - 80).astype(object)
        for a, b in edges:
            if any(z and (a, b) in s for z, s in zip(flags, second, strict=True)):
                continue
            h[a, a] += 17
            h[b, b] += 17
            h[a, b] -= 1
            h[b, a] -= 1
        operators[flags] = EnergyAverage(h)
    return operators


def scores(z, weights):
    a, b = z
    return F(sum(x * int(w) for x, w in zip(a.diagonal(), weights, strict=True)), b)


def traces(setup, flip, board, rng):
    _, cells, ports, ab, _, degree, edges, second = setup
    n = len(ports)
    exact, rho = (np.ones((n, n), dtype=object), n), np.ones((n, n), complex) / n
    ops = native_energies(degree, edges, second)
    restored = ops[0, 0]

    def available(j):
        return int(j not in ab and tuple(sorted((j, ab[0]))) in edges) - int(
            j not in ab and tuple(sorted((j, ab[1]))) in edges
        )

    score = [available(j) for j in range(n)]
    inside = [int(j // 12 == 0) for j in range(n)]
    forward, reverse, dispatch = parity(board, cells)
    if "joint_parity_and_stop" in setup[0]:
        for i, j in FLAGS:
            assert forward[i, j] == sum(map(F, setup[0]["joint_parity_and_stop"][i][j]))
    # pi < 22/7. Both stages contract trace norm. Limit total error
    # to one part in 10^9 before a source or a second probe is read.
    bound = max(F(22, 7) * (o.rank - 1) / (STEP * o.gap) for o in ops.values())
    bits = 0
    while 2 * bound / 2**bits > F(1, 10**9):
        bits += 1
    panel = motion.second_instruments(ports, ab, degree, edges)
    columns = np.eye(n, dtype=object)[:, ab]
    blocks = [columns]
    for _ in range(1, n):
        blocks.append(restored.h @ blocks[-1])
    kernel = (
        sp.Matrix(np.column_stack(blocks).T.tolist()).to_DM().nullspace().to_Matrix()
    )
    pd = kernel.T * (kernel * kernel.T).inv() * kernel
    rank = n - kernel.rows
    scale = lcm(*(F(c).denominator for c in pd))
    dark = np.array((scale * pd).tolist(), dtype=object).astype(object)
    dark = np.vectorize(int, otypes=[object])(dark)
    bright = np.eye(n, dtype=object) * scale - dark
    assert np.array_equal(dark @ restored.h, restored.h @ dark)
    observ = np.diag(score).astype(object)
    pieces = [m @ observ @ m for m in (dark, bright)]
    # Full position monitoring followed by energy averaging is a stochastic
    # matrix. A lower overlap bounds total-variation contraction on every pair.
    transition = []
    for j in range(n):
        unit = np.zeros((n, n), dtype=object)
        unit[j, j] = 1
        val, den = restored.pinch((unit, 1))
        transition.append([F(int(q) * 2**50 // den, 2**50) for q in val.diagonal()])
    contraction = 1 - min(
        sum(min(a, b) for a, b in zip(v, w, strict=True))
        for v, w in combinations(transition, 2)
    )
    quiet = restored.pinch(exact)
    outputs = {}
    for recorded in (False, True):
        per_flag = {}
        probs = {}
        positions = {}
        for flags, o in ops.items():
            first = o.pinch(exact)
            wave = o.finite_wait(rho)
            if recorded:
                first = (np.diag(first[0].diagonal()), first[1])
                wave = np.diag(wave.diagonal())
            final = restored.pinch(first)
            finite = restored.finite_wait(wave)
            z = np.array([float(F(v, final[1])) for v in final[0].diagonal()])
            observed = finite.diagonal().real
            assert np.abs(z - observed).sum() <= float(2 * bound / 2**bits) + 1e-12
            assert observed.min() >= 0 and abs(observed.sum() - 1) < 1e-12
            probs[flags] = observed
            positions[flags] = [F(v, final[1]) for v in final[0].diagonal()]
            per_flag[flags] = [
                scores(final, score),
                scores(final, inside),
                *(F(sum((op * final[0].T).flat), final[1] * scale**2) for op in pieces),
            ]
        arms = {}
        for name, law in (("forward", forward), ("reversed", reverse)):
            values = [
                sum(w * (per_flag[f][j] - per_flag[0, 0][j]) for f, w in law.items())
                / 144
                for j in range(4)
            ]
            accepted = sum(float(w) * probs[f] for f, w in law.items())
            reads = {}
            for tag, p in (("fresh", 25 / 36), ("shared", 101 / 144)):
                pos = p * accepted + (1 - p) * probs[0, 0]
                reads[tag] = motion.observe(
                    pos, panel, np.array(inside), score, ab, rng
                )
                qpass = F(25, 36) if tag == "fresh" else F(101, 144)
                reads[tag]["stationary_position_floor_2pow100"] = [
                    int(
                        2**100
                        * (
                            (1 - qpass) * positions[0, 0][k]
                            + qpass * sum(w * positions[f][k] for f, w in law.items())
                        )
                        // 1
                    )
                    for k in range(n)
                ]
            shared, fresh = reads["shared"], reads["fresh"]
            arms[name] = {
                "exact_available": str(values[0]),
                "exact_near": str(values[1]),
                "available_by_restored_modes": {
                    k: str(v)
                    for k, v in zip(
                        ("gateway_dark", "gateway_visible", "cross"),
                        (values[2], values[3], values[0] - values[2] - values[3]),
                        strict=True,
                    )
                },
                "available": local.difference(shared["available"], fresh["available"]),
                "near": local.difference(shared["near"], fresh["near"]),
                "reads": reads,
                "finite_flight": [
                    [
                        local.difference(x, y)
                        for x, y in zip(a["scores"], b["scores"], strict=True)
                    ]
                    for a, b in zip(
                        shared["held_flights"], fresh["held_flights"], strict=True
                    )
                ],
            }
        odd = (
            F(arms["forward"]["exact_available"])
            - F(arms["reversed"]["exact_available"])
        ) / 2
        outputs["position_read" if recorded else "unread"] = {
            "arms": arms,
            "exact_reversal_odd": str(odd),
            "release_odd": [str(c * odd) for c in (F(17, 2), F(65, 4), F(-31, 2))],
        }
    return {
        "sites": list(cells),
        "ports": ports,
        "gateway_indices": list(ab),
        "reflected_dispatch": dispatch,
        "cell_isomorphism": flip,
        "gateway_visible_dimension": rank,
        "dark_dimension": n - rank,
        "full_monitoring_contraction_upper": str(contraction),
        "spectrum": [
            {
                "flags": list(f),
                "squarefree_characteristic": o.p,
                "physical_gap_lower": str(o.gap),
                "fair_half_period_waits": len(o.waits),
                "maximum_total_duration": str(sum(o.waits)),
                "finite_trace_error_upper": str(o.timing_error),
            }
            for f, o in ops.items()
        ],
        "step": str(STEP),
        "bits_per_wait": bits,
        "uniform_clock_alternative_error": str(2 * bound / 2**bits),
        "combined_trace_error_upper": str(
            max(o.timing_error for o in ops.values()) + restored.timing_error
        ),
        "quiet_near": str(scores(quiet, inside)),
        "outcomes": outputs,
    }


def run():
    rng = np.random.default_rng(SEED)
    board = local.survey.bulk.federation(0)
    twin, mapping = symmetric_instrument(board)
    return {
        "seed": SEED,
        "original": traces(exp.flight.instrument(), None, board, rng),
        "mirror_matched": traces(twin, mapping, board, rng),
        "inputs": exp.save_inputs(
            "STATIONARY_POLARIZATION.md",
            "material_polarization.json",
            "routed_contacts.json",
        ),
    }


if __name__ == "__main__":
    result = run()
    exp.emit("stationary_polarization", result)
    for name in ("original", "mirror_matched"):
        row = result[name]
        print(name, row["sites"], row["bits_per_wait"], row["dark_dimension"])
        for k, r in row["outcomes"].items():
            print(
                k,
                [float(F(a["exact_available"])) for a in r["arms"].values()],
                float(F(r["exact_reversal_odd"])),
            )
