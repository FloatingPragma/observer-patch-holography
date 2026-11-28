"""Test every deterministic local summary of the two receiving write histories."""

from itertools import combinations, product
from math import exp, log, sqrt

from mpmath import iv
import networkx as nx
import numpy as np

from . import archived_overlap as archive

coherent = archive.coherent
local = archive.local
SEED = 28122923
N = coherent.N


def partitions(visible):
    # Canonical names by order of first occurrence give one representative
    # of every partition, regardless of which labels observers write down.
    inputs = [tuple(f[j] for j in visible) for f in archive.FLAGS]
    labels = list(dict.fromkeys(inputs))
    lookup = [labels.index(v) for v in inputs]
    for code in product(range(len(labels)), repeat=len(labels)):
        if code[0] or any(code[j] > max(code[:j]) + 1 for j in range(1, len(code))):
            continue
        yield tuple(code[j] for j in lookup)


def tomography(layouts, edges, rng):
    graph = nx.Graph(edges)
    raw = {}
    matrices = []
    for f in archive.FLAGS:
        x = layouts[f]
        pos = coherent.draws(x.diagonal(), rng)
        m = np.diag([k / N for k in pos]).astype(complex)
        pairs = []
        for i, j in combinations(range(len(x)), 2):
            # Route the two named ports together, add a controlled phase,
            # apply a balanced partial swap and undo the permutation.
            route = nx.shortest_path(graph, i, j)
            mean = (x[i, i] + x[j, j]) / 2
            trials = [
                coherent.draws([mean + d, mean - d, 1 - 2 * mean], rng)
                for d in (x[i, j].real, x[i, j].imag)
            ]
            re, im = [(k[0] - k[1]) / (2 * N) for k in trials]
            m[i, j] = complex(re, im)
            m[j, i] = m[i, j].conjugate()
            pairs.append({"path": route, "quadratures": trials})
        raw[str(f)] = {"positions": pos, "pairs": pairs}
        matrices.append(m)
    return raw, matrices


def spectral_entropy(values, error):
    lo = hi = iv.mpf(0)

    def h(x):
        v = iv.mpf(x)
        return -v * iv.log(v) if x > 0 else iv.mpf(0)

    for v in values:
        a = max(0, float(v) - error)
        b = min(1, float(v) + error)
        assert a <= b
        ha, hb = h(a), h(b)
        lo += min(float(ha.a), float(hb.a)) - 1e-15
        hi += (
            1 / exp(1) if a <= 1 / exp(1) <= b else max(float(ha.b), float(hb.b))
        ) + 1e-15
    return iv.mpf([float(lo.a) - 1e-14, float(hi.b) + 1e-14])


def entropies(sites, matrices, weights, entry_error, weight_error):
    d = len(sites)
    table = {}
    for size in range(1, 5):
        for key in combinations(range(4), size):
            blocks = []
            for w in weights:
                occupied = sum(w[i] * matrices[i][np.ix_(sites, sites)] for i in key)
                q = sum(w[i] for i in key)
                vacuum = q - occupied.trace().real
                eig, v = np.linalg.eigh(occupied)
                # Dimension <=24 and normalized matrices. The fixed allowance
                # covers the reconstruction and Gram products (gamma_24)
                # and the polar orthogonalization of this measured eigenbasis.
                assert max(abs(occupied - (v * eig) @ v.conj().T).sum(axis=1)) < 1e-12
                assert max(abs(v.conj().T @ v - np.eye(d)).sum(axis=1)) < 1e-12
                bound = q * (2 * d * entry_error) + size * weight_error + 2e-11
                blocks.append(spectral_entropy([*eig, vacuum], bound))
            table[key] = blocks
    return table


def survey_ball(sites, gates, matrices, weights, err, werr):
    table = entropies(sites, matrices, weights, err, werr)
    visible = [i for i, g in enumerate(gates) if g in sites]
    outcomes = {}
    for p in partitions(visible):
        groups = [tuple(i for i, v in enumerate(p) if v == u) for u in sorted(set(p))]
        minus, mid, plus = [sum(table[g][j] for g in groups) for j in range(3)]
        # Entropy is concave along the affine source preparation. Both endpoint
        # differences bound its slope at the fresh reference, including zero
        # eigenvalues because these three mixtures have identical support.
        ds = iv.mpf([float((plus - mid).a), float((mid - minus).b)])
        outcomes["".join(map(str, p))] = (mid, ds)
    return outcomes


def run():
    iv.dps = 45
    rng = np.random.default_rng(SEED)
    board = local.survey.bulk.federation(0)
    setup, _ = archive.steady.symmetric_instrument(board)
    layouts, _ = coherent.stationary_layouts(setup)
    raw, matrices = tomography(layouts, setup[6], rng)
    saved = local.inputs("restored_balls")["surveys"]
    wait = float(
        archive.F(
            local.inputs("stationary_polarization")["mirror_matched"][
                "combined_trace_error_upper"
            ]
        )
    )
    err = sqrt(log(2 * 100000 / 0.01) / (2 * N)) + 5e-13
    fwd, rev, _ = archive.steady.parity(board, setup[1])
    records = {}
    output = []
    for direction, law in (("forward", fwd), ("reversed", rev)):
        weights = []
        count = []
        for dose in (99 / 144, 100 / 144, 101 / 144):
            wf = [
                (1 - dose) * int(f == (0, 0)) + dose * float(law[f])
                for f in archive.FLAGS
            ]
            hits = coherent.draws(wf, rng)
            weights.append([v / N for v in hits])
            count.append(hits)
        cache = {}
        for row in saved:
            if row["direction"] != direction:
                continue
            inner, outer = [tuple(z) for z in row["ball_sites"]]
            for ball in (inner, outer):
                if ball not in cache:
                    cache[ball] = survey_ball(
                        list(ball), setup[3], matrices, weights, err + wait, err
                    )
            area = iv.mpf(row["record_source"]["area_tangent"])
            coin = iv.mpf(row["record_source"]["public_probability_tangent"])
            gains = {}
            for left, (s0, ds0) in cache[inner].items():
                for right, (s1, ds1) in cache[outer].items():
                    gain = -((ds0 + ds1) / 2 + coin * (s1 - s0)) / area
                    gains[left + "_" + right] = archive.bounds(gain)
            largest = max(gains, key=lambda p: gains[p][1])
            output.append(
                {
                    "center": row["center"],
                    "direction": direction,
                    "radii": row["hop_radii"],
                    "gains": gains,
                    "largest_upper": gains[largest][1],
                    "largest_keys": largest,
                }
            )
        records[direction] = count
    return {
        "seed": SEED,
        "starts_per_setting": N,
        "history_conditional_tomography": raw,
        "archive_counts_opposite_fresh_shared": records,
        "surveys": output,
        "inputs": archive.steady.exp.save_inputs(
            "ARCHIVED_OVERLAP.md",
            "restored_balls.json",
            "stationary_polarization.json",
            "shared_contact.json",
            "archived_overlap.json",
        ),
    }


if __name__ == "__main__":
    z = run()
    archive.steady.exp.emit("archive_partitions", z)
    print([(r["center"], r["radii"], r["largest_upper"]) for r in z["surveys"]])
