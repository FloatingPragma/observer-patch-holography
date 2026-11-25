"""Route-resolved neutral transport with source refusals and two local signed receivers."""

import json
from pathlib import Path
from functools import cache
import numpy as np
import networkx as nx
from scipy.linalg import eigh
from scipy.sparse.linalg import expm_multiply
from scipy.stats import beta

from ....apparatus import native_graph, fixture
from ...instruments import material
from .detectors import verify
from .geometry import configure, layout
from .flight import TIMES, frame, chart
from ..measure import gate, WORDS

HERE = Path(__file__).resolve().parent
SHORT = 2**38
LARGE = 2**46


def source_records():
    f = fixture()
    marks = layout()[1]
    heads = [list(map(tuple, f["marked"][3:])).index(e) for e in marks]
    src = material("hub")["character"]
    h = native_graph()
    y = np.zeros((h.shape[0], 3), complex)
    y[:8] = src
    gen = -32j * h
    received = expm_multiply(gen, y, traceA=gen.diagonal().sum())[64:128].reshape(
        8, 8, 3
    )
    # The two side heads store independent leaf currents during release.
    # Distinct input head sectors remain orthogonal on every charge-number read.
    prob = np.sum(abs(received) ** 2, axis=1)
    masks = np.arange(8)
    labels = sum(((masks >> j) & 1) << i for i, j in enumerate(heads))
    weights = np.array([prob[labels == i].sum(axis=0) for i in range(4)])
    np.testing.assert_allclose(np.sum(abs(src) ** 2, axis=0), 1 / 3, atol=1e-12)
    assert (h != native_graph(1)).nnz == 0
    return weights


def interval(k, n):
    tail = 0.01 / (2 * 4096)
    return [
        float(beta.ppf(tail, k, n - k + 1)) if k else 0.0,
        float(beta.isf(tail, k + 1, n - k)) if k < n else 1.0,
    ]


def bank(p, source, n, rng):
    full = np.r_[p * source, 1 / 3 - source, 2 / 3]
    assert min(full) > -1e-13 and abs(sum(full) - 1) < 2e-12
    return rng.multinomial(n, np.maximum(full, 0) / sum(np.maximum(full, 0))).tolist()


def joined(p):
    return float(p[2] + p[6])


def reread(p):
    a = p.reshape(3, 3)
    return np.outer(a.sum(axis=1), a.sum(axis=0)).ravel()


def collect(rows, near, headweights, source, sign, k, bench, rng):
    results = []
    for i, t in enumerate(TIMES):
        ok = float(sum(headweights[j] * near[j][i] for j in range(len(headweights))))
        close = bank(np.array([ok, 1 - ok]), source, SHORT, rng)
        for sensor in ("opposite", "second_owner"):
            p = sum(headweights[j] * rows[j][sensor][i] for j in range(len(rows)))
            results.append(
                [
                    sign,
                    k,
                    bench,
                    t,
                    sensor,
                    close,
                    bank(p, source, SHORT, rng),
                    ok,
                    joined(p),
                    joined(reread(p)),
                ]
            )
    return results


def held(p, source, n, rng):
    counts = bank(p, source, n, rng)
    k = counts[2] + counts[6]
    admitted = sum(counts[:9])
    return dict(
        counts=counts,
        coincidences=k,
        admitted=admitted,
        conditional_99=interval(k, admitted),
        original_99=interval(k, n),
    )


@cache
def one_defect(where, sign, mask, written):
    start = (where, -1, mask) if sign == -1 else (-1, where, mask)
    states, h = chart(start, written)
    ev, U = eigh(h, driver="evd")
    p = abs((np.exp(-1j * np.outer(TIMES, ev)) * U[0]) @ U.T) ** 2
    d = np.array(states)[:, 0 if sign == -1 else 1]
    return {v: p[:, d == v].sum(axis=1) for v in set(d)}


def independent_flights(sign, weight, written, sensor, idx):
    a, b = (207, 206) if sensor == "opposite" else (106, 96)
    # Independently prepared single-sign launches with the same inferred
    # marginal side-head distribution. Compare both orders at the far pair.
    singles = []
    for c, site in ((sign, 160), (-sign, 167)):
        singles.append(
            [
                sum(
                    p * one_defect(site, c, m, written)[v][idx]
                    for m, p in enumerate(weight)
                )
                for v in (a, b)
            ]
        )
    return float(singles[0][0] * singles[1][1] + singles[1][0] * singles[0][1])


def run():
    weights = source_records()
    ring, marks, leaves, owners, candidates = layout()
    rng = np.random.default_rng(35124)
    new = np.random.default_rng(35125)
    setup = []
    spectra = []
    primary = []
    controls = []
    separated = []
    graph = nx.Graph(zip(ring, ring[1:] + ring[:1]))
    graph.add_edges_from(marks)
    cuts = [e for e in graph.edges if (e[0] // 12 == 13) != (e[1] // 12 == 13)]
    assert len(cuts) == 2
    graph.remove_edges_from(cuts)
    assert all(not nx.has_path(graph, 160, r) for r in (207, 206, 106, 96))
    for sign in (-1, 1):
        c = WORDS.index((1 + sign, 1, 1 - sign))
        n = WORDS.index((1 + sign, 1 - sign, 1))
        assert abs(gate(np.pi, 1 - sign)[n, c]) ** 2 > 1 - 1e-12
        for written in (True, False):
            bench = "intact" if written else "retired"
            setup.append(configure(sign, not written))
            fields = []
            near = []
            for m in range(4 if written else 1):
                d, z, f = frame(sign, m, written=written)
                spectra.append(dict(sign=sign, heads=m, bench=bench, **d))
                fields.append(f)
                near.append(z)
            for k in range(3):
                prob = float(weights[:, k].sum())
                w = weights[:, k] / prob if written else np.ones(1)
                if written and sign > 0:
                    w = w[::-1]
                primary.extend(collect(fields, near, w, prob, sign, k, bench, rng))
                t = 2048 if written else 256
                i = TIMES.index(t)
                for sensor in ("opposite", "second_owner"):
                    p = sum(w[j] * r[sensor][i] for j, r in enumerate(fields))
                    blocked = np.eye(9)[4]
                    for shots in (SHORT, LARGE):
                        data = {
                            name: held(dist, prob, shots, new)
                            for name, dist in [
                                ("free", p),
                                ("independent_marginals", reread(p)),
                                ("two_gates_closed", blocked),
                            ]
                        }
                        controls.append(
                            dict(
                                sign=sign,
                                k=k,
                                bench=bench,
                                t=t,
                                sensor=sensor,
                                starts=shots,
                                pair_given_launch=joined(p),
                                repaired_independent_marginals_given_launch=joined(
                                    reread(p)
                                ),
                                independent_single_load_product=independent_flights(
                                    sign, w, written, sensor, i
                                ),
                                data=data,
                            )
                        )
            if not written:
                dd, close, ra = frame(sign, 0, (160, 152), False)
                separated.append(
                    dict(
                        sign=sign,
                        initial=[160, 152],
                        high_energy_weight=dd["high_energy_weight"],
                        adjacent=close.tolist(),
                        opposite=(ra["opposite"][:, 2] + ra["opposite"][:, 6]).tolist(),
                    )
                )
    return dict(
        local_receiver_meters=verify(),
        source="original three native character sources; free signed flight 32, full ladder at 160/167",
        seeds=[35124, 35125],
        ring=ring,
        owners=owners,
        unknown_side_heads=marks,
        leaves=leaves,
        geometry_candidates=candidates,
        blocked_seams=cuts,
        original_source_to_neutral_probabilities=weights.sum(axis=0).tolist(),
        unconditioned_input_head_weights=weights.tolist(),
        times=TIMES,
        grid_starts=SHORT,
        high_statistics_starts=LARGE,
        count_order=[
            "00",
            "01",
            "02",
            "10",
            "11",
            "12",
            "20",
            "21",
            "22",
            "no far-hub charge",
            "original source refused",
        ],
        grid_columns=[
            "charge",
            "winding",
            "apparatus",
            "time",
            "detector",
            "adjacency_bank_[yes,no,miss,refusal]",
            "signed_counts",
            "expected_adjacency",
            "expected_joint",
            "expected_independent_marginals",
        ],
        setup=setup,
        spectra=spectra,
        grid=primary,
        held=controls,
        separated_control=separated,
    )


if __name__ == "__main__":
    data = run()
    extra = {
        k: data.pop(k)
        for k in ("setup", "spectra", "grid", "held", "separated_control")
    }
    body = json.dumps(data, indent=2)[:-2]
    for k, rows in extra.items():
        body += (
            ",\n  "
            + json.dumps(k)
            + ": [\n    "
            + ",\n    ".join(json.dumps(r) for r in rows)
            + "\n  ]"
        )
    (HERE / "results.json").write_text(body + "\n}\n")
