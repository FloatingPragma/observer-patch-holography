"""Matched optical jets and registered stress with different translation force."""

import argparse
import json
from pathlib import Path
from fractions import Fraction

import networkx as nx
import numpy as np
import sympy as sp
from mpmath import iv
from scipy.sparse.linalg import expm_multiply

from ...action_response.four_ticks import raw_effect
from .. import bridge
from ..traction.certificate import compress, four_tick_columns
from ..translation import measure as translation
from ..translation.counts import big_binomial, selection_norm

HERE = Path(__file__).resolve().parent


def controls(board, route, patch, q, odd):
    h, ops = translation.op_table(board, route, patch, q, odd)
    s = translation.map_matrix(board, patch)
    a = route[0]
    b = int(s[:, 4 * route[1]].nonzero()[0][0]) // 4
    sites = bridge.beam.indices([a, b])

    def invisible(mat):
        return mat[4 * a : 4 * a + 4, 4 * b : 4 * b + 4].count_nonzero() == 0

    for tag in ("shift", "domain_cut", "image_cut", "boundary"):
        assert invisible(ops[tag][0])
    flux = translation.flux_on_edges(board, patch)
    assert all(invisible(v) for v in flux.values())
    # Q_u for every cell and every lifting of its fixed registered translation.
    for i in range(len(patch["cells"])):
        part = translation.map_matrix(
            board, {k: [patch[k][i]] for k in ("cells", "images", "lifts")}
        )
        assert invisible(part - part.T)
    # Independently measured native cell and contact energies, bridge work, and
    # energy currents have support on a single cell or adjacent cell dyads.
    assert (
        nx.shortest_path_length(
            nx.Graph((p // 12, t // 12) for p, t in board.external), a // 12, b // 12
        )
        == 2
    )
    graph = nx.Graph(board.seams)
    distance = nx.shortest_path_length(graph, a, b)
    lp = nx.single_source_shortest_path_length(graph, a)
    lt = nx.single_source_shortest_path_length(graph, b)
    assert min(lp[y] + lt[y] for y in graph) == distance
    # All translated chart partners preserve gateway type, even if a different
    # local port lifting was enrolled. Every translation flux uses only gateways.
    selector = np.repeat(board.degree == 6, 4).astype(int)
    for _, m in [("shift", ops["shift"][0]), *flux.items()]:
        co = m.tocoo()
        assert np.all(selector[co.row] == selector[co.col])
    graph.remove_edges_from((route[:2], route[2:]))
    assert nx.has_path(graph, a, b)
    assert distance == 9 and board.degree[a] == 5 and board.degree[b] == 6
    f4 = ops["marked_force"][0][sites, :][:, sites].toarray()
    pair = np.zeros((8, 8), dtype=np.int64)
    pair[:4, 4:] = f4[:4, 4:]
    return h, sites, (pair + pair.T), distance, ops["marked_force"][0]


def optical_jets(w, h, receivers):
    vel = (-0.5j) * (h @ w)
    acc = (-0.5j) * (h @ vel)
    a, b, c = [
        rows[receivers].reshape(10, 4, 10, 4).transpose(0, 2, 1, 3)
        for rows in (w, vel, acc)
    ]

    def adj(x):
        return x.conj().transpose(0, 1, 3, 2)

    return np.stack(
        (adj(a) @ a, adj(a) @ b + adj(b) @ a, adj(a) @ c + 2 * adj(b) @ b + adj(c) @ a)
    )


def run():
    iv.dps = 30
    atlas = translation.Atlas(1)
    board = atlas.board
    route = bridge.native.iso.pairing.formation.region(board)
    patch = translation.chosen(atlas)
    rings = bridge.beam.screens(board, route)
    launches = bridge.beam.indices(sum((list(r) for r in rings), []))
    source = bridge.action.front.inputs()
    coeff = bridge.native.prior.weights()[7]
    states = bridge.action.preparations()
    angles = {
        ("ordinary" if row["source"] == "ordinary" else "pairing"): row[
            "angle_per_unit"
        ]
        for row in json.loads(
            (bridge.HERE.parent / "signed_transfers.json").read_text()
        )["tests"]
    }
    rng = np.random.default_rng(17021926)
    n = 2**76
    output = {}
    for q in (-1, 1):
        contracted = []
        direct = []
        prepared = []
        sector_norm = []
        protocols = {"together": [], "opposite": []}
        optical = {"together": [], "opposite": []}
        for odd in (False, True):
            h, idx, cross, distance, force4 = controls(board, route, patch, q, odd)
            # The two phase lotteries differ by twice the off-diagonal block
            # of F; f4=4F, hence the registered contrast effect is cross/2.
            w = four_tick_columns(h, idx, launches, coeff)
            contracted.append(compress(w, cross, 2))
            c = sp.Matrix(cross[:4, 4:] + cross[4:, :4]) / 4
            prepared.append(c)
            ht = -2j * h
            start = np.zeros((h.shape[0], len(launches)))
            start[launches] = np.eye(len(launches))
            all_w = expm_multiply(ht, start, traceA=ht.diagonal().sum())
            z = all_w[idx]
            direct.append(
                np.mean([v.conj().T @ cross @ v / 2 for v in np.hsplit(z, 10)], axis=0)
            )
            for label, settings in (
                ("together", ((1, 1), (-1, -1))),
                ("opposite", ((-1, 1), (1, -1))),
            ):
                effect = np.zeros((4, 4), dtype=complex)
                jets_optical = 0
                for x, y in settings:
                    twisted = all_w.copy()
                    twisted[idx[:4]] *= x
                    twisted[idx[4:]] *= y
                    jets_optical += optical_jets(twisted, h, launches) / 2
                    acted = force4 @ twisted
                    effect += (
                        sum(
                            v.conj().T @ a
                            for v, a in zip(
                                np.hsplit(twisted, 10), np.hsplit(acted, 10)
                            )
                        )
                        / 80
                    )
                protocols[label].append(effect)
                optical[label].append(jets_optical)
            np.testing.assert_allclose(
                protocols["together"][-1] - protocols["opposite"][-1],
                direct[-1],
                atol=1e-15,
            )
            sector_norm.append(selection_norm(force4, 4))
        difference = np.abs(
            np.asarray(optical["together"]) - np.asarray(optical["opposite"])
        )
        assert difference.max() < 1e-10
        base = np.asarray(optical["together"])[:, 0]
        plane = bridge.beam.probabilities(base, states["equilibrium"])
        sigma = [
            np.linalg.svd(bridge.beam.image(sub)[0], compute_uv=False)[-1]
            for sub in (plane[5:, :5], plane[:5, 5:])
        ]
        assert min(sigma) > 1e-5
        # Independently sample the force meter on the two randomly phased arms.
        norm = max(sector_norm)
        expected = bridge.evaluate(direct, states)
        counted = {}
        for label, panel in protocols.items():
            individual = bridge.evaluate(panel, states)
            records = {}
            for name, x in individual.items():
                plus = big_binomial(rng, n, 0.5 + x / (2 * norm))
                records[name] = {
                    "plus": plus,
                    "mean": x,
                    "estimate": float(Fraction(2 * plus - n, n)) * norm,
                }
            counted[label] = records
        means = {
            key: {
                "estimate": counted["together"][key]["estimate"]
                - counted["opposite"][key]["estimate"]
            }
            for key in states
        }
        eps = 2 * norm * np.sqrt(2 * np.log(2 * 2 * 2 * len(states) / 0.01) / float(n))
        assert abs(means["equilibrium"]["estimate"] - expected["equilibrium"]) <= eps
        samples = {}
        for name, ang in angles.items():
            m = {key: v["estimate"] for key, v in means.items()}
            slope = m[name] - m["equilibrium"] + ang * (m["plus"] - m["minus"]) / 2
            width = eps * (2 + abs(ang))
            nominal = bridge.evaluate(direct, bridge.tangents())[name]
            assert abs(slope - nominal) < width and abs(slope) > width
            samples[name] = dict(estimate=slope, halfwidth=width)
        jets = bridge.action.source_jets(contracted, source)
        jets = {k: v + iv.mpf(["-1e-12", "1e-12"]) for k, v in jets.items()}
        val = raw_effect(
            contracted, source[0], sp.Rational("6e-7"), sp.Rational("1e-15")
        )
        pair = raw_effect(prepared, source[0], sp.Rational("6e-7"), 0)
        for key, bound in [("reference", val), *jets.items()]:
            assert (
                0 not in bound
                and bridge.evaluate(direct, bridge.tangents())[key] in bound
            )
        output[str(q)] = dict(
            after_original_flight=dict(
                reference=str(val),
                ordinary=str(jets["ordinary"]),
                pairing=str(jets["pairing"]),
            ),
            freshly_split_from_registered_heads=str(pair),
            synthetic_counts=counted,
            count_slopes=samples,
            binary_scale=norm,
            reference_count_halfwidth=eps,
            maximum_direct_optical_table_and_jet_difference=float(difference.max()),
            reference_minimum_image_singular_values=list(map(float, sigma)),
        )
    return dict(
        ports=[int(idx[0] // 4), int(idx[4] // 4)],
        native_seam_distance=distance,
        highest_common_position_diagonal_jet=distance - 1,
        tagged_dwell=4,
        equal_registered_stress_contrasts=0,
        simulated_trials_per_setting=n,
        sampling_seed=17021926,
        coverage=0.99,
        charges=output,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (HERE / "results.json").write_text(result)
    else:
        print(result, end="")
