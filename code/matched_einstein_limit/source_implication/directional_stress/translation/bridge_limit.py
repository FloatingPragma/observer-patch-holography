"""Support obstruction and independent smooth-wave and energy-shell bounds on the extended region."""

import argparse
import json

import networkx as nx
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from .. import bridge
from . import measure


def path_tables(board, route, patch):
    graph = nx.Graph(board.seams)
    graph.remove_edges_from((route[:2], route[2:]))
    tails = []
    heads = []
    rough = sparse.lil_matrix((board.n, board.n), dtype=np.int64)
    lengths = []
    for u, v, perm in zip(patch["cells"], patch["images"], patch["lifts"]):
        for i, j in enumerate(perm):
            path = nx.shortest_path(graph, 12 * u + i, 12 * v + j)
            m = len(path) - 1
            lengths.append(m)
            for a, b in zip(path, path[1:]):
                tails.append(a)
                heads.append(b)
                rough[a, a] += m - 1
                rough[b, b] += m - 1
                rough[a, b] -= m - 1
                rough[b, a] -= m - 1
    v = sparse.csr_matrix(
        (np.ones(len(tails), dtype=np.int64), (heads, tails)), shape=(board.n, board.n)
    )
    eye = sparse.eye(4, dtype=np.int64)
    return (
        sparse.kron(v - v.T, eye, format="csr") / (2j),
        sparse.kron(rough.tocsr(), eye, format="csr") / 2,
        lengths,
    )


def support_witness(board, patch, ops):
    q, den = ops["shift"]
    # Disjoint support from *every* active per-edge energy current follows
    # directly from the operator form: both cells involved are neighbors.
    a, b = q.nonzero()
    graph = nx.Graph((i // 12, j // 12) for i, j in board.external)
    dist = dict(nx.all_pairs_shortest_path_length(graph, cutoff=2))
    assert {dist[i // 48][j // 48] for i, j in zip(a, b)} == {2}
    s = measure.map_matrix(board, patch)
    p = 0
    dest = int(s[:, 4 * p].nonzero()[0][0]) // 4
    cut = [4 * p, 4 * dest]
    mat = q[cut, :][:, cut].toarray() / den
    reads = []
    for sign in (-1, 1):
        amp = np.array([1, sign * 1j]) / np.sqrt(2)
        reads.append(float(np.vdot(amp, mat @ amp).real))
    np.testing.assert_allclose(reads, [0.5, -0.5], atol=1e-14)
    return dict(
        ports=[p, dest],
        phase_order=["-i", "+i"],
        shift=reads,
        normalized_shift_contrast=1,
        energy_current_contrast_on_every_bridge=0,
    )


def irreducible_force(h, force4):
    """Cyclic-trace obstruction to absorbing this source at fixed boundary flux."""
    f = force4.tocoo()
    ix = np.unique(f.col)
    columns = np.zeros((h.shape[0], len(ix)), dtype=np.int64)
    columns[ix] = np.eye(len(ix), dtype=np.int64)
    row_bound = int(abs(h).sum(axis=1).max())
    trace = []
    for k in range(8):
        assert row_bound**k < 2**62
        trace.append(
            sum(
                int(v) * int(w)
                for v, w in zip(f.data, columns[f.row, np.searchsorted(ix, f.col)])
            )
        )
        if k < 7:
            columns = h @ columns
    assert trace == [0] * 7 + [64]
    assert row_bound < 64
    return dict(
        integer_traces=trace,
        first_nonzero_degree=7,
        native_trace_H7_force="1/8",
        dimension=h.shape[0],
        twice_generator_spectral_norm_upper=row_bound,
    )


def run():
    atlas = measure.Atlas(1)
    board = atlas.board
    region = measure.chosen(atlas)
    route = bridge.native.iso.pairing.formation.region(board)
    rings = bridge.beam.screens(board, route)
    src = bridge.beam.indices(list(rings[0]) + list(rings[1]))
    path, rough, lengths = path_tables(board, route, region)
    mask = np.repeat(
        np.isin(np.arange(board.cells), region["cells"]).astype(np.int64), 48
    )
    proj = sparse.diags(mask, dtype=np.int64, format="csr")
    diag = np.zeros(board.n, dtype=np.int64)
    for a, b in board.external:
        if mask[4 * a] != mask[4 * b]:
            diag[a] += 17 * (1 - 2 * mask[4 * a])
            diag[b] += 17 * (1 - 2 * mask[4 * b])
    delta4 = sparse.diags(np.repeat(diag, 4), dtype=np.int64, format="csr")
    states = bridge.action.preparations()
    row = {}
    irreducible = {}
    for charge in (-1, 1):
        meters = {}
        for odd in (False, True):
            h, ops = measure.op_table(board, route, region, charge, odd)
            irreducible[f"{charge}_{int(odd)}"] = irreducible_force(
                h, ops["marked_force"][0]
            )
            if not odd and charge == -1:
                witness = support_witness(board, region, ops)
            e4, energy_rate, flows = measure.regional_energy(
                board, route, region, charge, odd, h
            )
            # Form the mobile excitation H explicitly from *all* native seam
            # terms. This fixes its energy origin relative to the parked heads.
            cells = bridge.operators(board, route, charge, odd)
            h2 = bridge.masked(cells, range(board.cells)) + bridge.bridges(
                board, board.external
            )
            assert (e4 - h2 @ proj - proj @ h2 - delta4).count_nonzero() == 0
            qmat = ops["shift"][0] / ops["shift"][1]
            for j in flows.values():
                assert qmat.multiply(j).count_nonzero() == 0
            velocity = 0.5j * (h2 @ proj - proj @ h2)
            energy = 0.5 * h2
            storage_rate = 0.125j * (h2 @ delta4 - delta4 @ h2)
            operators = dict(
                shift=qmat,
                summed_path_orientation=path,
                path_error_bound=rough,
                energy=energy,
                energy_square=energy @ energy,
                load_inflow=velocity,
                load_inflow_square=velocity @ velocity,
                energy_inflow=0.125j * energy_rate,
                storage_rate=storage_rate,
            )
            # Entire operator identity before computing any covariance.
            assert (
                operators["energy_inflow"]
                - storage_rate
                - (energy @ velocity + velocity @ energy) / 2
            ).count_nonzero() == 0
            start = np.zeros((h.shape[0], len(src)))
            start[src] = np.eye(len(src))
            ht = -2j * h
            w = expm_multiply(ht, start, traceA=ht.diagonal().sum())
            for key, m in operators.items():
                current = m @ w
                meters.setdefault(key, []).append(
                    np.mean(
                        [
                            x.conj().T @ y
                            for x, y in zip(np.hsplit(w, 10), np.hsplit(current, 10))
                        ],
                        axis=0,
                    )
                )
        values = {k: bridge.evaluate(v, states) for k, v in meters.items()}
        received = {}
        for key in states:
            v = {k: x[key] for k, x in values.items()}
            error = v["shift"] - v["summed_path_orientation"]
            delta_e = np.sqrt(max(0, v["energy_square"] - v["energy"] ** 2))
            load_std = np.sqrt(max(0, v["load_inflow_square"] - v["load_inflow"] ** 2))
            err = (
                v["energy_inflow"] - v["storage_rate"] - v["energy"] * v["load_inflow"]
            )
            assert (
                abs(error) <= v["path_error_bound"] and abs(err) <= delta_e * load_std
            )
            received[key] = {
                **v,
                "actual_path_error": error,
                "energy_standard_deviation": delta_e,
                "energy_load_covariance": err,
                "energy_covariance_bound": delta_e * load_std,
            }
        row[str(charge)] = received
    return dict(
        irreducible_marked_force=irreducible,
        coherence_witness=witness,
        registered_path_lengths=dict(
            min=min(lengths), max=max(lengths), paths=len(lengths)
        ),
        charges=row,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (measure.HERE / "bridge_limit.json").write_text(data)
    else:
        print(data, end="")
