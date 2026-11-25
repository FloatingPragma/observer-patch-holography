"""Charge-independent leading-arrival clocks with a fixed hop order per local link."""

import argparse
from itertools import combinations
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components, shortest_path
from scipy.optimize import brentq

from .. import signed_transfers
from ..full_algebra import native_precision as native
from ..geometry_gain import simplex_area_volume
from ..null_optics import measure as optics


HERE = Path(__file__).resolve().parent


def hop(board, route, charge, odd):
    h = native.generator(board, route, charge, odd).tocoo()
    mask = (h.row // 4 != h.col // 4) & (h.data != 0)
    matrix = sparse.csr_matrix(
        (-h.data[mask], (h.row[mask], h.col[mask])), shape=h.shape
    )
    assert matrix.data.min() > 0
    return matrix


def path_blocks(board, route, order):
    """Exact integer sums on all shortest paths of this order, for both orientations."""
    hops = shortest_path(board.adjacency(), unweighted=True, directed=False)
    left, right = np.where(np.triu(hops == order, 1))
    ops = []
    labels = np.arange(4)
    at_i = 4 * left[:, None, None] + labels[None, :, None]
    at_j = 4 * right[:, None, None] + labels[None, None, :]
    indices = (at_i.ravel().repeat(4), np.broadcast_to(at_j, (len(left), 4, 4)).ravel())
    for odd in (False, True):
        blocks = []
        for q in (1, -1):
            h = hop(board, route, q, odd)
            degree = int(abs(h).sum(axis=1).max())
            assert 8 * degree ** (2 * order) < 2**63
            power = h**order
            blocks.append(np.asarray(power[indices]).reshape(-1, 4, 4))
        np.testing.assert_array_equal(blocks[0].transpose(0, 2, 1), blocks[1])
        m = blocks[0]
        adj = m.transpose(0, 2, 1)
        pair = m @ adj + adj @ m
        assert np.all(np.trace(pair, axis1=-1, axis2=-2) > 0)
        ops.append(pair)
    return left, right, np.stack(ops)


def graph(n, stations, weights, edges):
    a, b = edges
    matrix = sparse.csr_matrix(
        (np.r_[weights, weights], (np.r_[a, b], np.r_[b, a])), shape=(n, n)
    )
    d, prev = shortest_path(
        matrix, directed=False, indices=stations, return_predecessors=True
    )
    return d[:, stations], prev, matrix


def gram(distances):
    n = len(distances)
    j = np.eye(n) - np.ones((n, n)) / n
    return -j @ (distances**2) @ j / 2


class Ruler:
    def __init__(self, board, route, order):
        self.n = board.n
        self.route = route
        self.order = order
        a, b, self.ops = path_blocks(board, route, order)
        self.edges = a, b
        self.ids = {(int(i), int(j)): k for k, (i, j) in enumerate(zip(a, b))}
        # Multiple positions frequently have the identical coin effect.
        flattened = self.ops.transpose(1, 0, 2, 3).reshape(len(a), 32)
        unique, self.index = np.unique(flattened, axis=0, return_inverse=True)
        self.menu = unique.reshape(-1, 2, 4, 4)

    def blank_reference(self, board):
        adjacency = board.adjacency().astype(np.int64)
        assert np.all(adjacency.data == 1)
        remote = shortest_path(
            adjacency, directed=False, unweighted=True, indices=self.route
        ).min(axis=0)
        a, b = self.edges
        candidates = np.flatnonzero(remote[a] > self.order)
        counts = np.asarray(
            (adjacency**self.order)[a[candidates], b[candidates]]
        ).ravel()
        i = candidates[np.flatnonzero(counts == 1)[0]]
        np.testing.assert_array_equal(self.ops[:, i], np.tile(2 * np.eye(4), (2, 1, 1)))
        return [int(a[i]), int(b[i])]

    def coefficients(self, state):
        return np.einsum("kbij,bji->k", self.menu, state).real / 2

    def distances(self, state):
        weights = self.coefficients(state)[self.index] ** (-1 / (2 * self.order))
        return graph(self.n, (*self.route, 5), weights, self.edges)

    def observe(self, state, tangents):
        coeff = self.coefficients(state)
        weights = coeff ** (-1 / (2 * self.order))
        d, prev, matrix = self.distances(state)
        assert connected_components(matrix)[0] == 1
        paths = []
        derivatives = {k: np.zeros((4, 4)) for k in tangents}
        edge_derivative = {
            k: -weights * self.coefficients(direction) / (2 * self.order * coeff)
            for k, direction in tangents.items()
        }
        for i, j in combinations(range(4), 2):
            p, q = self.route[i], self.route[j]
            nodes = [q]
            while nodes[-1] != p:
                nodes.append(int(prev[i, nodes[-1]]))
            nodes.reverse()
            used = [
                self.index[self.ids[tuple(sorted(e))]] for e in zip(nodes, nodes[1:])
            ]
            paths.append(dict(stations=nodes))
            for name, changes in edge_derivative.items():
                derivatives[name][i, j] = derivatives[name][j, i] = changes[used].sum()
        return dict(
            order=self.order,
            edges=len(self.index),
            effect_types=len(self.menu),
            minimum_link=weights.min(),
            distances_with_held_fifth=d.tolist(),
            four_gram=np.linalg.eigvalsh(gram(d[:4, :4])).tolist(),
            five_gram=np.linalg.eigvalsh(gram(d)).tolist(),
            paths=paths,
            slopes={k: m.tolist() for k, m in derivatives.items()},
        )


def tetrahedron(ruler, state, rotation, tangents, panel):
    station = signed_transfers.station

    def read(x):
        d, _, _ = ruler.distances(x)
        return np.array([*simplex_area_volume(d[:4, :4]), station.entropy(x)])

    base = read(state)
    d = np.array(panel["distances_with_held_fifth"])[:4, :4]
    response = {}
    for k, v in panel["slopes"].items():
        response[k] = np.imag(simplex_area_volume(d + 1e-16j * np.array(v))) / 1e-16
    checked = {}
    for k, raw in tangents.items():
        alpha = -response[k][1] / response["rotation"][1]
        area = (response[k] + alpha * response["rotation"])[0]
        entropy = -station.pair(raw[None], station.log_blocks(state))[0]
        finite = []
        for step in (1e-4, 1e-5, 1e-6):
            secant = []
            for sign in (-1, 1):
                t = sign * step
                trial = state + t * raw

                def rotated(theta):
                    u = signed_transfers.rotation(rotation, theta)
                    return u @ trial @ u.conj().T

                # A narrow local bracket can avoid a distinct root for large
                # t. Include any failed bracket explicitly in the receipt.
                def error(theta):
                    return read(rotated(theta))[1] - base[1]

                lo, hi = sorted((0.0, 2 * alpha * t))
                if error(lo) * error(hi) > 0:
                    finite.append(dict(step=step, sign=sign, bracketed=False))
                    break
                theta = brentq(error, lo, hi, xtol=1e-14)
                y = rotated(theta)
                measures = read(y)
                secant.append(measures)
                finite.append(
                    dict(
                        step=step,
                        sign=sign,
                        bracketed=True,
                        correction_angle=theta,
                        volume_error=measures[1] - base[1],
                    )
                )
            if len(secant) == 2:
                finite.append(
                    dict(
                        step=step,
                        centered_area_slope=(secant[1][0] - secant[0][0]) / (2 * step),
                        centered_entropy_slope=(secant[1][2] - secant[0][2])
                        / (2 * step),
                    )
                )
        checked[k] = dict(
            angle=alpha,
            area_slope=area,
            entropy_slope=entropy,
            individual_gain=-entropy / area,
            finite=finite,
        )
    return dict(
        baseline=base.tolist(),
        ordinary=checked["ordinary"],
        pairing=checked["unpair"],
        pairing_transfer_at_ordinary_gain=checked["ordinary"]["individual_gain"]
        / checked["unpair"]["individual_gain"],
    )


def run():
    states = optics.source()
    x = states.pop("state")
    receipt = json.loads((HERE.parent / "signed_transfers.json").read_text())
    r = signed_transfers.right_coin(x)
    pauli = signed_transfers.phase.AXES[1:]
    bloch = np.array([np.trace(r @ p).real for p in pauli])
    axis = sum(v * p for v, p in zip(bloch, pauli)) / np.linalg.norm(bloch)
    h = np.kron(np.eye(2), axis / 2)
    rotation = -1j * (h @ x - x @ h)
    all_tangents = dict(states, rotation=rotation)
    for row in receipt["tests"]:
        key = "ordinary" if row["source"] == "ordinary" else "unpair"
        all_tangents["cap_six_matched_" + key] = (
            states[key] + row["angle_per_unit"] * rotation
        )
    board = optics.pairing.formation.federation(1)
    route = optics.pairing.formation.region(board)
    data = []
    for order in (1, 2, 3, 4):
        ruler = Ruler(board, route, order)
        p = ruler.observe(x, all_tangents)
        p["blank_reference"] = ruler.blank_reference(board)
        if order == 4:
            p["own_volume"] = tetrahedron(ruler, x, axis, states, p)
        data.append(p)
    return dict(ports=route, held_port=5, cells=80, panels=data)


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    body = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if cli.parse_args().write:
        (HERE / "results.json").write_text(body)
    else:
        print(body, end="")
