"""Spectra and owner-distance flights of the lowest isolated native neutral-pair cluster."""

import argparse
import json
from collections import Counter
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy import linalg, sparse
from scipy.sparse.csgraph import connected_components, shortest_path

from model.executable_model import federation
from model.executable_model.carrier import cell
from .. import neutral_defect as defect
from ..pair_recoil.measure import PLACEMENTS

HERE = Path(__file__).resolve().parent
CENTER = 36.115636
CUT = 36.1164
SHOTS = 2**26
TIMES = (0.0, 100000.0, 300000.0, 1000000.0)


def motifs(board):
    internal = shortest_path(cell().adjacency(), unweighted=True)
    labels = []
    for owner in range(board.cells):
        gateways = np.flatnonzero(board.degree[12 * owner : 12 * (owner + 1)] == 6)
        labels.append(
            "".join(
                map(
                    str,
                    sorted(int(internal[i, j]) for i, j in combinations(gateways, 2)),
                )
            )
        )
    return labels


def active_geometry(board):
    labels = motifs(board)
    active = np.flatnonzero(np.asarray(labels) == "222")
    edges = np.asarray([(a // 12, b // 12) for a, b in board.external])
    row, col = edges.T
    dual = sparse.coo_matrix(
        (np.ones(len(row) * 2), (np.r_[row, col], np.r_[col, row])),
        shape=(board.cells, board.cells),
    ).tocsr()
    sub = dual[active][:, active]
    count, ids = connected_components(sub, directed=False)
    parts = []
    lengths = shortest_path(sub, unweighted=True)
    for i in range(count):
        members = np.flatnonzero(ids == i)
        parts.append(
            dict(
                cells=len(members),
                diameter=int(lengths[np.ix_(members, members)].max()),
            )
        )
    return dict(motifs=dict(Counter(labels)), active_components=parts)


class PairReduction:
    """Eliminate complementary configurations with a certified spectral gap."""

    def __init__(self, h, cost):
        self.high = np.flatnonzero(cost == 32)
        self.other = np.flatnonzero(cost != 32)
        self.hh = h[self.high][:, self.high].toarray()
        self.qq = h[self.other][:, self.other].tocsr()
        self.b = h[self.other][:, self.high].tocsr()
        self.diagonal = self.qq.diagonal()
        self.off = self.qq - sparse.diags(self.diagonal)
        self.b_bound = np.sqrt(
            abs(self.b).sum(axis=0).max() * abs(self.b).sum(axis=1).max()
        )
        self.q_upper = float(
            np.max(self.diagonal + np.asarray(abs(self.off).sum(axis=1)).ravel())
        )
        if self.q_upper > 28:
            raise ValueError("Complementary block exceeds the native gap bound")

    def schur(self, energy):
        scale = sparse.diags(1 / (energy - self.diagonal))
        step = (scale @ self.off).tocsr()
        term = (scale @ self.b).tocsr()
        resolvent = term.copy()
        for _ in range(40):
            term = step @ term
            term.data[abs(term.data) < 2e-17] = 0
            term.eliminate_zeros()
            if not term.nnz:
                break
            resolvent = resolvent + term
        residual = self.b - (energy * resolvent - self.qq @ resolvent)
        error = float(
            self.b_bound * sparse.linalg.norm(residual) / (energy - self.q_upper)
        )
        estimate = self.hh + (self.b.T @ resolvent).toarray()
        return (estimate + estimate.T) / 2, resolvent, error

    def inertia(self, energy):
        schur, _, error = self.schur(energy)
        values = linalg.eigvalsh(schur) - energy
        clearance = float(np.min(abs(values)))
        if clearance < error + 2e-12:
            raise ValueError("Unresolved spectral threshold")
        return dict(
            energy=float(energy),
            upper_lines_below=int(np.sum(values < 0)),
            schur_error_bound=error,
            clearance=clearance,
        )

    def lowest(self, h, count):
        schur, lift, error = self.schur(CENTER)
        gram = np.eye(len(self.high)) + (lift.T @ lift).toarray()
        shifts, vectors = linalg.eigh(
            schur - CENTER * np.eye(len(gram)), gram, subset_by_index=(0, count - 1)
        )
        first = lift @ vectors
        scale = (1 / (CENTER - self.diagonal))[:, None]
        term = scale * first
        second = term.copy()
        for _ in range(32):
            term = scale * (self.off @ term)
            second += term
        # Taylor-correct the energy dependence in the eliminated amplitudes.
        basis = np.empty((h.shape[0], count))
        basis[self.high] = vectors
        basis[self.other] = first - second * shifts
        projected = basis.T @ (h @ basis)
        energies, mixing = linalg.eigh((projected + projected.T) / 2, basis.T @ basis)
        basis = basis @ mixing
        residual = h @ basis - basis * energies
        return energies, basis, float(linalg.norm(residual)), error


def transport(board, x, y, cost, basis, energies, residual, projector_error):
    edges = np.array([(a // 12, b // 12) for a, b in board.external])
    row, col = edges.T
    dual = sparse.coo_matrix(
        (np.ones(len(row) * 2), (np.r_[row, col], np.r_[col, row])),
        shape=(board.cells, board.cells),
    )
    distances = shortest_path(dual, unweighted=True).astype(int)
    x_owner, y_owner = x // 12, y // 12
    flights = []
    for seam in PLACEMENTS:
        for side, gateway in enumerate(board.external[seam]):
            owner = gateway // 12
            regular = min(
                j
                for j in board.neighbours[gateway]
                if j // 12 == owner and board.degree[j] == 5
            )
            launch = defect.ordered_index(regular, gateway, board.n)
            overlap = basis[launch]
            admission = float(overlap @ overlap)
            receipt = dict(
                seam=seam,
                side=side,
                owner=owner,
                ports=[regular, gateway],
                admission=admission,
                appointments=[],
            )
            flights.append(receipt)
            if admission < 0.001:
                continue
            coefficients = overlap / np.sqrt(admission)
            dist = distances[owner]
            diameter = int(dist.max())
            bins = np.where(x_owner == y_owner, dist[x_owner], diameter + 1)
            for it, t in enumerate(TIMES):
                state = basis @ (
                    np.exp(-1j * (energies - energies[0]) * t) * coefficients
                )
                probs = np.bincount(
                    bins, weights=abs(state) ** 2, minlength=diameter + 2
                )
                np.testing.assert_allclose(probs.sum(), 1, rtol=0, atol=5e-12)
                counts = np.random.default_rng([770412, seam, side, it]).multinomial(
                    SHOTS, probs / probs.sum()
                )
                # Norm distance bounds any individual bin or union of bins for normalized pure states.
                prob_error = (
                    2 * projector_error / np.sqrt(admission) + t * residual + 1e-10
                )
                receipt["appointments"].append(
                    dict(
                        t=t,
                        owner_distance_then_separated=probs.tolist(),
                        clicks=counts.tolist(),
                        probability_error_bound=prob_error,
                        shared_owner_two_or_more=float(probs[2:-1].sum()),
                    )
                )
                bound_probs = np.bincount(
                    bins[cost == 32],
                    weights=abs(state[cost == 32]) ** 2,
                    minlength=diameter + 1,
                )
                bound_probs = np.r_[bound_probs, 1 - bound_probs.sum()]
                bound_clicks = np.random.default_rng(
                    [770413, seam, side, it]
                ).multinomial(SHOTS, bound_probs)
                receipt["appointments"][-1]["bound"] = dict(
                    upper_address_distance_then_other=bound_probs.tolist(),
                    clicks=bound_clicks.tolist(),
                    two_or_more=float(bound_probs[2:-1].sum()),
                )
            baseline = receipt["appointments"][0]
            for a in receipt["appointments"][1:]:
                for first_read, last_read in (
                    (baseline, a),
                    (baseline["bound"], a["bound"]),
                ):
                    rates = np.array(
                        [
                            sum(k["clicks"][2:-1]) / SHOTS
                            for k in (first_read, last_read)
                        ]
                    )
                    se = float(np.sqrt((rates * (1 - rates)).sum() / SHOTS))
                    gain = float(rates[1] - rates[0])
                    last_read.update(
                        two_hop_gain=gain,
                        sampling_standard_error=se,
                        resolved=bool(
                            gain
                            > 5 * se
                            + baseline["probability_error_bound"]
                            + a["probability_error_bound"]
                        ),
                    )
    return flights


def measure():
    board = federation(0)
    x, y, _, _, h, cost = defect.generators(board)
    reduction = PairReduction(h, cost)
    cut = reduction.inertia(CUT)
    energies, basis, residual, schur_error = reduction.lowest(
        h, cut["upper_lines_below"]
    )
    np.testing.assert_allclose(
        basis.T @ basis, np.eye(len(energies)), rtol=0, atol=1e-12
    )
    margin = max(2e-9, 10 * residual)
    counts = [
        reduction.inertia(energies[0] - margin),
        reduction.inertia(energies[-1] + margin),
        cut,
        reduction.inertia(36.11715),
        reduction.inertia(36.11720),
    ]
    if [c["upper_lines_below"] for c in counts] != [0, 24, 24, 24, 32]:
        raise ValueError(
            "The independently counted subspaces differ from the enrolled cluster"
        )
    projector_error = residual / min(CUT - energies[-1], energies[0] - 30)
    high = reduction.high
    high_owner = x[high] // 12
    weight = np.bincount(
        high_owner, weights=np.sum(basis[high] ** 2, axis=1), minlength=board.cells
    )
    labels = motifs(board)
    bind = linalg.eigvalsh(basis[high].T @ basis[high])[0]
    return dict(
        carrier=board.name,
        ordered_states=len(x),
        native_high_rank=len(high),
        q_upper=reduction.q_upper,
        lines=energies.tolist(),
        full_generator_residual_frobenius=residual,
        calibrated_schur_error_bound=schur_error,
        cluster_projector_error=projector_error,
        threshold_counts=counts,
        high_address_binding_floor=float(bind - 2 * projector_error),
        sharp_spectral_weight_per_cell=weight.tolist(),
        gateway_motifs=labels,
        all_refinement_geometry=[active_geometry(federation(i)) for i in range(3)],
        shots=SHOTS,
        flights=transport(
            board, x, y, cost, basis, energies, residual, projector_error
        ),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = measure()
    payload = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(payload)
    else:
        print(payload, end="")
