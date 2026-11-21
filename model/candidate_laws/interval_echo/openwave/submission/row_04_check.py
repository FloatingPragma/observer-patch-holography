"""Independent sparse-exponential check of the freely released Interval Echo pair."""

import argparse
import base64
import hashlib
import itertools
import json
from pathlib import Path
import zlib

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply


HERE = Path(__file__).resolve().parent
BOARD_HASH = "99870398e3d516f899c600e9bf25ba53c7f513142ba8facbe1cedd94b23f915b"
EXPOSURES = (1, 7, 37)
SHOTS = 4096
MARGIN = np.sqrt(np.log(8 / 0.01) / (2 * SHOTS))


def screen(path):
    if hashlib.sha256(path.read_bytes()).hexdigest() != BOARD_HASH:
        raise ValueError("Unexpected 20-cell carrier file")
    with np.load(path, allow_pickle=False) as board:
        edges = np.sort(board["edges"], axis=1)
    edges = edges[np.lexsort((edges[:, 1], edges[:, 0]))]
    n = int(edges.max()) + 1
    if (n, len(edges)) != (240, 630):
        raise ValueError("The registered fixture has 240 ports and 630 contacts")
    ports = np.arange(n)
    adjacency = sparse.coo_matrix((np.ones(len(edges)), edges.T), shape=(n, n)).tocsr()
    adjacency += adjacency.T
    adjacency.sort_indices()
    pairs = np.array(list(itertools.combinations(ports, 2)))
    pair_index = np.full((n, n), -1)
    pair_index[pairs[:, 0], pairs[:, 1]] = np.arange(len(pairs))
    pair_index[pairs[:, 1], pairs[:, 0]] = np.arange(len(pairs))
    rows, columns = [], []
    for port in ports:
        neighbors = adjacency.indices[adjacency.indptr[port]:adjacency.indptr[port + 1]]
        for neighbor in neighbors:
            other = ports[(ports != port) & (ports != neighbor)]
            rows.extend(pair_index[port, other])
            columns.extend(pair_index[neighbor, other])
    moves = sparse.csr_matrix((np.ones(len(rows)), (rows, columns)), shape=(len(pairs),) * 2)
    exchange = (sparse.diags(np.asarray(moves.sum(axis=1)).ravel()) - moves) / 2
    adjacent = np.zeros(len(pairs))
    adjacent[pair_index[edges[:, 0], edges[:, 1]]] = 1
    return edges, pair_index, exchange, adjacent


def local_law(state, patch):
    probabilities = np.abs(state)**2
    if abs(probabilities.sum() - 1.0) > 1e-8:
        raise ArithmeticError("Free evolution lost norm")
    quiet = np.ones(len(state), dtype=bool)
    quiet[patch] = False
    return np.r_[probabilities[patch], probabilities[quiet].sum()]


def evolve(h, source, times):
    state = np.zeros(h.shape[0], complex)
    state[source] = 1
    previous = 0
    for exposure in times:
        state = expm_multiply(-1j * (exposure - previous) * h, state)
        yield state
        previous = exposure


def run(board, kappa):
    edges, index, exchange, adjacent = screen(board)
    old = edges[(edges[:, 0] < 12) & (edges[:, 1] < 12)]
    home = edges[(edges[:, 0] // 12 == 13) & (edges[:, 1] // 12 == 13)]
    gateways = np.unique(edges[edges[:, 0] // 12 != edges[:, 1] // 12])
    launch = home[np.isin(home, gateways).any(axis=1)][0]
    patch = index[home[:, 0], home[:, 1]]
    original = index[old[:, 0], old[:, 1]]
    # Decoder challenge uses genuine states in the same two-load sector.
    for state_index, outcome in [(patch[0], 0), (patch[1], 1), (index[0, 13], len(home))]:
        state = np.zeros(len(adjacent))
        state[state_index] = 1
        if local_law(state, patch)[outcome] != 1:
            raise ArithmeticError("Patch decoder missed a planted terminal outcome")
    report = {"home_cell": 13, "source": launch.tolist(), "addresses": home.tolist() + ["quiet"],
              "shots_per_setting": SHOTS, "native_kappa": kappa, "arms": []}
    for arm, strength in enumerate([kappa, 0.0]):
        h = exchange - 2 * strength * sparse.diags(adjacent)
        original_state = next(evolve(h, index[0, 1], [150]))
        original_probability = float(local_law(original_state, original)[:-1].sum())
        arm_result = {"strength": strength, "original_probability": round(original_probability, 9),
                      "readings": []}
        for j, state in enumerate(evolve(h, index[tuple(launch)], EXPOSURES)):
            p = local_law(state, patch)
            rng = np.random.default_rng([917421, arm, j])
            draws = rng.choice(len(p), size=SHOTS, p=p / p.sum()).astype(np.uint8)
            counts = np.bincount(draws, minlength=len(p))
            arm_result["readings"].append({"exposure": EXPOSURES[j], "counts": counts.tolist(),
                "probabilities": np.round(p, 9).tolist(),
                "trials_zlib_base64": base64.b64encode(zlib.compress(draws.tobytes(), 9)).decode()})
        report["arms"].append(arm_result)
    return report


def scientific_checks(report):
    native, control = report["arms"]
    hits = lambda setting: sum(setting["counts"][:-1]) / SHOTS
    if min(hits(s) for s in native["readings"]) - MARGIN <= 0.8:
        raise SystemExit("REJECT: freely released pair missed its local persistence gate")
    if hits(control["readings"][-1]) + MARGIN >= 0.65:
        raise SystemExit("REJECT: unchanged detector did not resolve loss under exchange alone")
    if max(abs(a["original_probability"] - p)
           for a, p in zip(report["arms"], [0.983881, 0.149559])) > 0.000002:
        raise SystemExit("Original headline could not be reproduced within published rounding")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--board", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, default=HERE / "row_04_receipt.json")
    parser.add_argument("--record", action="store_true", help="create a new receipt exclusively")
    parser.add_argument("--remove-contrast", action="store_true", help="exercise the real rejection path")
    args = parser.parse_args()
    result = run(args.board, 0.0 if args.remove_contrast else 8.0)
    if args.record:
        with args.receipt.open("x") as f:
            json.dump(result, f, ensure_ascii=False, indent=2)
            f.write("\n")
    scientific_checks(result)
    if result != json.loads(args.receipt.read_text()):
        raise SystemExit("Calculated records differ from the archived receipt")
    for arm in result["arms"]:
        print(arm["strength"], arm["original_probability"],
              [sum(setting["counts"][:-1]) for setting in arm["readings"]])
