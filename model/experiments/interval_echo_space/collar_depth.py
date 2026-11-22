"""Follow one mobile load and two distinct retained repair coins on the complete sphere."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.row_30_arrows.experiment import (
    local_check,
)
from .shell_source import Board, check, laplacian, local_source
from .shells import hands, protocol

HERE = Path(__file__).resolve().parent
BOOK = HERE / "collar_depth_protocol.json"


def lifted_operator(board, transfer):
    """Four coin words, each with one copy of the port space. Remove common vacuum energies."""
    c = board.car
    edge = min(e for e in board.active if 4 in e and e[1] < 12)
    sites = [edge, tuple(v + 12 * (c.cells // 2) for v in edge)]
    check(
        all(e in board.active for e in sites),
        "registered collars must be blank active contacts",
    )
    k = candidate.KAPPA
    h, _, shift = board.generator("native")
    spatial = h - shift * sparse.eye(c.n, format="csr") - laplacian(c.n, sites) / 2
    diagonal = np.zeros(c.n)
    for a, b in sites:
        diagonal[[a, b]] -= k
    spatial += sparse.diags(diagonal)
    full = sparse.kron(sparse.eye(4), spatial, format="csr")
    row, col, value = [], [], []
    for j, (a, b) in enumerate(sites):
        bit = 1 << j
        for mask in range(4):
            if mask & bit:
                continue
            left, right = mask * c.n + b, (mask | bit) * c.n + a
            row.extend([left, right])
            col.extend([left, right])
            value.extend([-k / 2, -k / 2])
            if transfer:
                row.extend([left, right])
                col.extend([right, left])
                value.extend([-k / 2, -k / 2])
    full += sparse.csr_matrix((value, (row, col)), shape=full.shape)
    check((full - full.T).nnz == 0, "non-Hermitian native restriction")
    return full, sites


def run(output):
    spec = json.loads(BOOK.read_text())
    matching = protocol("shell_protocol.json")["matching"]
    mark_error, occupation_error = local_check()
    hands_term, _ = local_source(1.0, candidate.KAPPA)
    clock = hands(hands_term, np.array(spec["appointments"]))
    result = dict(
        protocol_sha256=sha256(BOOK.read_bytes()).hexdigest(),
        native_mark_error=mark_error,
        occupation_read_error=occupation_error,
        boards=[],
    )
    result["source_sha256"] = {
        name: sha256((HERE / name).read_bytes()).hexdigest()
        for name in (
            "collar_depth.py",
            "shell_sources.json",
            "shell_source.py",
            "shells.py",
            "../../candidate_laws/interval_echo/openwave/row_30_arrows/experiment.py",
        )
    }
    saved = {}
    for level in (0, 1, 2):
        board = Board(level, matching)
        n = board.car.n
        depth = np.array([w.bit_count() for w in range(4)])
        experiments = []
        for transfer in (True, False):
            h, sites = lifted_operator(board, transfer)
            starts = [(i, word) for i in (0, 1) for word in (1, 2)]
            inputs = [word * n + sites[i][1] for i, word in starts]
            wave = np.zeros((4 * n, len(inputs)), complex)
            wave[inputs, np.arange(len(inputs))] = 1
            curvature = (
                (abs(h @ wave) ** 2).reshape(4, n, len(inputs))[depth == 2]
            ).sum(axis=(0, 1))
            for j, (i, word) in enumerate(starts):
                expected = (
                    (candidate.KAPPA / 2) ** 2
                    if transfer and not (word & (1 << i))
                    else 0.0
                )
                check(
                    abs(curvature[j] - expected) < 1e-12,
                    "depth gain disagrees with the native column",
                )
            previous = 0.0
            for index, t in enumerate(spec["appointments"]):
                action = (-1j * (t - previous)) * h
                wave = expm_multiply(action, wave, traceA=action.diagonal().sum())
                previous = t
                for j, (i, word) in enumerate(starts):
                    p = (abs(wave[:, j]) ** 2).reshape(4, n)
                    check(
                        abs(p.sum() - 1) < 1e-11, "incomplete mobile and memory records"
                    )
                    prob = (p[:, :, None] * clock[index]).reshape(-1)
                    rng = np.random.default_rng(
                        [spec["seed"], level, int(transfer), i, word, index]
                    )
                    counts = rng.multinomial(spec["shots"], prob / prob.sum()).reshape(
                        4, n, 4
                    )
                    saved[
                        f"l{level}_t{index}_site{i}_word{word}_transfer{int(transfer)}"
                    ] = counts
                    n2 = int(counts[3].sum())
                    check(
                        transfer or (p[3].sum() == 0 and n2 == 0),
                        "no-flip control moved in depth",
                    )
                    experiments.append(
                        dict(
                            transfer=transfer,
                            launch_site=i,
                            word=word,
                            hand_radians=t,
                            initial_depth_two_coefficient=float(curvature[j]),
                            depth_two_probability=float(p[3].sum()),
                            depth_two_count=n2,
                            complete_counts=int(counts.sum()),
                        )
                    )
        result["boards"].append(
            dict(cells=board.car.cells, marks=sites, experiments=experiments)
        )
    output.mkdir(parents=True, exist_ok=True)
    path = output / "collar_depth_counts.npz"
    np.savez_compressed(path, **saved)
    result["counts_sha256"] = sha256(path.read_bytes()).hexdigest()
    (output / "collar_depth.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args().output)
