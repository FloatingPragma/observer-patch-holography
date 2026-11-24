"""Select a disjoint pair using sector-ground variances and acquire held pulse records."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply
from scipy.stats import fisher_exact

from model.executable_model import federation
from model.experiments.material_inertia.neutral_defect import generators, ordered_index
from .measure import HERE, loop_minus, vacuum


def run():
    board = federation(0)
    x, y, _, _, h, _ = generators(board)
    energy, w = vacuum(h, board, x, y)
    at_gate = {}
    seen = set()
    for face in board.faces():
        for root in face:
            a, b = [v for v in face if v != root]
            label = (root, a, b)
            prob = (
                w[ordered_index(a, root, board.n)] - w[ordered_index(b, root, board.n)]
            ) ** 2 / 2
            for gate in face:
                if board.degree[gate] == 6:
                    at_gate.setdefault(gate, []).append((label, prob * (1 - prob)))
            seen.add(label)
    comparisons = []
    for u, v in board.external:
        for start, end in ((u, v), (v, u)):
            for a, va in at_gate[start]:
                for b, vb in at_gate[end]:
                    comparisons.append((-min(va, vb), a, b))
    negscore, home, remote = min(comparisons)
    hit = loop_minus(w, board, *home)
    p = float(hit @ hit)
    local = hit - p * w
    detector = loop_minus(w, board, *remote)
    reference = float(detector @ detector)
    generator = h - energy * sparse.eye(h.shape[0], format="csr")
    times = (0, 0.5, 1, 2, 4, 8)
    rng = np.random.default_rng(624702)
    shots = 2**30
    rows = []
    for t in times:
        state = w + expm_multiply((-1j * t) * generator, -2 * hit) if t else w - 2 * hit
        yes = float(np.linalg.norm(loop_minus(state, board, *remote)) ** 2)
        ns, ny = map(int, rng.binomial(shots, [reference, yes]))
        # SciPy's hypergeometric implementation uses integer arithmetic or double
        # precision according to the chosen regime. Zeros remain exact outcomes.
        test = fisher_exact([[ns, shots - ns], [ny, shots - ny]]).pvalue
        corr = (
            float(local @ expm_multiply(-t * generator, local))
            if t
            else float(local @ local)
        )
        rows.append(
            dict(
                t=t,
                sham=reference,
                pulse=yes,
                counts=[ns, ny],
                two_sided_p=float(test),
                resolved=bool(test < 0.01 / len(times)),
                home_connected=corr,
            )
        )
    return dict(
        protocol_sha256=hashlib.sha256(
            (HERE / "PLACEMENT.md").read_bytes()
        ).hexdigest(),
        rooted_triangles=len(seen),
        eligible_pairs=len(comparisons),
        home=home,
        remote=remote,
        smaller_variance=float(-negscore),
        home_mean=p,
        shots=shots,
        seed=624702,
        rows=rows,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with args.output.open("x") as f:
        json.dump(run(), f, indent=2, allow_nan=False)
        f.write("\n")
