"""Measure isolated affordable formation gaps against each strand's static vacuum."""

from hashlib import sha256
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from . import preparation, weak_strand as ws
from .measure import HERE
from ..passive_decay.measure import compare_records, hamiltonians

SHOTS = 2**22
SEED = 26092709


def acquire():
    rng = np.random.default_rng(SEED)
    core = eigh(hamiltonians(4)[0], eigvals_only=True)
    bound = core[-1] - core[2]
    designs = []
    for n in (6, 8):
        b, d, heads, _ = ws.ring_layout(n, True, False)
        full, _ = ws.constructed(d, b, heads, "source_cut")
        # Hold the first three logical bits fixed and remove their scalar.
        # With C,D physically cut this is the entire independent strand.
        h = full[: 2**n, : 2**n].tocsr()
        e, v = eigh(h.toarray())
        lines = np.arange(1, 1 + np.argmax(e - e[0] > bound))
        links, _ = preparation.mount(b, d, heads)
        edges = {
            (i - 3, j - 3): value for (i, j), value in links.items() if min(i, j) >= 3
        }
        designs.append((n, h, e, v, lines, edges))
    holds = 4.0 ** np.arange(7) / 64
    total = 2 * len(holds) * sum(len(d[4]) for d in designs)
    # Simultaneous Hoeffding bound on both coordinates at all unwrapping holds.
    radius = (4 * np.log(2 * total / 0.01) / SHOTS) ** 0.5
    rows = {}
    for n, h, e, v, lines, edges in designs:
        results = []
        for k in lines:
            ray = (v[:, 0] + v[:, k]) / np.sqrt(2)
            error, counts = preparation.synthesize(ray, edges)
            previous = 0
            points = []
            for t in holds:
                state = ws.propagate(h, ray, t)
                g, z = v[:, 0] @ state, v[:, k] @ state
                obs = []
                for phase in (1, 1j):
                    probability = abs(g + phase * z) ** 2 / 2
                    obs.append(int(rng.binomial(SHOTS, probability)))
                x, y = 2 * np.asarray(obs) / SHOTS - 1
                angle = np.arctan2(y, x)
                angle += 2 * np.pi * round((4 * previous - angle) / (2 * np.pi))
                previous = angle
                points.append(obs)
            line = float(previous / holds[-1])
            unc = float(np.arcsin(radius) / holds[-1])
            assert abs(line - (e[k] - e[0])) < unc
            results.append(
                dict(
                    level=int(k),
                    gap=line,
                    radius=unc,
                    preparation_error=float(error),
                    gates=sum(counts.values()),
                    counts=points,
                )
            )
        rows[str(n)] = results
    return dict(
        protocol_sha256=sha256((HERE / "STRAND_LINES.md").read_bytes()).hexdigest(),
        seed=SEED,
        shots_per_quadrature=SHOTS,
        holds=holds.tolist(),
        compared=rows,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--acquire", type=Path)
    args = p.parse_args()
    result = acquire()
    for n, rows in result["compared"].items():
        print(n, [(r["gap"], r["radius"]) for r in rows])
    if args.acquire:
        with args.acquire.open("x") as f:
            json.dump(result, f, allow_nan=False, separators=(",", ":"))
            f.write("\n")
    else:
        compare_records(json.loads((HERE / "strand_lines.json").read_text()), result)
