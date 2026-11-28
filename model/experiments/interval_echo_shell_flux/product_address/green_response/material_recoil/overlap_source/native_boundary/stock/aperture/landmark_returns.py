"""Does the unweighted native first return generate the entropy-selected source?"""

import argparse
from fractions import Fraction as Q
import json

import numpy as np
from scipy.sparse.csgraph import connected_components, shortest_path

from . import collective_retention as retention

HERE = retention.HERE


def ticket_tables(pulse):
    tables = []
    for mat in pulse:
        M = mat.tocsc()
        rows = [
            np.repeat(
                M.indices[M.indptr[i] : M.indptr[i + 1]],
                M.data[M.indptr[i] : M.indptr[i + 1]],
            )
            for i in range(M.shape[1])
        ]
        tables.append(np.stack(rows))
    return tables


def run():
    full, initial, _, service, _ = retention.setup()
    subset = np.flatnonzero(initial > 0)
    live = initial > 0
    n = len(full)
    cells = (retention.tidal.lab.energy.pair.stock.MAG[full[:, 4:6]] == 2).astype(int)
    labels = 2 * cells[:, 0] + cells[:, 1]
    pulse = [m for m, _ in service]
    size = 128**3
    cycles = {}
    for phase in range(3):
        sweep = pulse[(phase + 2) % 3] @ pulse[(phase + 1) % 3] @ pulse[phase]
        assert connected_components(sweep, directed=True, connection="strong")[0] == 1
        assert np.all(np.asarray(sweep.sum(0)) == size) and np.all(
            np.asarray(sweep.sum(1)) == size
        )
        # Column convention. Inflow to T from one uniform count per S
        # state equals the defect of one count per T state under the
        # killed kernel. The series of first returns therefore telescopes.
        ST = sweep[live][:, ~live]
        TT = sweep[~live][:, ~live]
        TS = sweep[~live][:, live]
        SS = sweep[live][:, live]
        assert np.all(TS.sum(1) + TT.sum(1) == size)
        assert np.all(SS.sum(1) + ST.sum(1) == size)
        dist = shortest_path(sweep, directed=True, unweighted=True, indices=subset).min(
            0
        )
        assert np.isfinite(dist).all()
        # Recurrence: every outside vertex can reach S with at least one
        # ticket per step in at most D sweeps.
        cycles[str(phase)] = dict(
            maximum_steps_to_enrollment=int(dist.max()),
            exact_return_blank=str(Q(int(sum(full[subset, -1] == 0)), len(subset))),
            exact_joint=list(
                map(
                    str,
                    [Q(int(sum(labels[subset] == k)), len(subset)) for k in range(4)],
                )
            ),
        )
    rng = np.random.default_rng(485127)
    tables = ticket_tables(pulse)
    trials = 2**18
    sites = rng.choice(subset, size=trials)
    arrived = np.zeros(n, int)
    steps = 0
    for steps in range(1, 4097):
        for lut in tables:
            sites = lut[sites, rng.integers(128, size=len(sites))]
        back = live[sites]
        arrived += np.bincount(sites[back], minlength=n)
        sites = sites[~back]
        if not len(sites):
            break
    eps = np.sqrt(np.log(2 * 5 / 0.01) / (2 * trials))
    counts = {}
    for feature, flag in [
        ("blank", full[:, -1] == 0),
        *[(f"cell_{k}", labels == k) for k in range(4)],
    ]:
        hits = int(arrived @ flag)
        band = [max(0, hits / trials - eps), min(1, (hits + len(sites)) / trials + eps)]
        assert band[0] <= initial @ flag <= band[1]
        counts[feature] = dict(seen=hits, bounds=band)
    return dict(
        cells=n,
        enrolled=len(subset),
        cycles=cycles,
        trials=trials,
        first_returns=int(arrived.sum()),
        unreturned=len(sites),
        last_sweep=steps,
        joint_coverage=0.99,
        held=counts,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    out = run()
    raw = json.dumps(out, indent=2, allow_nan=False) + "\n"
    path = HERE / "landmark_returns.json"
    if args.write:
        path.write_text(raw)
    else:
        assert path.read_text() == raw
    print(out["first_returns"], out["last_sweep"], out["cycles"]["0"])
