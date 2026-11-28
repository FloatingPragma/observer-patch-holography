"""Enumerate separated conversions against the parent's other registered heavy remnant."""

import argparse
import json
from hashlib import sha256
from pathlib import Path

from scipy.linalg import eigh

from . import joint_exterior as joint
from .measure import HERE
from ..passive_decay.measure import compare_records, hamiltonians


def acquire():
    exterior_path = HERE / "joint_exterior.json"
    panels = json.loads(exterior_path.read_text())["boards"]
    core = eigh(hamiltonians(4)[0], eigvals_only=True)
    budget = float(core[-1] - core[1])
    e, basis = eigh(joint.strand.strand(8).toarray())
    gaps = e - e[0]
    use = gaps <= budget
    levels = gaps[use]
    rates = {}
    for name, arm in joint.exterior.source_stocks(
        *joint.strand.ring_layout(8, True, False)[:3],
        basis,
        use,
        core_levels=slice(1, 2),
    ).items():
        if name != "known_D_phase":
            rates[name] = arm.diagonal().real.tolist()

    boards = {}
    for depth in range(3):
        b = joint.exterior.federation(depth)
        info = panels[str(b.cells)]
        rec = info["exterior"]["with_approach"]
        sites = sorted(
            {p for p in range(b.n) if p // 12 not in info["reserved_owners"]}
            | set(info["path"][1:-1])
        )
        h, _ = joint.single(b, sites)
        # The remaining candidate list is complete only if the independent
        # inertia difference equals the number of listed excited eigenvalues.
        count = joint.inertia(h, rec["ground"] + budget + 1e-7)[0] - 1
        remote = [x for x in rec["first_gaps"] if 1e-7 < x <= budget]
        assert len(remote) == count
        pairs = [
            [i, k, j, float(levels[i] - levels[k] - x)]
            for i in range(1, len(levels))
            for k in range(i)
            for j, x in enumerate(remote)
        ]
        best = min(pairs, key=lambda row: abs(row[-1]), default=None)
        boards[str(b.cells)] = dict(exterior_gaps=remote, tried=pairs, closest=best)
    return dict(
        protocol_sha256=sha256(
            (HERE / "LARGER_TRANSITION.md").read_bytes()
        ).hexdigest(),
        joint_archive_sha256=sha256(exterior_path.read_bytes()).hexdigest(),
        larger_parent_loss=budget,
        strand_levels=levels.tolist(),
        source_lower_remnant_and_level=rates,
        boards=boards,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquire", type=Path)
    args = parser.parse_args()
    data = acquire()
    for n, row in data["boards"].items():
        print(
            n,
            "exterior",
            row["exterior_gaps"],
            "closest index, residual, exterior, signed mismatch",
            row["closest"],
        )
    if args.acquire:
        with args.acquire.open("x") as out:
            json.dump(data, out, allow_nan=False, separators=(",", ":"))
            out.write("\n")
    else:
        compare_records(json.loads((HERE / "larger_transition.json").read_text()), data)
