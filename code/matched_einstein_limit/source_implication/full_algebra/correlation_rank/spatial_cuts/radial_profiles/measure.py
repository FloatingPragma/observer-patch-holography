"""Classify six-radius and cell-resolved weights of the native translation flux."""

import argparse
import json
from math import floor

import numpy as np
from mpmath import iv

from .. import measure as cuts
from . import effects, held, proofs

HERE = effects.HERE


def readings():
    raw = json.loads((HERE / "effects.json").read_text())
    if raw["scale"] != effects.SCALE:
        raise ValueError("radial transfer enclosure scale changed")
    by = {}
    for row in raw["rows"]:
        key = row["charge"], row["name"]
        if key in by:
            raise ValueError("duplicate regional or spatial reading")
        by[key] = np.array(
            [iv.mpf([lo, hi]) / raw["scale"] for lo, hi in row["row"]], dtype=object
        )
    cells = sorted(sum(raw["shells"], []))
    regional = np.array([by[q, n] for q in (-1, 1) for n in cuts.regional.NAMES])
    shells = np.array([[by[q, str(i)] for i in range(6)] for q in (-1, 1)])
    nodes = np.array([[by[q, f"cell-{u}"] for u in cells] for q in (-1, 1)])
    if len(raw["shells"]) != 6 or len(cells) != len(set(cells)) or len(cells) != 12:
        raise ValueError("unexpected chart partition")
    # The seventeen-constraint angular witness fixes all inner-cell transfers
    # by the indicated shell identities. Verify those exact memberships.
    if raw["shells"] != [[0], [3, 16], [2, 12, 19], [18], [28, 29, 66], [31, 70]]:
        raise ValueError("different chart cells require an independent angular witness")
    for charge in (-1, 1):
        for d, owners in enumerate(raw["shells"]):
            total = sum((by[charge, f"cell-{u}"] for u in owners))
            if any(0 not in term for term in total - by[charge, str(d)]):
                raise ValueError("cell transfer disagrees with shell intervention")
        boundary = by[charge, "domain_cut"] + by[charge, "image_cut"]
        if any(
            0 not in term for term in shells[int(charge == 1)].sum(axis=0) - boundary
        ):
            raise ValueError("chart Ward identity did not enclose the boundary sum")
    geo = json.loads((cuts.survey.HERE / "results.json").read_text())
    area, volume = [cuts.survey.unpack(geo[key]).ravel() for key in ("area", "volume")]
    return regional, shells, cells, nodes, area, volume


def angular_spread(equal, area, cells, nodes):
    # Protecting three more cells fixes both inner non-singleton shells
    # cell by cell: shell [3,16] and shell [2,12,19]. Only shells 4,5 can vary.
    net = nodes.sum(axis=0)
    protected = np.vstack((equal, net[[cells.index(i) for i in (12, 16, 19)]]))
    active = [cells.index(i) for i in (28, 29, 31, 66, 70)]
    outer_ent, changes = proofs.hidden_slope(protected, area, net[active])
    # Each active shell has zero total. Centre the prices at the
    # average of their maximum and minimum before bounding their sum.
    spread = 2 * outer_ent / sum(abs(x).b for x in changes)
    return outer_ent, floor(cuts.survey.dyadic(spread._mpi_[0]))


def certify():
    regional, shells, cells, nodes, area, volume = readings()
    equal = np.vstack((regional, volume, (shells[0] + shells[1])[1:]))
    ent, _ = proofs.hidden_slope(equal, area)
    outer_ent, minimum_spread = angular_spread(equal, area, cells, nodes)
    cones = [
        proofs.signed_cone(regional, shells, area, volume, signs)
        for signs in ((1, 1), (-1, -1), (1, -1), (-1, 1))
    ]
    if any(c is None for c in cones):
        raise ValueError("a monotone radial price pair lacks a separating direction")
    cell_reader = proofs.angular_profile(regional, cells, nodes)
    arms = held.run(regional, shells, cell_reader["regional_and_cell_prices"])
    return dict(
        scale=cuts.survey.SCALE,
        projective_charts=proofs.projective(regional, shells, area, volume),
        equal_probe_entropy_slope=cuts.survey.pack(ent),
        outer_angular_entropy_slope=cuts.survey.pack(outer_ent),
        minimum_within_shell_price_spread=minimum_spread,
        monotone_profiles=cones,
        independent_profiles=proofs.independent(regional, shells),
        outer_shells=proofs.outer_jump(regional, shells, area, volume),
        common_cell_profile=cell_reader,
        equal_probe_arms=arms,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    iv.dps = 54
    report = certify()
    if args.write:
        (HERE / "results.json").write_text(
            json.dumps(report, indent=2, allow_nan=False) + "\n"
        )
    print(json.dumps(report, allow_nan=False))
