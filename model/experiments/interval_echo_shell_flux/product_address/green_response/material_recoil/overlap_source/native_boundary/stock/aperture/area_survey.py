"""Test spherical area growth in the native swap-distance survey of the Poisson source."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json

from . import exchange_geometry as lab

HERE = lab.HERE


def ratios(areas, r):
    result = []
    for family in areas:
        a, b = sum(family[r : 2 * r]), sum(family[2 * r : 4 * r])
        if not a:
            if b:
                return None
            continue
        result.append(Q(b, 2 * a))
    return result


def run():
    source = HERE / "poisson_response.json"
    data = json.loads(source.read_bytes())
    lab.contact_swap()
    charts = []
    for level, record in enumerate(data["boards"]):
        board = lab.lab.lab.cold.federation(level)
        d = lab.lab.lab.cold.distances(board, record["written"])
        local = [0] * int(max(d))
        external = local.copy()
        for u, v in board.seams:
            if d[u] != d[v]:
                (local if u // 12 == v // 12 else external)[min(d[u], d[v])] += 1
        assert [x + y for x, y in zip(local, external, strict=True)] == record[
            "cut_edges"
        ]
        rows = []
        r = 1
        while 4 * r <= len(local):
            growth = ratios([local, external], r)
            # Permit the origin anywhere between the two written endpoints and
            # independently enlarge each inner radius by one hop for this lower bound.
            euclid = Q(
                sum(j * j for j in range(2 * r, 4 * r)),
                2 * sum((j + 1) ** 2 for j in range(r, 2 * r)),
            )
            total = ratios([[x + y for x, y in zip(local, external, strict=True)]], r)[
                0
            ]
            rows.append(
                dict(
                    first=r,
                    inner_to_outer_mean_area_ratio=str(total),
                    arbitrary_contact_areas=None
                    if growth is None
                    else [str(min(growth)), str(max(growth))],
                    square_law_lower=str(euclid),
                    excluded_for_both_types=growth is not None and max(growth) < euclid,
                    enclosed_ports=sum(record["sites"][: 4 * r]),
                )
            )
            r *= 2
        charts.append(
            dict(
                ports=board.n,
                within_owner=local,
                between_owners=external,
                comparisons=rows,
            )
        )
    return dict(source_sha256=sha256(source.read_bytes()).hexdigest(), charts=charts)


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    out = run()
    path = HERE / "area_survey.json"
    if args.write:
        path.write_text(json.dumps(out, indent=2) + "\n")
    else:
        assert json.loads(path.read_bytes()) == out
    print(out["charts"][-1]["comparisons"])
