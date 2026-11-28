"""Certify location and minimum reach of positive Poisson sources from observed rings."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
import math

from . import kernel_budget as budget
from . import newton_transport as transport

HERE = budget.HERE


def tail(distance):
    if distance <= 0:
        return Q(1)
    ratio = Q(180, 289)
    return min(Q(1), 2 * ratio**distance / (1 - ratio))


def reach(distance, needed, cap, min_degree, field=False):
    if cap < needed:
        return dict(inventory_insufficient=True)

    def upper(k):
        return min(Q(1), Q(7, 1 + min_degree) * tail(k - 1)) if field else tail(k)

    k = 0
    while cap * upper(k + 1) >= needed:
        k += 1
    assert cap * upper(k) >= needed
    return dict(
        inventory_insufficient=False,
        max_minimum_hops=k,
        length_floor_kpc=math.floor(distance * 10**6 / k) / 10**6 if k else None,
    )


def panel(rows, baryon_edge, baryon_mass, cap, min_degree):
    edge = Q(str(baryon_edge))
    # Printed processed inputs are verified as rational numbers.
    data = transport.mass_intervals(rows, baryon_mass)
    assert all(a[0] < b[0] for a, b in zip(data, data[1:]))
    exterior = [i for i, (r, _, _) in enumerate(data) if r >= edge]
    monotone = True
    prefix = Q(0)
    witness = None
    step = []
    for i, (_, lo, hi) in enumerate(data):
        prefix = max(prefix, lo)
        if prefix > hi and monotone:
            monotone = False
            witness = i
        step.append(prefix)
    ext_upper = min((data[i][2] for i in exterior), default=prefix)
    local = monotone and prefix <= ext_upper
    increases = []
    for i in exterior:
        r, _, hi = data[i]
        for j in exterior:
            outer, lo, _ = data[j]
            if outer <= r or lo <= hi:
                continue
            # For an inconsistent negative inner upper bound, the signed
            # increase is still valid, but a positive model has failed already.
            gap = lo - hi
            x = dict(
                first=i,
                last=j,
                radius_kpc=float(r),
                outer_kpc=float(outer),
                outside_map_kpc=math.floor((r - edge) * 10**6) / 10**6,
                minimum_net_ratio=math.floor(gap * 10**8) / 10**8,
            )
            if hi >= 0:
                for field in (False, True):
                    x["field" if field else "occupation"] = reach(
                        r - edge, gap, cap, min_degree, field
                    )
            increases.append(x)
    strongest = max(increases, key=lambda x: x["minimum_net_ratio"], default=None)
    farthest = max(increases, key=lambda x: x["outside_map_kpc"], default=None)
    tight = {}
    for k in ("occupation", "field"):
        eligible = [
            x for x in increases if x.get(k, {}).get("length_floor_kpc") is not None
        ]
        tight[k] = max(eligible, key=lambda x: x[k]["length_floor_kpc"], default=None)
    # Witness masses are placed at the recorded ring for interior rings
    # and at the map edge for any extra mass first needed farther out.
    deposits = []
    if local:
        last = Q(0)
        for k, (r, _, _) in enumerate(data):
            goal = prefix if r >= edge else step[k]
            if goal > last:
                deposits.append((min(r, edge), goal - last))
                last = goal
        assert all(
            lo <= sum(w for d, w in deposits if d <= r) <= hi for r, lo, hi in data
        )
    return dict(
        edge_kpc=float(edge),
        exterior_rings=len(exterior),
        unlimited_positive_anywhere=monotone,
        positive_within_map=local,
        first_positive_conflict=witness,
        deposits=[[str(d), str(w)] for d, w in deposits],
        strongest_exterior_pair=strongest,
        farthest_exterior_pair=farthest,
        tightest_hop=tight,
    )


def run():
    inputs, charts, _, masses = budget.load_inputs()
    native = budget.local.lab.HERE / "poisson_response.json"
    minimum = min(
        b["degree_range"][0] for b in json.loads(native.read_bytes())["boards"]
    )
    cosmic = HERE / "mixture.json"
    cap = Q(str(json.loads(cosmic.read_bytes())["highest_galaxy_cap_at_95"]))
    panels = {
        n: panel(r, charts[n][0][-1], masses[n], cap, minimum)
        for n, r in inputs.items()
    }
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                native,
                cosmic,
                HERE / "EXTERIOR_MARKS.md",
                HERE / "source_profile.json",
                budget.local.profiles.maps.HERE / "leroy_table7.dat",
                budget.local.profiles.maps.HERE / "leroy_sample.dat",
                budget.local.profiles.galaxies.sparc.DATA / "table1.dat",
            )
        },
        cosmic_ceiling=str(cap),
        minimum_live_degree=minimum,
        galaxies=panels,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = run()
    path = HERE / "exterior_marks.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        assert json.loads(path.read_bytes()) == result
    for k, v in result["galaxies"].items():
        w = v["strongest_exterior_pair"] or {}
        t = v["tightest_hop"]["occupation"] or {}
        print(
            k,
            v["positive_within_map"],
            v["unlimited_positive_anywhere"],
            v["edge_kpc"],
            w.get("minimum_net_ratio"),
            t.get("occupation"),
        )
