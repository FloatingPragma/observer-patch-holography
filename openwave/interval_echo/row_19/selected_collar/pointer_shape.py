"""Test reconstructed spatial shape independently of sender clock gains and common scale."""

import argparse
from hashlib import sha256
import json

import numpy as np

from . import charge_selected_flight as chosen
from . import pointer_mass_equivalence as cone
from . import regional_density as regions

mat, lab, reg, iv = chosen.material, chosen.lab, chosen.reg, chosen.iv
HERE = mat.HERE


def run():
    iv.dps = 35

    def read(n):
        return json.loads((HERE / (n + ".json")).read_bytes())

    data = read("pointer_mass_equivalence")
    six = read("six_meter_pointer")
    original = read("regional_density")["windows"]
    control = read("ordinary_service_distribution")["windows"]
    selected = read("charge_selected_flight")["windows"]
    times = list(selected)
    routes = lab.six.survey.ROUTES
    clocks = [f"clock/{a}/{b}" for a, b in routes]
    graph = mat.c.clock_map()
    contrasts = graph[1:] - graph[0]
    labels = [f"length_ratio/{a}/{b}" for a, b in mat.c.EDGES[1:]]
    changes = np.array([[int(i == source) for source in range(6)] for i, j in routes])
    np.testing.assert_array_equal(contrasts @ changes, np.zeros((14, 6)))
    dense = {}
    for t, r in original.items():
        timings, _ = reg.timing_slopes(
            r["radar_counts"], lab.SHOTS, (("sparse", "dense", 1),)
        )
        dense[t] = timings["dense"]

    def band(v):
        return iv.mpf(v)

    bounds = []
    multipliers = []
    for j, n in enumerate(clocks):
        arms = []
        for t in times:
            u, v = map(band, selected[t]["meter_quadratures"][n])
            # Both endpoints are inside this convex hull for any time lottery.
            span = iv.sqrt(u**2 + v**2)
            arms.extend(
                [
                    dense[t][j] + reg.band(-1, 1) * span,
                    dense[t][j] + reg.band(0, 1) * band(control[t]["held_ordinary"][n]),
                ]
            )
        interval = iv.mpf([min(a.a for a in arms), max(a.b for a in arms)])
        assert interval.a > 0
        bounds.append(lab.six.bounds(interval))
        multipliers.append(1 / interval)

    def project(panel):
        scaled = [v * band(panel[n]) for v, n in zip(multipliers, clocks, strict=True)]
        out = {k: panel[k] for k in data["groups"]["material"]}
        for n, row in zip(labels, contrasts, strict=True):
            out[n] = lab.six.bounds(
                sum(float(v) * k for v, k in zip(row, scaled, strict=True) if v)
            )
        return out

    q = {t: [project(p) for p in pair] for t, pair in data["quadratures"].items()}
    ordinary = {t: project(z) for t, z in data["ordinary"].items()}
    groups = {
        "shape": labels,
        "shape_and_material": labels + data["groups"]["material"],
    }
    cases = {t: cone.comparison(q, ordinary, [t], groups) for t in times}
    cases["mixture"] = cone.comparison(q, ordinary, times, groups)
    totals = {
        t: {
            g: {
                stat: sum(v["status"] == stat for v in tests.values())
                for stat in ("excluded", "compatible_relaxation", "unresolved")
            }
            for g, tests in rows.items()
        }
        for t, rows in cases.items()
    }
    beta = iv.mpf(six["common_ratio"][0]) / six["common_ratio"][1]
    # Both endpoint ranges share the independently acquired dense baseline.
    start = np.array(
        [
            sum(
                dense[t][j] + beta * band(control[t]["held_ordinary"][n]) for t in times
            )
            / 3
            for j, n in enumerate(clocks)
        ],
        dtype=object,
    )
    delta = np.array(
        [band(six["clock_minus_equal_ordinary"][n]) for n in clocks], dtype=object
    )
    balls = regions.ball_change(start, start + delta, delta=delta)
    rows = [
        {
            "triangle": list(tri),
            "vertex": int(vertex),
            "area_sparse_minus_equal_ordinary": lab.six.bounds(v),
        }
        for (tri, vertex), v in balls.items()
    ]
    return {
        "reference_edge": list(mat.c.EDGES[0]),
        "edges": [list(e) for e in mat.c.EDGES],
        "clock_ranges": bounds,
        "ratio_matrix": contrasts.tolist(),
        "quadratures": q,
        "ordinary": ordinary,
        "groups": groups,
        "counts": totals,
        "certificates": cases,
        "matched_proper_volume_balls": rows,
        "inputs": {
            n: sha256((HERE / n).read_bytes()).hexdigest()
            for n in (
                "POINTER_SHAPE_PROTOCOL.md",
                "pointer_mass_equivalence.json",
                "six_meter_pointer.json",
                "regional_density.json",
                "ordinary_service_distribution.json",
                "charge_selected_flight.json",
            )
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    r = run()
    path = HERE / "pointer_shape.json"
    if args.write:
        path.write_text(json.dumps(r, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), r)
    print(r["counts"])
    print(
        "resolved fixed-volume differences",
        sum(
            a * b > 0
            for a, b in (
                row["area_sparse_minus_equal_ordinary"]
                for row in r["matched_proper_volume_balls"]
            )
        ),
    )
