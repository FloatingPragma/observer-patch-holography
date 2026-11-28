"""Measure nonlinear geometry from independently tagged clock means."""

import argparse
from hashlib import sha256
import json

import numpy as np

from . import six_meter_pointer as source
from . import regional_density as regions

HERE, mat, reg, lab, iv = source.HERE, source.mat, source.reg, source.lab, source.iv


def run():
    iv.dps = 35
    p = HERE / "six_meter_pointer.json"
    prior = json.loads(p.read_bytes())
    beta = iv.mpf(prior["common_ratio"][0]) / prior["common_ratio"][1]
    times = prior["times"]
    settings = {t: complex(*prior["central_settings"][t]) for t in times}
    _, prep, *_ = mat.c.trials()
    rng = np.random.default_rng(28112119)
    ports, bins = lab.six.survey.return_bins()
    _, _, heads = mat.density.history(
        prep["reference"], ports, rng, {"held_" + t: g for t, g in settings.items()}
    )
    _, _, ordinary = mat.density.history(prep["ordinary"], ports, rng)
    states = {
        "reference": prep["reference"],
        "dense": sum(heads[int(t)]["dense"].sum(0) for t in times) / 3,
        "held": sum(heads[int(t)]["held_" + t].sum(0) for t in times) / 3,
        "ordinary": sum(ordinary[int(t)]["dense"].sum(0) for t in times) / 3,
    }
    radar = {
        k: [lab.counts(np.r_[p, 1 - p.sum()], rng) for p in mat.c.station.pair(bins, y)]
        for k, y in states.items()
    }
    means, _ = reg.timing_slopes(radar, lab.SHOTS, (("held", "dense", 1),))
    r = {k: np.array(v, dtype=object) for k, v in means.items()}
    dose = mat.density.DOSE
    start = (1 - dose) * r["reference"] + dose * r["dense"]
    held = dose * (r["held"] - r["dense"])
    paid = np.array([beta * v for v in dose * (r["ordinary"] - r["dense"])])
    delta = held - paid
    anchor = [i for i, route in enumerate(lab.six.survey.ROUTES) if max(route) < 4]

    def geometry(a, d):
        return source.selected.geometry_change(a, d, anchor)

    graph = mat.c.clock_map()
    log_difference = [iv.log1p(v / u) for v, u in zip(delta, start + paid, strict=True)]
    shape = {
        str(edge): lab.six.bounds(
            sum(float(c) * z for c, z in zip(row, log_difference, strict=True) if c)
        )
        for edge, row in zip(mat.c.EDGES[1:], graph[1:] - graph[0], strict=True)
    }
    balls = regions.ball_change(start + paid, start + held, delta=delta)
    return {
        "seed": 28112119,
        "starts": lab.SHOTS,
        "clock_counts": radar,
        "held_geometry": geometry(start, held),
        "ordinary_geometry": geometry(start, paid),
        "held_minus_ordinary_geometry": geometry(start + paid, delta),
        "log_edge_ratio_differences": shape,
        "proper_volume_area_differences": [
            {"triangle": list(tri), "vertex": i, "difference": lab.six.bounds(z)}
            for (tri, i), z in balls.items()
        ],
        "inputs": {
            f.name: sha256(f.read_bytes()).hexdigest()
            for f in (p, HERE / "SIX_METER_CLOCKS.md")
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    r = run()
    path = HERE / "six_meter_clocks.json"
    if args.write:
        path.write_text(json.dumps(r, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), r)
    print(r["held_geometry"], r["ordinary_geometry"], r["held_minus_ordinary_geometry"])
