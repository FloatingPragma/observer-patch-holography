"""Match the complete local transfer law and volume with an opposite modular source."""

import argparse
import json
from hashlib import sha256

import numpy as np

from . import finite_dose_material as beam
from . import position_source as pos

c, lab, iv, reg = pos.c, pos.lab, pos.iv, pos.reg
HERE = pos.HERE


def preparations(x, held, signed, hit, ranges):
    basis = c.dose.heads.measure.joint_basis()

    def coeff(E):
        return c.station.pair(basis, E)

    r = c.station.pair(ranges, held)
    v = pos.geometry_jet([iv.mpf(float(z)) for z in r])[1]
    vol = np.tensordot([float(w.mid) for w in v.g], ranges, axes=1)
    K = -c.station.log_blocks(x)
    raw = np.array([coeff(e) for e in [signed, hit, vol, K]])
    normalized = raw / np.linalg.norm(raw, axis=1)[:, None]
    dual = np.linalg.pinv(normalized)
    directions = (dual / np.linalg.norm(dual, axis=0)).T
    directions = np.tensordot(directions, basis, axes=1)
    for exponent in range(6, 18):
        step = 2.0**-exponent
        center = held - step * directions[3]
        sources = {}
        for i, side in enumerate((-1, 1)):
            for j, (a, b) in enumerate(((1, 0), (0, 1), (-1, -1))):
                sources[f"{i}{j}"] = center + step * (
                    side * directions[2] / 4
                    + (a * directions[0] + b * directions[1]) / 8
                )
        if all(np.linalg.eigvalsh(y).min() > 0 for y in sources.values()):
            try:
                recipes = {k: c.station.preparation(y) for k, y in sources.items()}
            except ArithmeticError:
                continue
            if (
                min(min(v["min_kept"], v["min_dephased"]) for v in recipes.values()) > 0
                and max(v["max_error"] for v in recipes.values()) < 1e-12
            ):
                return sources, recipes, step
    raise ValueError("No open reachable fabric for the registered competitors")


def det(a, b):
    return a[0] * b[1] - a[1] * b[0]


def containment(three, target):
    a, b, z = three
    denominator = det(b - a, z - a)
    if denominator.a <= 0 <= denominator.b:
        return None
    return [
        lab.six.bounds(det(t - target, v - target) / denominator)
        for t, v in ((b, z), (z, a), (a, b))
    ]


def hull_geometry(means, ref, labels, component):
    start = np.asarray(means[ref], dtype=object)
    delta = [np.asarray(means[k], dtype=object) - start for k in labels]
    box = np.array(
        [
            iv.mpf(
                [
                    min([v.a] + [means[k][i].a for k in labels]),
                    max([v.b] + [means[k][i].b for k in labels]),
                ]
            )
            for i, v in enumerate(start)
        ]
    )
    grad = pos.geometry_jet(box)[component].g
    # On every segment from start, the gradient stays in the same interval
    # box. Its pairing with a convex combination is bounded at the vertices.
    candidates = [grad @ d for d in delta]
    return lab.six.bounds(
        iv.mpf([min(z.a for z in candidates), max(z.b for z in candidates)])
    )


def run():
    iv.dps = 35
    _, original, *_ = c.trials()
    x, held = (original[k] for k in ("reference", "held"))
    edge = (1, 0, 3)
    effects, scales, occupancy = beam.flight()
    scale = scales[edge]
    hit = np.asarray(occupancy[edge]) / scale
    signed = np.array(effects[edge])
    _, all_bins = lab.six.survey.return_bins()
    ids = [i for i, r in enumerate(lab.six.survey.ROUTES) if max(r) < 4]
    bins = all_bins[ids]
    ranges = 3 * np.eye(4) + np.einsum("rtbij,t->rbij", bins, reg.SCORE[:-1] - 3)
    alt, recipes, step = preparations(x, held, signed, hit, ranges)
    keylist = list(alt)
    states = {"reference": x, "held": held, **alt}
    seed = 28112102
    rng = np.random.default_rng(seed)
    n = 2**68
    readings, raw, time, head = {}, {}, {}, {}
    for name, y in states.items():
        hits, transfer = beam.port_read(
            rng,
            n,
            c.station.pair(hit, y),
            c.station.pair(signed, y) / scale,
            scale,
            np.log(40000),
        )
        prob = (hits[0] + hits[1]) / n
        err = reg.src.probability_error(prob, n) + 2e-12
        raw[name] = hits
        readings[name] = np.array(
            [transfer, reg.band(prob - err, prob + err)], dtype=object
        )
        time[name] = [
            lab.counts(np.append(p, 1 - p.sum()), rng) for p in c.station.pair(bins, y)
        ]
        head[name] = lab.pairs(y, rng, False)
    means, _ = reg.timing_slopes(time, lab.SHOTS, (("held", "reference", 1),))
    # All comparison mixtures have the same marginals and parity as x.
    # Thus the source charge equals the frozen full modular score.
    estimated, ex = reg.src.observed(head["reference"], lab.SHOTS, False)
    K = -c.station.log_blocks(estimated)
    floors = np.linalg.eigvalsh(estimated).min(axis=1) - 2 * ex - 1e-12
    assert min(floors) > 0
    charge = {}
    for name in ["held", *keylist]:
        y, ey = reg.src.observed(head[name], lab.SHOTS, False)
        change = y - estimated
        center = c.station.pair(K, change)
        error = 1e-12
        for a in range(2):
            pauli = abs(np.einsum("ijkl,lk->ij", c.dose.heads.PRODUCTS, K[a])).sum() / 4
            error += (
                pauli * (ey[a] + ex[a])
                + (np.linalg.norm(change[a], ord="nuc") + 4 * (ex[a] + ey[a]))
                * 2
                * ex[a]
                / floors[a]
            )
        charge[name] = lab.six.bounds(reg.band(center - error, center + error))
    tri = []
    for group in (keylist[:3], keylist[3:]):
        tri.append(
            {
                "sources": group,
                "barycentric_weights": containment(
                    [readings[k] for k in group], readings["held"]
                ),
                "volume_minus_held": hull_geometry(means, "held", group, 1),
            }
        )
    return {
        "seed": seed,
        "link_shots": n,
        "radar_and_head_shots": lab.SHOTS,
        "edge": list(edge),
        "step": step,
        "rail_lotteries": recipes,
        "link_counts": raw,
        "head_counts": head,
        "time_counts": time,
        "modular_change": charge,
        "triangles": tri,
        "all_alternatives_area_minus_reference": hull_geometry(
            means, "reference", keylist, 0
        ),
        "held_area_minus_reference": hull_geometry(means, "reference", ["held"], 0),
        "inputs": {
            z: sha256((HERE / z).read_bytes()).hexdigest()
            for z in (
                "LINK_SOURCE_ALTERNATIVE.md",
                "finite_pairing_dose.json",
                "link_source_contraction.json",
            )
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    result = run()
    path = HERE / "link_source_alternative.json"
    if cli.parse_args().write:
        path.write_text(json.dumps(result, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), result)
    print(
        result["triangles"],
        result["all_alternatives_area_minus_reference"],
        result["modular_change"],
    )
