"""Compare a measured spatial distribution of modular source with fixed-volume ball areas."""

import argparse
import json
from fractions import Fraction
from hashlib import sha256
from itertools import combinations

import numpy as np

from . import erasure_balls as balls
from . import observation_density as density

pos, c, lab, iv, reg = density.pos, density.c, density.lab, density.iv, density.reg
HERE = pos.HERE
SEED = 28112104


def joint_tape(local, rng):
    rows = []
    for a in c.dose.heads.AXES[1:]:
        row = []
        for b in c.dose.heads.AXES[1:]:
            E = [
                np.kron((np.eye(2) + s * a) / 2, (np.eye(2) + t * b) / 2)
                for s in (1, -1)
                for t in (1, -1)
            ]
            p = np.einsum("kij,sbji->sbk", E, local).real
            row.append(np.asarray(lab.counts(p.ravel(), rng)).reshape(7, 8).tolist())
        rows.append(row)
    return rows


def at_port(tape, site):
    selected = np.asarray(tape)[:, :, site]
    missed = lab.SHOTS - selected.sum(axis=-1, keepdims=True)
    return reg.src.observed(
        np.concatenate([selected, missed], axis=-1), lab.SHOTS, True
    )


def costs(reference):
    x, errors = reference
    rad = max(2 * errors)
    left, right = c.dose.heads.measure.marginals(x)

    def logerr(m, e):
        floor = np.min(m) - e - 1e-12
        assert floor > 0
        return e / floor

    full = -c.station.log_blocks(x)
    m = -c.station.log_blocks(c.dose.heads.product_reference(x))
    em = (
        logerr(np.linalg.eigvalsh(left), 2 * rad)
        + logerr(np.linalg.eigvalsh(right), 2 * rad)
        + logerr(np.trace(x, axis1=1, axis2=2).real, 2 * rad)
    )
    ef = logerr(np.linalg.eigvalsh(x), rad)
    energy = np.stack([np.zeros((4, 4)), np.eye(4)])
    prices = {}
    for label, K, e in (
        ("connected", full - m, ef + em),
        ("marginal", m, em),
        ("energy", energy, 0.0),
    ):
        price = K - c.station.pair(K, x) * np.eye(4)
        error = 2 * e + 4 * errors.sum() * np.max(abs(np.linalg.eigvalsh(K))) + 1e-11
        prices[label] = (price, error)
    return prices


def local_contrasts(prices, a, b):
    high, eh = a
    low, el = b
    d = high - low
    size = sum(np.linalg.norm(z, ord="nuc") for z in d) + 4 * sum(eh + el)
    result = {}
    for label, (K, logerr) in prices.items():
        center = c.station.pair(K, d)
        e = logerr * size + 1e-12
        for z, w in zip(K, eh + el, strict=True):
            e += abs(np.einsum("ijkl,lk->ij", c.dose.heads.PRODUCTS, z)).sum() * w / 4
        result[label] = density.DOSE * reg.band(center - e, center + e)
    return result


def ball_change(start, end, *, delta=None):
    if delta is None:
        delta = end - start
    path = np.array(
        [
            iv.mpf([min(a.a, b.a), max(a.b, b.b)])
            for a, b in zip(start, end, strict=True)
        ]
    )
    return ball_area_rates(path, delta)


def ball_area_rates(path, direction, *, tangent=False):
    """Bound the fixed-volume area rate for balls of initial radius one half."""
    multiple = direction.ndim == 2
    square = reg.six.squared_range_jets(path, direction if multiple else direction[:, None])
    raw = np.vectorize(
        lambda z: z.v if hasattr(z, "v") else iv.mpf(0), otypes=[object]
    )(square)
    omega = {}
    alt_min = np.inf
    for tet in combinations(range(6), 4):
        for i in tet:
            rest = tuple(v for v in tet if v != i)
            omega[tet, i] = balls.solid_angle_jets(square, i, rest)
            vol = balls.geom.solid_angle(raw, i, rest)[1]
            alt_min = min(alt_min, float((3 * vol / balls.geom.face_area(raw, rest)).a))
    result = {}
    for tri in c.FILLINGS:
        cells = [t for t in c.TETRA if len(set(t) & set(tri)) == 2]
        for i in range(6):
            solid = sum(omega[t, i] for t in cells if i in t)
            ratio = iv.mpf(1) if tangent else solid.v / solid.v
            assert 0.5 * ratio.b ** (iv.mpf(1) / 3) < alt_min
            factor = iv.mpf(1) / 12 * ratio ** (iv.mpf(2) / 3)
            result[tri, i] = (
                [factor * value for value in solid.g] if multiple else factor * solid.g[0]
            )
    return result


def admissible(a, b):
    # Outward intervals on intercept and slope. Allow each combination of
    # their endpoints, which only enlarges the set consistent with zero.
    lo, hi = Fraction(0), None
    for u, v in ((a[0], b[0]), (-a[1], -b[1])):
        u, v = Fraction(u), Fraction(v)
        if v == 0 and u > 0:
            return []
        if v > 0:
            top = -u / v
            hi = top if hi is None else min(hi, top)
        if v < 0:
            lo = max(lo, -u / v)
    if hi is not None and lo > hi:
        return []
    return [
        max(0.0, float(np.nextafter(float(lo), -np.inf))),
        None if hi is None else float(np.nextafter(float(hi), np.inf)),
    ]


def overlap(sets):
    if any(not v for v in sets):
        return []
    lo = max(a for a, b in sets)
    hi = min((b for a, b in sets if b is not None), default=None)
    return [] if hi is not None and lo > hi else [lo, hi]


def run():
    iv.dps = 35
    _, inputs, *_ = c.trials()
    x, o = (inputs[k] for k in ("reference", "ordinary"))
    ports, bins = lab.six.survey.return_bins()
    rng = np.random.default_rng(SEED)
    _, _, flight = density.history(x, ports, rng)
    _, _, ordinary = density.history(o, ports, rng)
    ref = lab.pairs(x, rng, False)
    prices = costs(reg.src.observed(ref, lab.SHOTS, False))
    receipts = {}
    for tau in density.WINDOWS:
        local = {k: flight[tau][k] for k in ("sparse", "dense")}
        local["ordinary"] = ordinary[tau]["dense"]
        heads = {k: joint_tape(v, rng) for k, v in local.items()}
        sources = {
            k: (1 - density.DOSE) * x + density.DOSE * z.sum(0)
            for k, z in local.items()
        }
        times = {
            k: [lab.counts(np.r_[p, 1 - p.sum()], rng) for p in c.station.pair(bins, y)]
            for k, y in sources.items()
        }
        means, _ = reg.timing_slopes(times, lab.SHOTS, (("sparse", "dense", 1),))
        measured = {
            k: [
                local_contrasts(prices, at_port(t, i), at_port(heads["dense"], i))
                for i in range(7)
            ]
            for k, t in heads.items()
            if k != "dense"
        }
        area = {
            k: ball_change(
                np.array(means["dense"], dtype=object), np.array(means[k], dtype=object)
            )
            for k in ("sparse", "ordinary")
        }
        rows = []
        for tri in c.FILLINGS:
            for i, port in enumerate(ports):
                s, o = measured["sparse"][i], measured["ordinary"][i]
                a, h = area["sparse"][tri, i], area["ordinary"][tri, i]
                determinant = {k: lab.six.bounds(s[k] * h - o[k] * a) for k in s}
                rows.append(
                    {
                        "triangle": list(tri),
                        "port": int(port),
                        "sparse_area_change": lab.six.bounds(a),
                        "ordinary_area_change": lab.six.bounds(h),
                        "determinants": determinant,
                        "allowed_charge_weight": {
                            k: admissible(determinant[k], determinant["connected"])
                            for k in ("marginal", "energy")
                        },
                    }
                )
        receipts[str(tau)] = {
            "joint_terminal_counts": heads,
            "radar_counts": times,
            "source_by_receiver": {
                k: [{j: lab.six.bounds(v) for j, v in row.items()} for row in m]
                for k, m in measured.items()
            },
            "balls": rows,
        }
    rows = [r for t in receipts.values() for r in t["balls"]]
    return {
        "ports": list(ports),
        "receiver_complement_index": 6,
        "seed": SEED,
        "starts_per_setting": lab.SHOTS,
        "source_fraction": density.DOSE,
        "reference_counts": ref,
        "windows": receipts,
        "allowed_by_filling": {
            str(tri): {
                k: overlap(
                    [
                        r["allowed_charge_weight"][k]
                        for r in rows
                        if r["triangle"] == list(tri)
                    ]
                )
                for k in ("marginal", "energy")
            }
            for tri in c.FILLINGS
        },
        "inputs": {
            n: sha256((HERE / n).read_bytes()).hexdigest()
            for n in (
                "REGIONAL_DENSITY_PROTOCOL.md",
                "OBSERVATION_DENSITY_PROTOCOL.md",
                "finite_pairing_dose.json",
                "ordinary_stationarity.json",
            )
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    d = run()
    path = HERE / "regional_density.json"
    if cli.parse_args().write:
        path.write_text(json.dumps(d, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), d)
    print(d["allowed_by_filling"])
