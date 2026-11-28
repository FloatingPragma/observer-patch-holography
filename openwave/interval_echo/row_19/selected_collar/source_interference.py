"""Predict held finite geometry from the click/miss interference measured on the source."""

import argparse
import json
from hashlib import sha256
from os.path import relpath
from pathlib import Path

import numpy as np

from . import observation_density as density

pos, c, lab, reg, iv = density.pos, density.c, density.lab, density.reg, density.iv
HERE = pos.HERE
PHASES = {"plus": 1, "minus": -1, "left": 1j, "right": -1j}
HELD = {
    "omitted": 1,
    "partial_a": 0.3 + 0.4j,
    "partial_b": -0.4 + 0.3j,
    "orthogonal": 0,
}
SEED = 28112105


def trace_bound(pair):
    x, e = pair
    return sum(np.linalg.norm(m, ord="nuc") for m in x) + 4 * sum(e) + 1e-12


def difference(a, b):
    return a[0] - b[0], a[1] + b[1]


def read(price, observation, centered=True):
    m, e = observation
    # The contrast has trace zero. Centering the physical clock effect
    # avoids charging uncertainty on its common constant.
    scalar = np.trace(price, axis1=1, axis2=2).real.sum() / 8 if centered else 0
    K = price - scalar * np.eye(4)
    value = c.station.pair(K, m)
    coeff = (
        abs(np.einsum("ijkl,blk->bij", c.dose.heads.PRODUCTS, K)).sum(axis=(1, 2)) / 4
    )
    err = coeff @ e + 2e-12
    return reg.band(value - err, value + err)


def volume_from_source(baseline, delta, ranges):
    start = [read(R, baseline, False) for R in ranges]
    dr = [read(R, delta) for R in ranges]
    path = [x + reg.band(0, 1) * d for x, d in zip(start, dr, strict=True)]
    g = pos.geometry_jet(path)[1].g
    mid = np.array([float(v.mid) for v in g])
    rem = sum(abs(v - a) * abs(d) for a, v, d in zip(mid, g, dr, strict=True))
    val = read(np.tensordot(mid, ranges, axes=1), delta)
    return val + reg.band(-1, 1) * (rem + 1e-12 * (1 + abs(mid).sum()))


def run():
    iv.dps = 35
    rng = np.random.default_rng(SEED)
    source, states, *_ = c.trials()
    x = states["reference"]
    ports, bins = lab.six.survey.return_bins()
    ids = [i for i, z in enumerate(lab.six.survey.ROUTES) if max(z) < 4]
    bins = bins[ids]
    ranges = 3 * np.eye(4) + np.einsum("rtbij,t->rbij", bins, reg.SCORE[:-1] - 3)
    pointers = {
        **PHASES,
        "partial_a": HELD["partial_a"],
        "partial_b": HELD["partial_b"],
    }
    menus, _, _ = density.history(x, ports, rng, pointers)
    source_ref = lab.pairs(x, rng, False)
    ref = reg.src.observed(source_ref, lab.SHOTS, False)
    base = HERE / "link_source_contraction.json"
    key = Path(c.dose.enclosure.__file__).with_name("certificate.json")
    cert = json.loads(key.read_bytes())
    eta = iv.mpf(cert["gain_interval"])
    lam = iv.mpf(json.loads(base.read_bytes())["ordinary_volume_multiplier"])
    hs = cert["state_HS_error_upper"] + 1e-12
    floor = np.linalg.eigvalsh(ref[0]).min() - 2 * max(ref[1]) - hs - 1e-12
    assert floor > 0
    L = cert["global_survey"]["trace_norm_curvature_limit"]
    reference_error = hs / floor + L * np.sqrt(8) * hs + 2e-12
    windows = {}
    for tau, menu in menus.items():
        heads = {k: lab.pairs(menu[k], rng, False) for k in PHASES}
        obs = {k: reg.src.observed(t, lab.SHOTS, False) for k, t in heads.items()}
        quadratures = {
            k: [
                density.price_contrast(ref, obs[a], obs[b], k) / 2
                for a, b in (("plus", "minus"), ("left", "right"))
            ]
            for k in ("full", "connected", "marginals")
        }
        dense = tuple(
            (a + b) / 2 for a, b in zip(obs["plus"], obs["minus"], strict=True)
        )
        f, ef = difference(obs["plus"], obs["minus"])
        j, ej = difference(obs["left"], obs["right"])
        dense_offset = trace_bound(difference(dense, ref))
        predictions = {}
        for k, v in HELD.items():
            dx = (v.real * f + v.imag * j) / 2
            ex = (abs(v.real) * ef + abs(v.imag) * ej) / 2
            pred = (dense[0] + dx, dense[1] + ex)
            length = trace_bound((dx, ex))
            if v == 0:
                length = 0.0
            radius = max(dense_offset, trace_bound(difference(pred, ref)))
            remainder = (reference_error + L * radius) * length + 1e-12
            full = v.real * quadratures["full"][0] + v.imag * quadratures["full"][1]
            curved = full + reg.band(-remainder, remainder)
            volume = volume_from_source(dense, (dx, ex), ranges)
            area = (lam * volume - curved) / eta
            predictions[k] = {
                "pointer": [v.real, v.imag],
                "chord_trace_norm_upper": float(length),
                "max_offset_from_reference": float(radius),
                "geometric_remainder_nats": float(remainder),
                "source_full_price": lab.six.bounds(full),
                "adjusted_contraction_nats": lab.six.bounds(curved),
                "source_volume": lab.six.bounds(volume),
                "source_area": lab.six.bounds(area),
            }
        # No held clock has been acquired at the point the forecasts are fixed.
        check = {
            "dense": menu["dense"],
            "omitted": menu["sparse"],
            "partial_a": menu["partial_a"],
            "partial_b": menu["partial_b"],
            "orthogonal": menu["erased"],
        }
        radar = {
            k: [lab.counts(np.r_[p, 1 - p.sum()], rng) for p in c.station.pair(bins, z)]
            for k, z in check.items()
        }
        means, _ = reg.timing_slopes(radar, lab.SHOTS, (("omitted", "dense", 1),))
        rb = np.array(means["dense"], dtype=object)
        for k, r in predictions.items():
            ra = np.array(means[k], dtype=object)
            d = ra - rb
            area, volume = pos.geometry_jet(
                rb + np.array([reg.band(0, 1) * v for v in d])
            )
            obs_area = area.g @ d
            obs_adj = (volume.g * lam - area.g * eta) @ d
            r["radar_area"] = lab.six.bounds(obs_area)
            r["radar_adjusted_contraction_nats"] = lab.six.bounds(obs_adj)
        windows[str(tau)] = {
            "head_quadrature_counts": heads,
            "measured_price_quadratures": {
                k: list(map(lab.six.bounds, v)) for k, v in quadratures.items()
            },
            "predictions": predictions,
            "held_clock_histograms": radar,
        }
    return {
        "seed": SEED,
        "starts_per_setting": lab.SHOTS,
        "reference_head_counts": source_ref,
        "reference_gradient_error_upper": reference_error,
        "mixed_geometry_Hessian_upper": L,
        "windows": windows,
        "inputs": {
            relpath(p, HERE): sha256(p.read_bytes()).hexdigest()
            for p in (
                base,
                key,
                HERE / "SOURCE_INTERFERENCE_PROTOCOL.md",
                HERE / "finite_pairing_dose.json",
            )
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    data = run()
    path = HERE / "source_interference.json"
    if cli.parse_args().write:
        path.write_text(json.dumps(data, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), data)
    for tau, window in data["windows"].items():
        print(
            tau,
            {
                k: [r["source_area"], r["radar_area"]]
                for k, r in window["predictions"].items()
            },
        )
