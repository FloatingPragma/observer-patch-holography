"""Run two ordinary histories whose calibrated source scores agree at every receiver."""

import argparse
from hashlib import sha256
import json

import numpy as np

from . import distributed_source as dist

HERE, lab, reg, iv = dist.HERE, dist.lab, dist.reg, dist.iv
region, mat = dist.regions, dist.service.mat
SEED = 28112121


def run():
    iv.dps = 35
    calibration = dist.read("regional_source_calibration")["calibrations"][
        "energy+connected+marginal"
    ]
    data = dist.read("distributed_source")
    held, record = next(iter(calibration["ordinary"].items()))
    n, den = record["common_physical_dose"]
    scale = n / den
    coefficients = {
        held: 1.0,
        **{k: -(a + b) / 2 for k, (a, b) in record["coefficients"].items()},
    }
    p = {
        k: {
            arm: 36 * scale * max(0, sign * v)
            for arm, sign in (("plus", 1), ("minus", -1))
        }
        for k, v in coefficients.items()
    }
    assert all(0 <= v <= 1 for z in p.values() for v in z.values())
    rng = np.random.default_rng(SEED)
    ports, bins = lab.six.survey.return_bins()
    _, src, *_ = region.c.trials()
    x = src["reference"]
    states = {k: np.zeros((7, *x.shape), complex) for k in ("plus", "minus", "dense")}
    for home in ports:
        for sign in (-1, 1):
            menus = {
                k: region.density.history(
                    src[k], (home,), rng, signs=(sign,), read_ports=ports
                )[2]
                for k in ("ordinary", "reference")
            }
            for t in region.density.WINDOWS:
                k = f"{t}/{home}/{sign}"
                low = menus["reference"][t]["dense"]
                delta = menus["ordinary"][t]["dense"] - low
                states["dense"] += low / 36
                for arm in ("plus", "minus"):
                    states[arm] += (low + p.get(k, {}).get(arm, 0.0) * delta) / 36
    arms = ("plus", "minus")
    ref = lab.pairs(x, rng, False)
    tape = {k: region.joint_tape(states[k], rng) for k in arms}
    prices = region.costs(reg.src.observed(ref, lab.SHOTS, False))
    locations = [str(z) for z in ports] + ["complement"]
    source = {
        n: {
            key: lab.six.bounds(v)
            for key, v in region.local_contrasts(
                prices, *(region.at_port(tape[k], i) for k in arms)
            ).items()
        }
        for i, n in enumerate(locations)
    }
    aggregate = {
        k: reg.src.observed(np.asarray(tape[k]).sum(axis=2), lab.SHOTS, False)
        for k in arms
    }
    dose = region.density.DOSE
    pauli = {}
    for parity in range(2):
        for i in range(4):
            for j in range(4):
                op = region.c.dose.heads.PRODUCTS[i, j]
                diff = np.trace(
                    op @ (aggregate["plus"][0][parity] - aggregate["minus"][0][parity])
                ).real
                e = sum(aggregate[k][1][parity] for k in arms) + 1e-13
                pauli[f"{parity}/{i}/{j}"] = lab.six.bounds(
                    dose * reg.band(diff - e, diff + e)
                )
    radar = {
        k: [
            lab.counts(np.r_[v, 1 - v.sum()], rng)
            for v in region.c.station.pair(bins, z)
        ]
        for k, z in {"reference": x, **{k: a.sum(0) for k, a in states.items()}}.items()
    }
    means, _ = reg.timing_slopes(radar, lab.SHOTS, (("plus", "minus", 1),))
    start = [
        (1 - dose) * a + dose * b
        for a, b in zip(means["reference"], means["dense"], strict=True)
    ]
    negative = [
        (1 - dose) * a + dose * b
        for a, b in zip(means["reference"], means["minus"], strict=True)
    ]
    diff = [dose * (a - b) for a, b in zip(means["plus"], means["minus"], strict=True)]
    effects = np.einsum("rtbij,t->rbij", bins, reg.SCORE[:-1] - 3)
    head = [
        dose * mat.forecast(effect, aggregate["plus"], aggregate["minus"])
        for effect in effects
    ]
    shape = {}
    columns = data["result"]["columns"]
    for edge, row in zip(
        region.c.EDGES[1:],
        region.c.clock_map()[1:] - region.c.clock_map()[0],
        strict=True,
    ):
        key = f"ratio/{edge[0]}/{edge[1]}"

        def project(v):
            return sum(float(c) * val for c, val in zip(row, v, strict=True) if c)

        exact = scale * iv.mpf(record["observed_minus_prediction"][key])
        offset = sum(
            (iv.mpf(w) + coefficients[k]) * iv.mpf(columns[k][key])
            for k, w in record["coefficients"].items()
        )
        tangent = project([dr / r for dr, r in zip(diff, start, strict=True)])
        forecast = project([dr / r for dr, r in zip(head, start, strict=True)])
        finite = project(
            [iv.log1p(dr / r) for dr, r in zip(diff, negative, strict=True)]
        )
        # Desired programmed direction minus exact direction is sum((w-center)*column).
        shape[key] = {
            "exact_equal_source_forecast": lab.six.bounds(exact),
            "center_forecast_with_setting_error": lab.six.bounds(
                exact + scale * offset
            ),
            "center_coincidence_forecast": lab.six.bounds(forecast),
            "center_clock_tangent": lab.six.bounds(tangent),
            "center_finite_log_ratio": lab.six.bounds(finite),
        }
    return {
        "seed": SEED,
        "starts": lab.SHOTS,
        "held_history": held,
        "dose": [n, den],
        "ordinary_probabilities": p,
        "reference_heads": ref,
        "address_counts": tape,
        "clock_counts": radar,
        "regional_source": source,
        "summed_pauli_contrasts": pauli,
        "shape": shape,
        "inputs": {
            n: sha256((HERE / n).read_bytes()).hexdigest()
            for n in (
                "MATCHED_REGIONAL_FLIGHT.md",
                "regional_source_calibration.json",
                "distributed_source.json",
            )
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    r = run()
    path = HERE / "matched_regional_flight.json"
    if args.write:
        path.write_text(json.dumps(r, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), r)
    print(r["held_history"], r["dose"], r["shape"]["ratio/1/2"])
