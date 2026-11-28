"""Read the impulse and instantaneous force of both probes after omitting an interrogation."""

import argparse
import json
from hashlib import sha256

import numpy as np

from . import finite_dose_material as beam
from . import observation_density as density

c, reg, lab, iv = density.c, density.reg, density.lab, density.iv
HERE = density.HERE
METERS = ("marked_force", "free_rate", "shift")
POWERS = (68, 72)
SEED = 28112109


def forecast(effect, high, low):
    a, ea = high
    b, eb = low
    co = np.einsum("pqij,bji->bpq", c.dose.heads.PRODUCTS, effect).real / 4
    radius = abs(co).sum(axis=(1, 2)) @ (ea + eb) + 1e-12
    mu = c.station.pair(effect, a - b)
    return reg.band(mu - radius, mu + radius)


def acquire(rng, n, source, effect, occupancy, scale):
    p = c.station.pair(occupancy, source) / scale
    mu = c.station.pair(effect, source) / scale
    return beam.port_read(rng, n, p, mu, scale, np.log(2 * 4000 / 0.01))


def run(tagged=False, repeat=False):
    iv.dps = 35
    seed = SEED + int(tagged) + int(repeat)
    rng = np.random.default_rng(seed)
    prior = (
        json.loads((HERE / "density_material.json").read_bytes()) if tagged else None
    )
    _, sources, *_ = c.trials()
    x = sources["reference"]
    ports, _ = lab.six.survey.return_bins()
    arms, _, active = density.history(x, ports, rng)
    _, _, ordinary = density.history(sources["ordinary"], ports, rng)
    effects, scales, hits = beam.flight(METERS, include_edges=False)
    energy = np.stack([np.zeros((4, 4)), np.eye(4)])
    windows = {}
    windows_to_read = (2,) if repeat else density.WINDOWS
    powers = (76,) if repeat else (72,) if tagged else POWERS
    for tau in windows_to_read:
        states = dict(arms[tau])
        states["ordinary"] = (1 - density.DOSE) * x + density.DOSE * ordinary[tau][
            "dense"
        ].sum(0)
        if tagged:
            # Forecast from the original, independently diluted head ensembles.
            heads = prior["windows"][str(tau)]["source_heads"]
            states = {k: z.sum(0) for k, z in active[tau].items()}
            states["ordinary"] = ordinary[tau]["dense"].sum(0)
        else:
            heads = {k: lab.pairs(z, rng, False) for k, z in states.items()}
        src = {k: reg.src.observed(t, lab.SHOTS, False) for k, t in heads.items()}
        ordinary_cost = forecast(energy, src["ordinary"], src["dense"])
        predicted = {
            (key, arm): forecast(np.array(effect), src[arm], src["dense"])
            for key, effect in effects.items()
            for arm in ("sparse", "ordinary", "erased")
        }
        exposures = {}
        for power in powers:
            n = 2**power
            tape, readings = {}, {}
            # Freeze both probes' ordinary calibrations before any held read.
            for arm in ("dense", "ordinary", "sparse", "erased"):
                for key, effect in effects.items():
                    tape[arm, key], readings[arm, key] = acquire(
                        rng, n, states[arm], effect, hits[key], scales[key]
                    )
            meters = {}
            for meter in METERS:
                result, changes = {}, {}
                for sign in (-1, 1):
                    key = (sign, meter)
                    diff = {
                        a: (readings[a, key] - readings["dense", key])
                        * (density.DOSE if tagged else 1)
                        for a in ("sparse", "ordinary", "erased")
                    }
                    denominator = diff["ordinary"]
                    detail = {
                        "analyzer_scale": scales[key],
                        "counts": {a: tape[a, key] for a in states},
                        "contrasts": {a: lab.six.bounds(z) for a, z in diff.items()},
                        "source_forecasts": {
                            a: lab.six.bounds(predicted[key, a]) for a in diff
                        },
                    }
                    if not denominator.a <= 0 <= denominator.b:
                        detail["ordinary_energy_equivalent"] = lab.six.bounds(
                            ordinary_cost * diff["sparse"] / denominator
                        )
                    result[str(sign)] = detail
                    changes[sign] = diff
                a, b = changes[-1], changes[1]
                determinant = a["sparse"] * b["ordinary"] - b["sparse"] * a["ordinary"]
                meters[meter] = {
                    "probes": result,
                    "common_response_determinant": lab.six.bounds(determinant),
                }
            exposures[str(power)] = meters
        windows[str(tau)] = {
            "source_heads": heads,
            "ordinary_energy": lab.six.bounds(ordinary_cost),
            "sparse_energy": lab.six.bounds(
                forecast(energy, src["sparse"], src["dense"])
            ),
            "exposures": exposures,
        }
    return {
        "seed": seed,
        "windows": windows,
        "inputs": {
            n: sha256((HERE / n).read_bytes()).hexdigest()
            for n in (
                "DENSITY_MATERIAL_PROTOCOL.md",
                "observation_density.json",
                "finite_pairing_dose.json",
                "ordinary_stationarity.json",
                *(
                    ("DENSITY_MATERIAL_ACTIVE.md", "density_material.json")
                    if tagged
                    else ()
                ),
                *(
                    ("DENSITY_MATERIAL_IMPULSE.md", "density_material_active.json")
                    if repeat
                    else ()
                ),
            )
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    cli.add_argument(
        "--active", action="store_true", help="Acquire the tagged active flights only."
    )
    cli.add_argument(
        "--repeat",
        action="store_true",
        help="Repeat the contracting window at higher tagged exposure.",
    )
    args = cli.parse_args()
    result = run(args.active or args.repeat, args.repeat)
    filename = "density_material.json"
    if args.repeat:
        filename = "density_material_impulse.json"
    elif args.active:
        filename = "density_material_active.json"
    path = HERE / filename
    if args.write:
        path.write_text(json.dumps(result, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), result)
    for tau, window in result["windows"].items():
        print(
            tau,
            {
                name: [
                    r.get("ordinary_energy_equivalent")
                    for r in meter["probes"].values()
                ]
                for name, meter in window["exposures"][
                    "76" if args.repeat else "72"
                ].items()
            },
        )
