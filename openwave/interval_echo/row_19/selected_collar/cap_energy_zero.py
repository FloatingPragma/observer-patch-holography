"""Test local cap remainders against arbitrary sector energy origins."""

import argparse
from hashlib import sha256
from itertools import product
import json

import numpy as np
from scipy.optimize import linprog

from . import distributed_source as dist
from .pointer_mass_equivalence import negate

HERE, reg, iv, lab, region = dist.HERE, dist.reg, dist.iv, dist.lab, dist.regions
SEED = 28112126


def acquire():
    rng = np.random.default_rng(SEED)
    ports = dist.read("distributed_source")["ports"]
    settings = {
        int(t): complex(*v)
        for t, v in dist.read("six_meter_pointer")["central_settings"].items()
    }
    _, prep, *_ = region.c.trials()
    counts = {}
    for sign in (-1, 1):
        _, _, local = region.density.history(
            prep["reference"],
            ports,
            rng,
            {"p" + str(t): v for t, v in settings.items()},
            signs=(sign,),
        )
        counts[str(sign)] = {}
        for arm in ("dense", "overlap"):
            states = sum(
                local[t]["dense" if arm == "dense" else "p" + str(t)] for t in settings
            ) / len(settings)
            prob = np.trace(states, axis1=-2, axis2=-1).real
            counts[str(sign)][arm] = np.reshape(
                lab.counts(prob.ravel(), rng), (7, 2)
            ).tolist()
    return {"ports": ports, "seed": SEED, "starts": lab.SHOTS, "counts": counts}


def robust_negative(cols, target):
    # For each fixed orthant find one nonnegative setting whose least
    # available shifted bulk energy strictly exceeds every modular reading.
    keys = list(cols)
    names = list(target)
    matrix = np.array([[cols[k][n][0] for k in keys] for n in names])
    right = np.array([target[n][1] for n in names])
    unit = np.maximum(np.maximum(abs(matrix).max(axis=1), abs(right)), 1e-28)
    res = linprog(
        np.zeros(len(keys)),
        A_ub=-matrix / unit[:, None],
        b_ub=-right / unit - 1e-4,
        bounds=(0, None),
        method="highs",
    )
    if not res.success:
        return None
    from fractions import Fraction as F

    weights = {k: F(float(max(0, v))) for k, v in zip(keys, res.x, strict=True)}
    margin = min(
        sum(w * F(cols[k][n][0]) for k, w in weights.items()) - F(target[n][1])
        for n in names
    )
    if margin <= 0:
        return None
    return {
        "weights": {k: [w.numerator, w.denominator] for k, w in weights.items() if w},
        "strict_margin": [margin.numerator, margin.denominator],
    }


def decode(data):
    iv.dps = 35
    cap = dist.read("local_cap")["balls"]
    energy = dist.read("flight_energy")["result"]["groups"]["mixture"][
        "late_minus_dense"
    ]["overlap"]
    # At fixed radius the six fillings share the same local source and denominator.
    first = {str(row["port"]): row for row in cap[:6]}
    raw = {}
    for sign, table in data["counts"].items():
        arms = {}
        for arm, tape in table.items():
            x = np.array(tape, dtype=float) / data["starts"]
            error = reg.src.probability_error(x, data["starts"])
            arms[arm] = [
                [reg.band(a - b, a + b) for a, b in zip(row, e, strict=True)]
                for row, e in zip(x, error, strict=True)
            ]
        for b in range(2):
            raw[f"{sign}/{b}"] = {
                str(site): region.density.DOSE
                / 2
                * (arms["overlap"][i][b] - arms["dense"][i][b])
                for i, site in enumerate(data["ports"])
            }
    cols = {k: {n: lab.six.bounds(v) for n, v in z.items()} for k, z in raw.items()}

    def sums(keys):
        return {n: lab.six.bounds(sum(raw[k][n] for k in keys)) for n in first}

    modes = {
        "common": {"all": sums(list(raw))},
        "parity": {
            str(b): sums([k for k in raw if k.endswith("/" + str(b))]) for b in (0, 1)
        },
        "sign_parity": cols,
    }
    output = {}
    for scale in ("fixed", "free"):
        target = {
            n: lab.six.bounds(
                iv.mpf(z["readings"]["overlap"]["modular"])
                - (iv.pi / 2 * iv.mpf(energy[n]) if scale == "fixed" else 0)
            )
            for n, z in first.items()
        }
        groups = {}
        for name, axes in modes.items():
            tests = {}
            for signs in product((1, -1), repeat=len(axes)):
                c = {
                    k: (z if sg == 1 else {n: negate(v) for n, v in z.items()})
                    for (k, z), sg in zip(axes.items(), signs, strict=True)
                }
                if scale == "free":
                    c["flight_energy"] = {n: energy[n] for n in first}
                exact = dist.service.separate(c, target, list(first))
                nonpos = {
                    **c,
                    **{
                        "surplus/" + n: {
                            j: [-1.0, -1.0] if n == j else [0.0, 0.0] for j in first
                        }
                        for n in first
                    },
                }
                tests["".join("+" if sg == 1 else "-" for sg in signs)] = {
                    "zero": exact,
                    "nonpositive": dist.service.separate(nonpos, target, list(first)),
                    "strictly_negative_witness": robust_negative(c, target),
                }
            groups[name] = tests
        output[scale] = {"target": target, "groups": groups}
    return {
        "offset_columns": modes,
        "energy": {n: energy[n] for n in first},
        "tests": output,
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    r = acquire()
    r.update(decode(r))
    r["inputs"] = {
        n: sha256((HERE / n).read_bytes()).hexdigest()
        for n in (
            "CAP_ENERGY_ZERO.md",
            "distributed_source.json",
            "six_meter_pointer.json",
            "local_cap.json",
            "flight_energy.json",
        )
    }
    path = HERE / "cap_energy_zero.json"
    if args.write:
        path.write_text(json.dumps(r, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), r)
    for scale, result in r["tests"].items():
        for name, tests in result["groups"].items():
            print(
                scale,
                name,
                len(tests),
                sum(z["zero"]["status"] == "excluded" for z in tests.values()),
                sum(z["nonpositive"]["status"] == "excluded" for z in tests.values()),
                sum(z["strictly_negative_witness"] is not None for z in tests.values()),
            )
