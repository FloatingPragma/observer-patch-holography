"""Compare calibrated changes of native seam exposure time to the distributed two-word contrast."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json

import numpy as np

from model.executable_model import federation
from model.experiments.interval_echo_einstein_station import shared_variance, volume

from . import resolved_shells as shells

HERE, iv, bounds = shells.HERE, shells.iv, shells.bounds
STARTS = 2**48
SEED = 28112130
H = 1 / 1024


def response(energies, law, t):
    x = [iv.cos(t * h) for h in energies]
    y = [-iv.sin(t * h) for h in energies]
    moment = {}
    for k, left, right in (("xx", x, x), ("yy", y, y), ("xy", x, y), ("yx", y, x)):
        moment[k] = (
            sum(int(n) * left[a] * right[b] for (a, b), n in np.ndenumerate(law)) / 2592
        )
    return moment


def fringe(means):
    return [means["xx"] - means["yy"], means["xy"] + means["yx"]]


def acquired(quad, admission, rng):
    tape, decoded = {}, {}
    for key, moment in quad.items():
        hits = int(rng.binomial(STARTS, float(admission)))
        yes = int(rng.binomial(hits, (1 + float(moment.mid)) / 2))
        tape[key] = [yes, hits - yes, STARTS - hits]
        p = np.array(tape[key][:2]) / STARTS
        err = shells.reg.src.probability_error(p, STARTS).sum()
        z = (yes - (hits - yes)) / STARTS
        decoded[key] = shells.reg.band(z - err, z + err)
    return tape, fringe(decoded)


def run():
    iv.dps = 35
    rng = np.random.default_rng(SEED)
    src = volume.radar_source_and_hand(federation(0))[-1]
    energies = np.asarray(src["source_word_energies"], dtype=int)
    laws, _ = shared_variance.source_words(src)
    reverse, fresh, shared = np.rint(2592 * laws).astype(int)
    prior = json.loads((HERE / "resolved_shells.json").read_bytes())
    rows = {}
    for weight in (1, 2):
        e = (energies[:, None] + weight * energies).ravel().tolist()
        t = iv.pi / 3
        forecasts = {
            k: response(e, law, time)
            for k, law, time in (
                ("shared", shared, t),
                ("fresh", fresh, t),
                ("early", reverse, t - H),
                ("late", reverse, t + H),
            )
        }
        df = fringe(response(e, shared - fresh, t))
        energy = np.add.outer(e, e)
        # The derivative of exp(-itE) is -i E exp(-itE).
        d_real, d_imag = fringe(response(e, reverse * energy, t))
        rate = [d_imag, -d_real]
        views = []
        for panel in prior["levels"]:
            p = F(panel["eligible_pair_probability"])
            prob = shells.interval(p)
            tape, read = {}, {}
            for name, quad in forecasts.items():
                tape[name], read[name] = acquired(quad, p, rng)
            data = []
            for j, n in enumerate(("real", "imaginary")):
                source = prob * df[j]
                predicted = prob * rate[j]
                held = read["shared"][j] - read["fresh"][j]
                # Real and imaginary derivatives of order three share the
                # spectral norm bound. All misses have identically zero score.
                error = prob * int(energy.max()) ** 3 * H**2 / 6
                measured = (read["late"][j] - read["early"][j]) / (2 * H) + iv.mpf(
                    [-error.b, error.b]
                )
                data.append(
                    {
                        "quadrature": n,
                        "source": bounds(source),
                        "from_measured_modular": bounds(
                            iv.mpf(panel["direct_entropy_tangent"])
                            * df[j]
                            / iv.mpf(prior["conditional_entropy_slope"])
                        ),
                        "held": bounds(held),
                        "native_time_derivative": bounds(predicted),
                        "measured_time_derivative": bounds(measured),
                        "time_equivalent": bounds(source / predicted)
                        if 0 not in predicted
                        else None,
                        "held_time_equivalent": bounds(held / measured)
                        if 0 not in measured
                        else None,
                    }
                )
            views.append({"level": panel["level"], "tape": tape, "readings": data})
        rows[str(weight)] = {
            "conditional_source": list(map(bounds, df)),
            "conditional_time_derivative": list(map(bounds, rate)),
            "boards": views,
        }
    return {
        "seed": SEED,
        "starts": STARTS,
        "exposure": "pi/3",
        "time_step": H,
        "energies": energies.tolist(),
        "clocks": rows,
        "inputs": {
            n: sha256((HERE / n).read_bytes()).hexdigest()
            for n in ("resolved_shells.json", "WORD_PHASE_READ.md")
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    data = run()
    path = HERE / "word_phase_read.json"
    if args.write:
        path.write_text(json.dumps(data, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), data)
    for k, v in data["clocks"].items():
        print(
            k,
            v["conditional_source"],
            v["conditional_time_derivative"],
            [x["held_time_equivalent"] for x in v["boards"][-1]["readings"]],
        )
