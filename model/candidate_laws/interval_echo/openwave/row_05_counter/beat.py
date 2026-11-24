"""Independent static return-phase calibration and a held precession read."""

import argparse
import json
from pathlib import Path
import numpy as np

from model.candidate_laws.interval_echo.openwave.electron_field.replay import compare
from .dispersive import HERE, SHOTS, prepare, release, probe_counts
from .experiment import draw

CAL = np.arange(0, 16385, 512)
HELD = (20000, 32768, 40000, 65536, 96000)


def spectral_products(basis):
    branches = []
    for m in (2, -2):
        i = 5 + m % 5
        sel, angle, v = basis[(1 + m) % 5]
        w = abs(v[np.where(sel == i)[0][0]]) ** 2
        branches.append((angle, w))
    return branches


def calibration_counts(rng, a, b):
    ident = (abs(a) ** 2 + abs(b) ** 2) / 2
    z = a * b.conjugate()
    return [
        draw(
            np.array(
                [(ident + x) / 2, (ident - x) / 2, (1 - ident) / 2, (1 - ident) / 2]
            ),
            rng,
        )
        for x in (z.real, z.imag)
    ]


def calibrate(basis, rng):
    branches = spectral_products(basis)
    # A controlled product of the two evolutions, with inverse phase on the
    # second reference copy, gives the difference of bare return phases.
    a, b = [np.exp(1j * np.outer(CAL, x)) @ w for x, w in branches]
    counts = np.array(
        [calibration_counts(rng, x, y) for x, y in zip(a, b, strict=True)]
    )
    z = (counts[:, :, 0] - counts[:, :, 1] + counts[:, :, 2] - counts[:, :, 3]) / SHOTS
    slope = float(np.polyfit(CAL, np.unwrap(np.angle(z[:, 0] + 1j * z[:, 1])), 1)[0])
    dominant = tuple(x[np.argmax(w)] for x, w in branches)
    return counts.tolist(), slope, float(dominant[0] - dominant[1])


def run():
    rng = np.random.default_rng(511243)
    out = []
    for p in (-1, 1):
        for d in (1 / 512, 1 / 1024):
            basis, _, bound = prepare(p, d)
            cal, frequency, exact = calibrate(basis, rng)
            held = []
            for n in HELD:
                arms = []
                for sign in (1, -1):
                    a = release(basis, sign, n)
                    bins = probe_counts(a, rng)
                    xy = (bins[:, :, 0] - bins[:, :, 1]).sum(axis=1) / SHOTS
                    z = complex(*xy)
                    forecast = n * frequency * sign
                    arms.append(
                        dict(
                            sign=sign,
                            counts=bins.tolist(),
                            phase_error=np.angle(z * np.exp(-1j * forecast)),
                            visibility=abs(z),
                            source_retention=bins[:, sign, :].sum() / (2 * SHOTS),
                        )
                    )
                held.append(dict(periods=n, arms=arms))
            out.append(
                dict(
                    probe=p,
                    duration=d,
                    calibration_counts=cal,
                    measured_frequency=frequency,
                    spectral_frequency=exact,
                    product_survival_bound=bound,
                    held=held,
                )
            )
    return dict(shots=SHOTS, calibration_periods=CAL.tolist(), studies=out)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    data = run()
    if args.output:
        with args.output.open("x") as f:
            json.dump(data, f, indent=2, allow_nan=False)
            f.write("\n")
    else:
        compare(
            json.loads(json.dumps(data)),
            json.loads((HERE / "beat.json").read_text()),
            "beat",
        )
        print(
            "Verified calibrated signed frequencies and all held standing-wave rotations."
        )


if __name__ == "__main__":
    main()
