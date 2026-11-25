"""Complex transfer to independently prepared opposite-parity transverse receivers."""

import argparse
import json

import numpy as np

from .measure import (
    HERE,
    SHOTS,
    chart,
    draw,
    preparation,
    prepare_target,
    spectrum,
    unit,
)


def acquire():
    rng = np.random.default_rng(540431)
    xyz, axes = chart()
    dipole = preparation("port", 0, True)
    # Poles have zero transverse amplitude. The sign appointment acts on
    # the lower ring; any phase at the two zero-amplitude poles is immaterial.
    even = prepare_target(unit(dipole * np.where(xyz @ axes["port"] < 0, -1, 1)))
    templates = np.stack([dipole, even])
    rows = []
    for scale in (1, 3):
        levels, modes = spectrum("cell", scale)
        for q in (0, 1, 2, 3, 4):
            state = preparation("port", q, True)
            weights = (templates.conj() @ modes) * (modes.T @ state)
            times = np.array([0, 1 / 6, 1 / 2])
            signals = weights @ np.exp(-1j * np.outer(levels, times))
            samples = []
            for template in signals:
                readings = []
                for z in template:
                    readings.append(
                        dict(
                            cos=draw(rng, [(1 + z.real) / 2, (1 - z.real) / 2]),
                            sin=draw(rng, [(1 - z.imag) / 2, (1 + z.imag) / 2]),
                        )
                    )
                samples.append(readings)
            rows.append(dict(q=q, scale=scale, holds=times.tolist(), receivers=samples))
    return dict(shots=SHOTS, banks=rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    dest = HERE / "transverse.json"
    data = json.dumps(acquire(), separators=(",", ":")) + "\n"
    if args.write:
        dest.write_text(data)
    elif dest.read_text() != data:
        raise SystemExit("Transverse clock records differ")
