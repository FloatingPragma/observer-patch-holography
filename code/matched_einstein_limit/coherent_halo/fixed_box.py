"""Normalized, fixed-volume refinement at the independently chosen cubic scattering pole."""

import json
from math import pi

import numpy as np

from .measure import binding, cubic


def panel(side, critical):
    mass = binding(side, critical)
    modes = np.rint(np.fft.fftfreq(side) * side).astype(int)
    split = 4 * np.sin(pi * modes / side) ** 2
    symbol = split[:, None, None] + split[None, :, None] + split[None, None, :]
    spectrum = 1 / (side**2 * (symbol + mass**2))
    field = np.fft.ifftn(spectrum).real * side**3
    # Remove the coincident site on the real finite guide before normalizing.
    removed = field[0, 0, 0] ** 2 / side**3
    field[0, 0, 0] = 0.0
    norm2 = np.mean(field**2)
    np.testing.assert_allclose(norm2 + removed, np.sum(spectrum**2), atol=2e-13)
    r2 = modes**2
    distance_squared = r2[:, None, None] + r2[None, :, None] + r2[None, None, :]
    prob = field**2 / (side**3 * norm2)
    return dict(
        side=side,
        mass_times_side=side * mass,
        field_norm_squared=norm2,
        norm_at_deleted_site=removed,
        enclosed_probabilities={
            str(r): float(prob[distance_squared <= (side * r) ** 2].sum())
            for r in (1 / 16, 1 / 8, 1 / 4)
        },
    )


def run():
    critical = 1 / (12 * (cubic.green(0.0) - 1 / 6))
    return [panel(n, critical) for n in (16, 32, 64, 128)]


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
