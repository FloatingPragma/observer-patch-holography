"""Terminal signed material probe, including all zero-weight and refused trials."""

import numpy as np

from model.experiments.material_inertia.pair_spectrum.native_source_probe import circuit
from model.experiments.material_inertia.pair_spectrum.native_source_probe.measure import (
    routes,
)
from model.executable_model import federation
from .material import vertices

SHOTS = 2**24


def calibrate(rng):
    meters, _ = circuit.audit()
    # Signs are clocked in the 20-versus-02 exit convention.
    response = np.array(
        [meters[1][q]["Y"][6] - meters[1][q]["Y"][2] for q in (-1, 0, 1)]
    )
    assert abs(response[1]) < 1e-11 and abs(response[0] + response[2]) < 1e-11
    direction = 1 if response[2] > 0 else -1
    response *= direction
    standards = [draw(rng, [(1 + x) / 2, (1 - x) / 2]) for x in response]
    path = routes(federation(0), {0})
    return float(response[2]), {
        "original_pointer_sign": direction,
        "counts": standards,
        "physical_route": path,
    }


def draw(rng, weights, shots=SHOTS):
    p = np.maximum(np.real(weights), 0)
    assert abs(sum(p) - 1) < 1e-9
    ticks = np.rint(p * 2**40).astype(np.int64)
    ticks[np.argmax(p)] += 2**40 - sum(ticks)
    return rng.multinomial(shots, ticks / 2**40).tolist()


class Quadrupole:
    def __init__(self, words, occupancy, axes, base, gain):
        v = vertices()
        center = (occupancy @ v) / occupancy.sum(axis=1, keepdims=True)
        self.effects = {}
        q = words - 1  # actual target relative to unit calibration
        for frame in ("cell", "centroid"):
            xyz = (
                v[None, :, :] - (0 * center if frame == "cell" else center)[:, None, :]
            )
            r2 = (xyz**2).sum(axis=2)
            for i, axis in enumerate(axes):
                weight = 3 * (xyz @ axis) ** 2 - r2
                scale = float(np.max(abs(weight)))
                accepted = abs(weight).mean(axis=1) / scale
                data = (weight * q).mean(axis=1) * gain / scale
                blank = weight.mean(axis=1) * base * gain / scale
                self.effects[f"{frame}_{i}"] = (scale, accepted, data, blank)
        self.gain = gain

    def bank(self, p, rng, shots=SHOTS):
        result = {}
        for key, (scale, accepted, data, blank) in self.effects.items():
            a = p @ accepted
            rows = []
            # Two independent arms use matching material centroid and acceptances.
            for diff in (data, blank):
                s = p @ diff
                rows.append(
                    draw(
                        rng, [(a + s) / 2, (a - s) / 2, p.sum() - a, 1 - p.sum()], shots
                    )
                )
            result[key] = {
                "weight": round(scale, 12),
                "material": rows[0],
                "host": rows[1],
            }
        return result


def decode(bank, gain):
    def signed(t):
        return (t[0] - t[1]) / max(1, sum(t[:3]))

    return (
        12 * bank["weight"] * (signed(bank["material"]) - signed(bank["host"])) / gain
    )


def sampler_radius(bank, gain_cal):
    # Hoeffding with a simultaneous budget of 4096 effects. Sign lies in [-1,1].
    def eps(n):
        return np.sqrt(2 * np.log(2 * 4096 / 0.01) / max(1, n))

    low = gain_cal - eps(SHOTS)
    value = decode(bank, gain_cal)
    return float(
        (
            12
            * bank["weight"]
            * (eps(sum(bank["material"][:3])) + eps(sum(bank["host"][:3])))
            + abs(value) * eps(SHOTS)
        )
        / low
    )
