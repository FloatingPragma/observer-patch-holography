"""Blind short-time test of the scalar graph continuation with a calibrated timer."""

import argparse
import json

import numpy as np
from scipy.sparse.linalg import expm_multiply

from model.executable_model import federation
from .measure import HERE, SHOTS, acquire, band, events, fixture, material


def wall_laplacian(links):
    ends = [{a, b} for a, _, b, _ in links]
    adjacency = np.array(
        [
            [bool(a & b) and i != j for i, a in enumerate(ends)]
            for j, b in enumerate(ends)
        ],
        dtype=float,
    )
    return np.diag(adjacency.sum(axis=0)) - adjacency


def run():
    rng = np.random.default_rng(790302)
    meters = json.loads((HERE.parent / "coupling_meter.json").read_text())["results"]
    mass, b = [meters[key]["value"] for key in ("predicted_m", "predicted_c_squared")]
    result = {"calibrated_mass": mass, "calibrated_c_squared": b, "settings": {}}
    for length in (2, 3):
        _, _, links = fixture(federation(0), length)
        h, walls = material(links, 8)
        spatial, v = np.linalg.eigh(wall_laplacian(links))
        lmax = spatial.max()
        energies = np.sqrt(mass**2 + b * spatial) - mass
        for r in (0, 1):
            word = int(next(x for x in np.flatnonzero(walls == 1) if x & 1 == r))
            initial = np.zeros(len(walls))
            initial[word] = 1.0
            counts = acquire(expm_multiply(-0.0625j * h, initial), rng)
            timer = rng.binomial(SHOTS, np.sin(0.0625 / 8) ** 2)
            t = 8 * np.arcsin(np.sqrt(timer / SHOTS))
            tlo, thi = 8 * np.arcsin(np.sqrt(band(timer)))
            dt = max(thi - t, t - tlo)
            # Uniform phase is immaterial. Bound the gradient of sqrt(m²+bλ)-m
            # throughout the whole parameter rectangle, then Duhamel in time.
            mlo, bhi = 0.99 * mass, 1.01 * b
            variation = 0.01 * mass * bhi * lmax / (2 * mlo**2) + 0.01 * b * lmax / (
                2 * mlo
            )
            error = thi * variation + dt * max(energies)
            psi = v @ (np.exp(-1j * t * energies) * v[0])
            report = {}
            for name in ("along_A", "turn_B", "turn_C"):
                mask = events(length)[name]
                amplitude = abs(psi[mask.bit_length() - 1])
                forecast = [
                    max(0.0, amplitude - error) ** 2,
                    min(1.0, amplitude + error) ** 2,
                ]
                k = int(counts[walls == mask].sum())
                observed = band(k)
                report[name] = {
                    "count": k,
                    "observed": observed,
                    "forecast": forecast,
                    "compatible": bool(
                        max(forecast[0], observed[0]) <= min(forecast[1], observed[1])
                    ),
                }
            result["settings"][f"L{length}_r{r}"] = {
                "counts": counts.tolist(),
                "timer": timer,
                "time_band": [tlo, thi],
                "reads": report,
            }
    return json.loads(json.dumps(result))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    output = HERE / "scalar.json"
    result = run()
    if args.write:
        output.write_text(json.dumps(result, indent=2) + "\n")
    else:
        from model.candidate_laws.interval_echo.openwave.electron_field.replay import (
            compare,
        )

        compare(result, json.loads(output.read_text()), "fork scalar extrapolation")
        print("Replayed the fresh scalar-junction comparison.")


if __name__ == "__main__":
    main()
