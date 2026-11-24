"""Material flux spectroscopy of a free density wall on the native belt."""

import argparse
import json
from pathlib import Path

import numpy as np

from model.candidate_laws.contrast_memory_spherical_v1.sectors import BinarySector
from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.electron_field.replay import compare
from model.candidate_laws.interval_echo.openwave.quark_winding.belt_walls import belt

HERE = Path(__file__).parent
LENGTH = 10
TIMES = np.arange(33)
SHOTS = 262144
CHARGES = [-1 / 3, 2 / 3, 0, -1, 1]


class Ring:
    def __init__(self, amount, contrast):
        sec = BinarySector(belt.ring(LENGTH), amount)
        self.h = sec.hamiltonian(contrast)
        self.wall = 0 if amount == 1 else int(belt.selected(sec)[0])
        sites = list(map(set, sec.states))
        self.d = np.zeros_like(self.h)
        rows, cols = np.nonzero(np.triu(self.h, k=1))
        for i, j in zip(rows, cols):
            (new,) = sites[i] - sites[j]
            (old,) = sites[j] - sites[i]
            step = (new - old + LENGTH // 2) % LENGTH - LENGTH // 2
            if abs(step) not in (1, 2):
                raise ValueError("Transition is not a belt contact")
            self.d[i, j], self.d[j, i] = step, -step
        self.orbits = []
        unseen = set(range(len(sites)))
        while unseen:
            word = min(unseen)
            orbit = []
            while word in unseen:
                unseen.remove(word)
                orbit.append(word)
                shifted = tuple(sorted((j + 1) % LENGTH for j in sites[word]))
                word = sec.index[shifted]
            self.orbits.append(orbit)
        self.momenta = [self.sector(m) for m in range(LENGTH)]
        self.unitary = np.exp(2j * np.pi * np.array([sum(s) for s in sites]) / LENGTH)

    def sector(self, m):
        vectors = []
        source = None
        for orbit in self.orbits:
            size = len(orbit)
            if (m * size) % LENGTH:
                continue
            vec = np.zeros(len(self.h), complex)
            vec[orbit] = np.exp(2j * np.pi * m * np.arange(size) / LENGTH) / np.sqrt(
                size
            )
            if self.wall in orbit:
                source = len(vectors)
            vectors.append(vec)
        if source is None:
            raise ValueError("Missing translated wall")
        return np.column_stack(vectors), source

    def generator(self, flux):
        return self.h * np.exp(1j * flux * self.d / LENGTH)

    def line(self, flux, m):
        basis, source = self.momenta[m]
        block = basis.conj().T @ self.generator(flux) @ basis
        energy, vec = np.linalg.eigh(block)
        weights = abs(vec[source]) ** 2
        starts = np.r_[0, np.flatnonzero(np.diff(energy) > 1e-10) + 1]
        line_weights = np.add.reduceat(weights, starts)
        line = int(np.argmax(line_weights))
        i = starts[line]
        centered = energy - self.h[self.wall, self.wall]
        waves = (np.exp(-1j * np.outer(TIMES, centered)) * vec[source].conj()) @ vec.T
        return float(centered[i]), float(line_weights[line]), waves, source


def ramsey(amplitude, rng):
    quadratures = np.array([amplitude.real, amplitude.imag])
    prob = np.clip((1 + quadratures) / 2, 0, 1)
    counts = rng.binomial(SHOTS, np.round(prob, 12))
    real, imag = 2 * counts / SHOTS - 1
    phase = np.unwrap(np.angle(real + 1j * imag))
    frequency = -np.polyfit(TIMES, phase, 1)[0]
    return {
        "real_imag_plus_counts": counts.tolist(),
        "frequency": float(frequency),
        "minimum_visibility": float(min(abs(real + 1j * imag))),
    }


def interpolate(values, step, angle):
    """Degree-five trigonometric interpolation with a real Nyquist coefficient."""
    idx = (step * np.arange(LENGTH)) % LENGTH
    ordered = np.empty(LENGTH)
    ordered[idx] = values
    coeff = np.fft.rfft(ordered) / LENGTH
    modes = np.arange(1, LENGTH // 2)
    return float(
        coeff[0].real
        + 2 * np.real(coeff[modes] @ np.exp(1j * modes * angle))
        + coeff[-1].real * np.cos(LENGTH * angle / 2)
    )


def study(load, kappa, rng):
    ring = Ring(load, kappa)
    gauge_residual = np.max(
        abs(
            ring.generator(2 * np.pi)
            - ring.unitary[:, None] * ring.h * ring.unitary.conj()
        )
    )
    if gauge_residual > 1e-12:
        raise ArithmeticError("Flux quantum differs from a unitary change of gauge")
    step = 1 if load == 1 else 3
    calibration = []
    for m in range(LENGTH):
        energy, weight, waves, site = ring.line(0, m)
        read = ramsey(waves[:, site], rng)
        calibration.append(
            dict(m=m, dominant_energy=energy, dominant_weight=weight, read=read)
        )
    zero_read = [r["read"]["frequency"] for r in calibration]
    zero_exact = [r["dominant_energy"] for r in calibration]
    panels = []
    for m in (1, 2, 3, 4):
        center = step * 2 * np.pi * m / LENGTH
        for phi in (np.pi / 4, np.pi / 2, 2 * np.pi / 3):
            ep, wp, plus, _ = ring.line(phi, m)
            em, wm, minus, _ = ring.line(-phi, m)
            overlap = np.sum(minus.conj() * plus, axis=1)
            read = ramsey(overlap, rng)
            forecasts = []
            for q in CHARGES:
                shift = step * q * phi / LENGTH
                estimates = [
                    interpolate(d, step, center - shift)
                    - interpolate(d, step, center + shift)
                    for d in (zero_read, zero_exact)
                ]
                forecasts.append(
                    {
                        "charge": q,
                        "forecast": estimates[0],
                        "noiseless_calibration": estimates[1],
                    }
                )
            panels.append(
                dict(
                    m=m,
                    flux=phi,
                    dominant_splitting=ep - em,
                    plus_minus_weights=[wp, wm],
                    read=read,
                    forecasts=forecasts,
                )
            )
    errors = []
    for i, q in enumerate(CHARGES):
        delta = [v["forecasts"][i]["forecast"] - v["read"]["frequency"] for v in panels]
        exact = [
            v["forecasts"][i]["noiseless_calibration"] - v["dominant_splitting"]
            for v in panels
        ]
        errors.append(
            {
                "charge": q,
                "maximum_measured_error": max(abs(np.array(delta))),
                "rms_measured_error": float(np.sqrt(np.mean(np.square(delta)))),
                "maximum_spectral_error": max(abs(np.array(exact))),
            }
        )
    return dict(
        load=load,
        contrast=kappa,
        gauge_residual=float(gauge_residual),
        calibration=calibration,
        held=panels,
        errors=errors,
    )


def run():
    axes = belt.belts()
    rng = np.random.default_rng(120436)
    cases = (
        (1, candidate.KAPPA),
        (3, candidate.KAPPA),
        (4, candidate.KAPPA),
        (3, 0),
        (4, 0),
    )
    return {
        "belt": next(iter(axes.values())),
        "equivalent_axes": len(axes),
        "shots_per_quadrature": SHOTS,
        "times": TIMES.tolist(),
        "studies": [study(n, k, rng) for n, k in cases],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    data = run()
    if args.output:
        with args.output.open("x") as f:
            json.dump(data, f, indent=2, allow_nan=False)
            f.write("\n")
        return
    expected = json.loads((HERE / "results.json").read_text())
    compare(json.loads(json.dumps(data, allow_nan=False)), expected, "material_flux")
    print(
        "Reproduced single-unit calibration, both wall spectra and the changed-law readings."
    )


if __name__ == "__main__":
    main()
