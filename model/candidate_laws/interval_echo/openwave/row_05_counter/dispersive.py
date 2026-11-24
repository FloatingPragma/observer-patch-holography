"""Can weak native interring contact rotate a probe while the source retains its motion?"""

import argparse
import json
from pathlib import Path
import numpy as np
from scipy.linalg import expm, schur

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.electron_field.replay import compare
from .experiment import HERE, D, F, SHOTS, board_h, draw


def apparatus(probe):
    h, *_ = board_h(probe, candidate.KAPPA, "all")
    blank, *_ = board_h(probe, candidate.KAPPA, "cut")
    contact = h - blank
    contact -= np.eye(D * D) * np.trace(contact) / (D * D)
    eye = np.eye(D)
    rim = eye - (np.roll(eye, 1, axis=0) + np.roll(eye, -1, axis=0)) / 2
    angle = np.diag(F.conj().T @ rim @ F).real
    phases = (angle[:, None] / 4 + angle[None, :] / 8).ravel()
    free = expm(-1j * np.kron(rim, eye) / 4) @ expm(-1j * np.kron(eye, rim) / 8)
    return phases, free, contact


def prepare(probe, duration):
    phases, free, contact = apparatus(probe)
    pos = np.arange(D**2)
    mom = (pos // D + pos % D) % D
    fourier = np.kron(F, F)
    v = fourier.conj().T @ contact @ fourier
    u = fourier.conj().T @ expm(-1j * duration * contact) @ free @ fourier
    assert np.max(np.abs(u) * (mom[:, None] != mom[None, :])) < 1e-12
    basis = []
    for m in range(D):
        keep = np.flatnonzero(mom == m)
        t, s = schur(u[np.ix_(keep, keep)], output="complex")
        assert np.max(abs(t - np.diag(np.diag(t)))) < 2e-12
        basis.append((keep, np.angle(np.diag(t)), s))
    forecasts = {}
    bounds = []
    for source in (1, -1):
        values = []
        for target in (2, -2):
            i = (source % D) * D + target % D
            m = mom[i]
            rest = np.flatnonzero((pos != i) & (mom == m))
            delta = phases[i] - phases[rest]
            assert np.min(abs(np.sin(delta / 2))) > 1e-5
            v2 = np.sum(abs(v[rest, i]) ** 2 / np.tan(delta / 2)) / 2
            values.append(phases[i] + duration * v[i, i].real + duration**2 * v2)
            keep, _, s = basis[m]
            weight = max(abs(s[np.where(keep == i)[0][0], :]) ** 2)
            bounds.append(max(2 * weight - 1, 0) ** 2)
        forecasts[source] = values[0] - values[1]
    return basis, forecasts, float(min(bounds))


def release(basis, source, n):
    initial = np.zeros((D, D), complex)
    initial[source, [2, -2]] = 1 / np.sqrt(2)
    out = np.zeros(D * D, complex)
    for keep, angles, s in basis:
        out[keep] = s @ (np.exp(1j * n * angles) * (s.conj().T @ initial.ravel()[keep]))
    return out.reshape(D, D)


def probe_counts(a, rng):
    plus, minus = a[:, 2], a[:, -2]
    other = np.sum(abs(a[:, [0, 1, 4]]) ** 2, axis=1)
    return np.array(
        [
            draw(
                np.stack(
                    (
                        abs(plus + c * minus) ** 2 / 2,
                        abs(plus - c * minus) ** 2 / 2,
                        other,
                    ),
                    axis=1,
                ),
                rng,
            )
            for c in (1, 1j)
        ]
    )


def run():
    rng = np.random.default_rng(511242)
    studies = []
    for probe in (-1, 1):
        for denominator in (512, 1024):
            duration = 1 / denominator
            basis, forecasts, bound = prepare(probe, duration)
            for active in (True, False):
                eigenspaces = basis if active else prepare(probe, 0.0)[0]
                enrolled = []
                for target in (np.pi / 4, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi):
                    n = round(target / abs(forecasts[1]))
                    signs = []
                    for source in (1, -1):
                        a = release(eigenspaces, source, n)
                        x = probe_counts(a, rng)
                        diff = (x[:, :, 0] - x[:, :, 1]).sum(axis=1) / SHOTS
                        z = complex(*diff)
                        prediction = -n * forecasts[source] if active else 0.0
                        signs.append(
                            dict(
                                source=source,
                                counts=x.tolist(),
                                predicted_phase=prediction,
                                phase_error=np.angle(z * np.exp(-1j * prediction)),
                                observed_phase=np.angle(z),
                                visibility=abs(z),
                                retention=float(x[:, source, :].sum() / (2 * SHOTS)),
                            )
                        )
                    enrolled.append(dict(target=target, periods=n, signs=signs))
                studies.append(
                    dict(
                        probe=probe,
                        active=active,
                        contact_duration=duration,
                        period=3 / 8 + (duration if active else 0),
                        forecast_per_cycle=forecasts[1],
                        all_period_product_survival=bound if active else 1.0,
                        times=enrolled,
                    )
                )
    return dict(shots=SHOTS, studies=studies)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    data = run()
    if args.output:
        with args.output.open("x") as stream:
            json.dump(data, stream, indent=2, allow_nan=False)
            stream.write("\n")
    else:
        compare(
            json.loads(json.dumps(data)),
            json.loads((HERE / "dispersive.json").read_text()),
            "dispersive",
        )
        print(
            "Verified signed weak-contact phase forecasts, leakage bounds and held detector counts."
        )


if __name__ == "__main__":
    main()
