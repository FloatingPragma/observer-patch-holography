"""Ramsey interference of an independent positive or negative material probe on written heads."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.interval_echo import candidate, currents
from . import polarization as pol


def contact_calibration():
    h0, h1 = law.local_term_parts()
    h = h0 + candidate.KAPPA * h1
    cols = []
    for c in (0, 1):
        head = np.zeros(11)
        for d in (-1, 1):
            head[1 + law.LABELS.index((c, d))] = 1 / np.sqrt(2)
        for sign in (-1, 1):
            for direction in (1, -1):
                pair = (1 + sign, 1) if direction == 1 else (1, 1 + sign)
                idle = np.kron(head, np.eye(9)[law.SRC.index((1, 1))])
                active = np.kron(head, np.eye(9)[law.SRC.index(pair)])
                permitted = int(c == int(sign * direction > 0))
                for t in (0.03125, 0.125, np.pi / 8):
                    u = expm(-1j * t * h)
                    amp = np.vdot(active, u @ active) / np.vdot(idle, u @ idle)
                    if t == np.pi / 8:
                        expected_head = np.zeros(11)
                        for d in (-1, 1):
                            expected_head[
                                1 + law.LABELS.index((1 - c if permitted else c, d))
                            ] = 1 / np.sqrt(2)
                        target = np.kron(
                            expected_head,
                            np.eye(9)[law.SRC.index(pair[::-1] if permitted else pair)],
                        )
                        assert abs(abs(np.vdot(target, u @ active)) - 1) < 1e-12
                    predicted = np.exp(4j * t) * np.cos(4 * t) if permitted else 1
                    assert abs(amp - predicted) < 1e-12
                    cols.append(
                        dict(
                            coin=c,
                            sign=sign,
                            direction=direction,
                            t=t,
                            permitted=permitted,
                            stay_real=amp.real,
                            stay_imag=amp.imag,
                            departure=1 - abs(amp) ** 2,
                        )
                    )
    return cols


def interferometer(before, after, path, f, theta, sign, cycles=8):
    grad = np.diff(f)
    old = np.where(np.diff(path) > 0, before, 1 - before)
    new = np.where(np.diff(path) > 0, after, 1 - after)
    return edge_response(old, new, grad, theta, sign, cycles)


def edge_response(old, new, grad, theta, sign, cycles=8):
    """Nested native return amplitudes for independently addressed oriented edges."""
    if sign < 0:
        old, new = 1 - old, 1 - new
    sense_a = np.where(grad >= 0, old, new)
    sense_b = np.where(grad >= 0, new, old)
    angle = theta * abs(grad) / cycles
    active = np.cos(angle) * np.exp(1j * angle)
    a = np.prod(np.where(sense_a, active, 1)) ** cycles
    b = np.prod(np.where(sense_b, active, 1)) ** cycles
    coh = a * b.conjugate()
    returns = (abs(a) ** 2 + abs(b) ** 2) / 2
    loss_a = (1 - abs(a) ** 2) / 2
    loss_b = (1 - abs(b) ** 2) / 2
    panels = {}
    for k, v in [("x", coh.real), ("y", coh.imag)]:
        panels[k] = [(returns + v) / 2, (returns - v) / 2, loss_a, loss_b]
        np.testing.assert_allclose(sum(panels[k]), 1, atol=1e-11)
    return dict(phase=np.angle(coh), visibility=abs(coh), quadratures=panels)


def run():
    path, _, _ = pol.layout()
    initial = pol.common_heads(path)
    cases = {}
    for shift in (1, -2):
        _, _, flux, original = pol.translated(path, shift)
        final = initial ^ (flux != 0).astype(int)
        geom = np.sign(np.diff(path))
        np.testing.assert_array_equal(geom * (initial - final), flux)
        panels = []
        for name, centre, plateau in [
            ("left", 0, 8),
            ("right", 120, 8),
            ("outside", -54, 8),
            ("pair", 60, 70),
        ]:
            for offset in (0, 1):
                for width in (0, 1, 2, 4, 8, 16, 32):
                    f = pol.potential(centre + offset, width, plateau)
                    for probe in (1, -1):
                        rows = [
                            interferometer(initial, final, path, f, theta, probe)
                            for theta in (np.pi / 2, 2 * np.pi)
                        ]
                        panels.append(
                            dict(
                                region=name,
                                offset=offset,
                                width=width,
                                probe=probe,
                                quarter_turn=rows[0],
                                full_turn=rows[1],
                            )
                        )
        cases[str(shift)] = dict(
            flipped_heads=int(np.sum(initial != final)), reads=panels
        )
    # One ordinary material unit, initially outside the electrode, ends on its plateau.
    f = pol.potential(0, 32)
    unit_path = (pol.X[:-1] >= -68) & (pol.X[:-1] < 0)
    unit_initial = (np.diff(path) > 0).astype(int)
    unit_final = unit_initial ^ unit_path.astype(int)
    calibrations = {}
    for name, b, a in [
        ("whole", unit_initial, unit_final),
        ("sham", unit_initial, unit_initial),
        ("head_independent_wall", initial, initial),
    ]:
        calibrations[name] = [
            interferometer(b, a, path, f, t, 1) for t in (np.pi / 2, 2 * np.pi)
        ]
    mixtures = {}
    for fraction in (1 / 3, 2 / 3):
        mixtures[str(fraction)] = {
            k: (
                (1 - fraction) * np.array(calibrations["sham"][1]["quadratures"][k])
                + fraction * np.array(calibrations["whole"][1]["quadratures"][k])
            ).tolist()
            for k in ("x", "y")
        }
    # Independent shot ensembles at an aperture family absent from the width scan.
    shoulder = np.clip((abs(pol.X + 1) - 9) / 23, 0, 1)
    held = 1 - 3 * shoulder**2 + 2 * shoulder**3
    settings = {"unit": (unit_initial, unit_final)}
    for shift in (1, -2):
        _, _, flux, _ = pol.translated(path, shift)
        settings[str(shift)] = (initial, initial ^ (flux != 0).astype(int))
    rng = np.random.default_rng(1201313)
    shots = 262144
    count = 2 * 2 * len(settings)
    delta = np.sqrt(2 * np.log(2 * count / 0.01) / shots)
    samples = {}
    for name, (before, after) in settings.items():
        for sign in (1, -1):
            r = interferometer(before, after, path, held, np.pi / 2, sign)
            hist = {
                ax: rng.multinomial(shots, bins)
                for ax, bins in r["quadratures"].items()
            }
            z = complex(*[(bins[0] - bins[1]) / shots for bins in hist.values()])
            err = np.arcsin(min(1.0, np.sqrt(2) * delta / abs(z)))
            samples[f"{name}_{sign}"] = dict(
                counts={k: v.tolist() for k, v in hist.items()},
                phase=np.angle(z),
                phase_radius=err,
                expected_phase=r["phase"],
            )
    for sign in (1, -1):
        unit = samples[f"unit_{sign}"]
        lo, hi = sorted(
            [
                sign * (unit["phase"] + d)
                for d in (-unit["phase_radius"], unit["phase_radius"])
            ]
        )
        for name in ("1", "-2"):
            entry = samples[f"{name}_{sign}"]
            ratios = [
                sign * (entry["phase"] + d) / den
                for d in (-entry["phase_radius"], entry["phase_radius"])
                for den in (lo, hi)
            ]
            entry["fraction_of_measured_unit"] = entry["phase"] / unit["phase"]
            entry["simultaneous_interval"] = [min(ratios), max(ratios)]
    return currents._portable(
        dict(
            contact=contact_calibration(),
            cycles=8,
            ports=["bright", "dark", "escape_a", "escape_b"],
            held=dict(shots=shots, seed=1201313, score_error=delta, samples=samples),
            common_heads=initial.tolist(),
            cases=cases,
            calibrations=calibrations,
            intermittent_whole_full_turn=mixtures,
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    encoded = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    with args.output.open("x") as f:
        f.write(encoded)
