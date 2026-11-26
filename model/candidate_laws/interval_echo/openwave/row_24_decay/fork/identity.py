"""Read the fork's charged levels after a held scalar-potential drive and in free return."""

import argparse
import json

import numpy as np
from scipy.linalg import eigh

from .measure import (
    HERE,
    SHOTS,
    SFRAME,
    XFRAME,
    core_pulse,
    local,
    real_counts,
    compare_records,
)
from .material import mount, generator


def run():
    board, d, heads, _ = mount()
    posts = set(sum(d, ()))
    pins = [
        next((s, k) for s, k in sorted(heads.items()) if port in s and set(s) - posts)
        for port in d[0]
    ]
    # Native sign-even/odd active columns at the physical exterior unit.
    word = np.arange(8) >> 2
    costs = [pins[0][1] * (1 - word), pins[1][1] * word]
    cut = generator(board, d, heads, (1, 2))
    e, v = eigh(cut)
    basis = np.kron(SFRAME, np.kron(XFRAME, XFRAME))
    cal = json.loads((HERE / "calibration.json").read_text())
    p = cal["raw"][0]["word_counts"][0] / SHOTS
    calibration = json.loads((HERE / "records.json").read_text())
    gap = calibration["source_clock"]["per_native_time"]
    fitted = np.array([[np.sqrt(p), -np.sqrt(1 - p)], [np.sqrt(1 - p), np.sqrt(p)]])
    rng = np.random.default_rng(241112)

    def receive(start, pulses):
        actual = basis[:, start * 4].astype(complex)
        forecast = fitted[:, start].astype(complex)
        for phase, rail, wait in pulses:
            sites = range(2) if rail == "both" else [rail]
            for k in sites:
                actual *= np.exp(-1j * phase / pins[k][1] * costs[k])
                charged = np.eye(2)[:, k]
                forecast *= np.exp(-1j * phase * charged)
            actual = v @ (np.exp(-1j * e * wait) * (v.T @ actual))
            forecast = fitted @ (
                np.exp(-1j * np.array([0.0, gap]) * wait) * (fitted.T @ forecast)
            )
        hist = real_counts(rng, abs(basis.T @ actual) ** 2)
        flip = abs((fitted.T @ forecast)[1 - start]) ** 2
        return dict(
            start=start,
            pulses=pulses,
            all_core_left_right=hist,
            predicted_conversion=float(flip),
        )

    banks = []
    for start in (0, 1):
        banks.append(receive(start, []))
        for angle in (np.pi / 3, np.pi / 2, np.pi):
            for rail in (0, 1, "both"):
                banks.append(receive(start, [(angle, rail, 0)]))
    for wait in (0, 1, 2, 4, 8):
        for second in (0, 1):
            banks.append(receive(0, [(np.pi / 2, 0, wait), (np.pi / 2, second, 0)]))
    h = generator(board, d, heads)
    e, v = eigh(h)
    ground = v[:, 0]
    pulse = local(core_pulse(SFRAME), 0)
    life = []
    for active in (True, False):
        start = pulse @ ground if active else ground
        for wait in (0, 12, 16, 20, 24, 28, 32, 36):
            state = v @ (np.exp(-1j * wait * e) * (v.T @ start))
            life.append(
                dict(
                    active=active,
                    wait=wait,
                    core_left_right=real_counts(rng, abs(basis.T @ state) ** 2),
                )
            )
    return dict(
        shots=SHOTS,
        seed=241112,
        source_pins=pins,
        forecast_primary_vacancy=p,
        forecast_gap=gap,
        charge_only=banks,
        lifetime=life,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--save", action="store_true")
    args = parser.parse_args()
    result = run()
    dest = HERE / "identity.json"
    if args.save:
        dest.write_text(json.dumps(result, indent=1) + "\n")
    else:
        compare_records(json.loads(dest.read_text()), json.loads(json.dumps(result)))
    print("Verified charge-only conversion and held parent-survival records.")
