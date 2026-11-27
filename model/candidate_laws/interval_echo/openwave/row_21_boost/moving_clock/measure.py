"""Original launch counts for packet motion, elapsed phase and returned signal."""

import argparse
import json
from pathlib import Path

import numpy as np
from . import apparatus as lab

HERE = Path(__file__).resolve().parent
SHOTS = 2**28
TARGETS = (0, 0.25, -0.25, 0.65, -0.65, 1.05, -1.05)


def histories(z, gates):
    bank = [z]
    for _ in range(lab.T):
        z = lab.step(z, gates)
        bank.append(z)
    return np.array(bank)


def phases(history, shots, rng, mixture=False):
    result = []
    for a, b in history:
        rows = []
        for q in (1, 1j, -1, -1j):
            prob = (
                (abs(a) ** 2 + abs(b) ** 2) / 2 if mixture else abs(a + q * b) ** 2 / 2
            )
            rows.append(lab.draw(prob, shots, rng).reshape(prob.shape))
        result.append(rows)
    return np.array(result, dtype=np.int32)


def panel(n, h, exchange, rng, bank, key):
    layers = lab.gates(exchange, 8.0)
    z = histories(lab.arms(n, h, layers, exchange, 8.0), layers)
    pilots = np.array(
        [
            [lab.draw(abs(ket) ** 2, SHOTS, rng).reshape(n, 2) for ket in z[t]]
            for t in (0, lab.T)
        ],
        dtype=np.int32,
    )
    bank[key + "_pilot"] = pilots
    centers = [lab.address(np.sum(t, axis=(0, 2))) for t in pilots]
    if n // 2 < abs(centers[1] - centers[0]):
        centers[1] += n if centers[1] < 0 else -n
    stations = np.rint(np.linspace(*centers, lab.T + 1)).astype(int)
    bank[key + "_packet"] = phases(z, SHOTS, rng)
    plane = histories(lab.arms(n, h, layers, exchange, 8.0, plane=True), layers)
    bank[key + "_plane"] = phases(plane, SHOTS, rng)
    bank[key + "_incoherent"] = phases(z, SHOTS, rng, mixture=True)
    radar = lab.radar(n, *centers, layers, SHOTS, rng)
    hand = lab.reference_hand(n, [radar[a] for a in ("emission", "turn", "receipt")])
    return dict(
        key=key,
        n=n,
        h=h,
        rail=int(h < 0),
        exchange=exchange,
        pilot_centers=centers,
        stations=stations.tolist(),
        radar=radar,
        hand_positions=hand,
        dispersion=lab.exact(n, h, exchange),
    )


def acquisition():
    rng = np.random.default_rng(211942)
    banks = {}
    panels = []
    routes = {}
    for level in (1, 2):
        route = lab.mount(level)
        n = len(route) // 2
        routes[level] = route
        for exchange, name in ((1.0, "native"), (0.8, "changed")):
            for target in TARGETS:
                h = round(target * n / (2 * np.pi))
                panels.append(panel(n, h, exchange, rng, banks, f"{level}_{name}_{h}"))
    # Null law: both local gates are the same up to their one-double phase,
    # independently of state, momentum and owner count.
    errors = {}
    for x in (1.0, 0.8):
        errors[str(x)] = [
            float(np.abs(a - np.exp(-32j * t) * b).max())
            for a, b, t in zip(
                lab.gates(x, 8.0), lab.gates(x, 0.0), lab.HOLDS, strict=True
            )
        ]
    return dict(
        shots=SHOTS,
        seed=211942,
        routes=routes,
        forecast_periods=lab.T,
        rows=panels,
        zero_contrast_gate_errors=errors,
    ), banks


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    path = parser.parse_args().output
    data, banks = acquisition()
    path.with_suffix(".json").write_text(
        json.dumps(data, indent=1, allow_nan=False) + "\n"
    )
    np.savez_compressed(path.with_suffix(".npz"), **banks)
