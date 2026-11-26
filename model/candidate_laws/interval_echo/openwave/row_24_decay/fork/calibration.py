"""Compare the two charged levels to the filled source by summing term readings."""

import argparse
import json
import numpy as np

from .measure import SHOTS, HERE, SFRAME, compare_records
from .material import edge, mount, formation


def acquire():
    board, dimers, heads, _ = mount()
    for wing in (1, 2):
        heads[edge(dimers[0][wing - 1], dimers[wing][0])] = 32
    rng = np.random.default_rng(241111)
    materials = []
    for choice in (0, 1, None):
        readings = []
        words = [(0, 1), (1, 0)] if choice is not None else [(1, 1)]
        for a, b in words:
            rails = np.ones(board.n, dtype=int)
            rails[list(dimers[0])] = [a, b]
            for first, _ in dimers[1:]:
                rails[first] = 0
            cost = int(np.dot(8 * rails**2, 6 - board.degree))
            cost += sum(
                k * (1 - int(rails[p] * rails[q])) for (p, q), k in heads.items()
            )
            # Receiver blank energies are omitted on both stocks and cancel.
            readings.append(cost + 8 * (a - b) ** 2)
        if choice is None:
            n = [SHOTS]
            odd = 0
        else:
            n = rng.multinomial(SHOTS, SFRAME[:, choice] ** 2).tolist()
            odd = int(
                rng.binomial(SHOTS, (SFRAME[0, choice] - SFRAME[1, choice]) ** 2 / 2)
            )
        materials.append(
            dict(
                word_counts=n,
                source_blank_odd=odd,
                word_cost=readings,
                mean=float((np.dot(n, readings) + odd) / SHOTS),
            )
        )
    baseline = materials[2]["mean"]
    measured = [r["mean"] - baseline for r in materials[:2]]
    radius = float(3 * np.sqrt(np.log(8 / 0.01) / (2 * SHOTS)))
    np.testing.assert_allclose(
        measured, np.linalg.eigvalsh(formation(board, dimers, heads)), atol=radius
    )
    charges = []
    for d in dimers[1:]:
        assert all(board.degree[p] == 5 for p in d)
        records = {}
        for a, b in ((0, 1), (1, 0), (0, 0), (1, 1)):
            phase = np.exp(-2j * np.pi * np.arange(8) * (a * a + b * b - 1) / 8)
            amps = np.fft.ifft(phase)
            counts = rng.multinomial(SHOTS, abs(amps) ** 2).tolist()
            records[str(a) + str(b)] = counts
        charges.append(dict(ports=d, charge_bins=records))
    return dict(
        shots=SHOTS,
        raw=materials,
        source_energy_above_filled=measured,
        simultaneous_radius=radius,
        distal_charge_calibration=charges,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--save", action="store_true")
    args = parser.parse_args()
    result = acquire()
    dest = HERE / "calibration.json"
    if args.save:
        dest.write_text(json.dumps(result, indent=1) + "\n")
    else:
        compare_records(json.loads(dest.read_text()), json.loads(json.dumps(result)))
    print("Verified independent contact-energy and charge records.")
