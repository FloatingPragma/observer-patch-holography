"""Acquire randomized preparations through the native seam-sign receiver."""

import argparse
import json
from pathlib import Path

import numpy as np

from ..measure import carrier, gates, preparation, prepare_target

HERE = Path(__file__).parent
N = 2**34
BLOCKS = 8
GRID = 2**32
DIRECTIONS = ("zero", "port", "edge", "face", "oblique")


def exits(psi, edge, decoder):
    return abs((decoder @ psi[list(edge)])[1]) ** 2


def native_probabilities():
    decoder = gates("native")[1][:2, :2]
    for i, j in carrier.CELL_SEAMS:
        for amps, expected in [((1, 1), 0), ((1, -1), 1), ((1, 0), 0.5)]:
            target = np.zeros(12, dtype=complex)
            target[[i, j]] = amps / np.linalg.norm(amps)
            observed = exits(prepare_target(target), (i, j), decoder)
            np.testing.assert_allclose(observed, expected, atol=1e-12, rtol=0)
    signals = {}
    for direction in DIRECTIONS:
        psi = (
            preparation("port", 0, False)
            if direction == "zero"
            else preparation(direction, 4, False)
        )
        signals[direction] = [
            round(exits(psi, edge, decoder) * GRID) / GRID
            for edge in carrier.CELL_SEAMS
        ]
    return signals


def acquire_bank(rng, signal, seam, blind):
    # Each batch has a fresh concealed arm on every shot. Only records from
    # completed batches can influence the next receiver response.
    feedback = 0
    ledger = []
    for batch in range(BLOCKS):
        offset = (1 + (seam + 3 * batch + feedback) % 6) / 16
        visibility = 0 if blind else (3 + (seam + batch + feedback) % 6) / 16
        response = offset + visibility * np.array([0, 1, signal, 0.5])
        probabilities = np.column_stack([response, 1 - response]).ravel() / 4
        tally = rng.multinomial(N // BLOCKS, probabilities).reshape(4, 2)
        ledger.append(tally.tolist())
        feedback = int(tally[:, 0].sum()) % 5
    return ledger


def run():
    rng = np.random.Generator(np.random.PCG64(751129))
    signals = native_probabilities()
    banks = []
    for direction, probabilities in signals.items():
        for seam, p in enumerate(probabilities):
            banks.append(
                dict(
                    axis=direction,
                    seam=seam,
                    blind=False,
                    batches=acquire_bank(rng, p, seam, False),
                )
            )
    banks.append(
        dict(
            axis="port",
            seam=0,
            blind=True,
            batches=acquire_bank(rng, signals["port"][0], 0, True),
        )
    )
    return {
        "events_per_bank": N,
        "arms": ["dark", "bright", "signal", "half"],
        "outcomes": ["click", "miss"],
        "banks": banks,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = json.dumps(run(), separators=(",", ":")) + "\n"
    receipt = HERE / "records.json"
    if args.write:
        receipt.write_text(result)
    elif receipt.read_text() != result:
        raise SystemExit("Acquisition differs from the stored records")
