"""Decode the complete clock/word banks without using the source Hamiltonian."""

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def decode(data, banks):
    out = []
    for panel in data["panels"]:
        n, mode = panel["n"], panel["mode"]
        base = f"{mode}_{n}"
        u = (np.arange(1 << n)[:, None] >> np.arange(n)) & 1
        rises = (np.diff(np.pad(u, ((0, 0), (1, 1))), axis=1) > 0).sum(axis=1)
        peaks = {}
        for src in panel["input_words"]:
            clock = banks[base + "_clock_" + src]
            if clock.sum() != data["shots"]:
                raise ValueError("Missing clock outcomes")
            b = int(np.argmax(clock))
            frequency = (
                (b if b < len(clock) // 2 else b - len(clock))
                * data["energy_period"]
                / len(clock)
            )
            peaks[src] = frequency
        flights = {}
        for length in (1, 2):
            expected = panel["input_words"]["one" if length == 1 else "two"]
            matching = (rises == 1) & (u.sum(axis=1) == length)
            moved = matching.copy()
            moved[expected] = False
            rows = []
            for t in (0, 2, 16, 32, 128):
                x = banks[f"{base}_flight_{length}_{t}"]
                if x.sum() != data["shots"] or len(x) != 3 + (1 << n):
                    raise ValueError("Missing signed wall readings")
                rows.append([int(x[3:][matching].sum()), int(x[3:][moved].sum())])
            flights[str(length)] = rows
        outcome = {}
        for length in (1, 2):
            for control in ("near", "parked", "far"):
                counts = banks[f"{base}_convert_{length}_{control}"]
                if counts.sum() != data["shots"]:
                    raise ValueError("Missing conversion refusal")
                outcome[f"{length}_{control}"] = int(counts[3:].sum())
        out.append(
            dict(
                n=n,
                mode=mode,
                dominant_ordinary_bins=peaks,
                flights=flights,
                conversion_both=outcome,
            )
        )
    return out


def read_colour(data):
    out = {}
    for key, panel in data["records"].items():
        out[key] = {}
        for t, rec in panel.items():
            counts = np.array(rec["counts"])
            if counts.shape != (4, 10) or np.any(
                counts.sum(axis=1) != data["shots_per_setting"]
            ):
                raise ValueError("Missing colour admission or joint refusal")
            value = (counts[:, [0, 4, 8]].sum() - counts[:, :9].sum() / 4) / (
                3 * data["shots_per_setting"]
            )
            if abs(value - rec["inferred_singlet"]) > 1e-14:
                raise ValueError("Inconsistent singlet estimate")
            out[key][t] = float(value)
    return out


if __name__ == "__main__":
    with np.load(HERE / "records.npz") as bank:
        states = decode(json.loads((HERE / "records.json").read_text()), bank)
    print(
        json.dumps(
            dict(
                bands=states,
                colour=read_colour(json.loads((HERE / "colour.json").read_text())),
            ),
            indent=1,
        )
    )
