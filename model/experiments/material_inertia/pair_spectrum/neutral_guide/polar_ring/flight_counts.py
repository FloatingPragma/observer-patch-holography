"""Read the full occupancy-word banks for all five deposited neutral polar rings."""

import json
from pathlib import Path

import numpy as np


def main():
    data = Path(__file__).resolve().parent
    packet = json.loads((data / "records.json").read_bytes())
    shots = packet["shots_per_bank"]
    assert packet["flight_sources"][:2] == ["local single upper", "all lower"]
    assert shots == 2**22 and len(packet["boards"]) == 6
    assert sum(bool(b["stock"]) for b in packet["boards"]) == 5
    meter = packet["meter"]
    for key in ("01", "10"):
        assert (meter[key]["signed_change"], meter[key]["squared_change"]) == (0, 0)
    assert (meter["00"]["signed_change"], meter["02"]["signed_change"]) == (-1, 1)
    panels = [p for p in packet["panels"] if p["kappa"] == 8]
    assert len(panels) == 5
    with np.load(data / "word_counts.npz", allow_pickle=False) as words:
        for panel in panels:
            n = panel["dimer_count"]
            counts = words[panel["prefix"] + "_t64"]
            assert counts.shape == (4, 2**n) and counts.min() >= 0
            assert np.all(counts.sum(axis=1) == shots)
            # Dimers 0 and 1 belong to the launch owner. A single raised dimer
            # has exactly one set bit. Every other word remains in the denominator.
            away = counts[:, [1 << site for site in range(2, n)]].sum(axis=1)
            saved = next(s for s in panel["flight_summary"] if s["time"] == 64)
            assert list(away) == saved["unique_upper_other_owner"]
            gap = panel["clocks"][0]
            assert gap["k"] == 0 and abs(gap["energy"] - panel["rest_gap"]) < 0.0002
            assert panel["isolation_lower_above"] > 0
            print(
                f"{panel['cells']} cells, {n} dimers: local flip {away[0] / shots:.5%}, "
                f"sham {away[1] / shots:.5%}, rest clock {gap['energy']:.6f}."
            )


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run without -O to enable the scientific checks.")
    main()
