"""Read all independent whole-region momentum banks on the two retained controls."""

import json
from pathlib import Path

from mpmath import iv
from mpmath.libmp import to_int, round_floor, round_ceiling


HERE = Path(__file__).parent
SOURCE = HERE.parents[3]


def evidence():
    iv.dps = 45
    base = json.loads((HERE / "results.json").read_text())
    source = json.loads((SOURCE / "full_algebra/certificate.json").read_text())[
        "signed_transfers"
    ]
    angles = json.loads(
        (SOURCE / "signed_transfers.json").read_text(), parse_float=str
    )["tests"]
    a, c = [
        iv.mpf(next(y["angle_per_unit"] for y in angles if y["source"] == name))
        for name in ("ordinary", "collar")
    ]
    material = SOURCE / "directional_stress/translation"
    certified = json.loads((material / "certificate.json").read_text())
    banks = json.loads((material / "counts.json").read_text(), parse_float=str)
    n = banks["simulated_trials_per_setting"]
    assert n == 2**76

    def bound(text):
        return iv.mpf(text.strip("[]").split(","))

    def outward(x):
        lower, upper = (x * 10**12)._mpi_
        return [to_int(lower, round_floor), to_int(upper, round_ceiling)]

    energy_ratios = base["tracers"]["-1"]["twelve_cell_pairing_over_ordinary"]
    ratios = {
        "fixed_full_entropy": bound(source["pairing_entropy"])
        / bound(source["ordinary_entropy"]),
        "fixed_vacancy_energy": iv.mpf([str(t) for t in energy_ratios]) / 10**8,
    }
    comparisons = {}
    for name, r in ratios.items():
        comparison = {}
        for q in ("-1", "1"):
            rows = {}
            for channel in ("boundary", "marked_force", "free_rate"):
                b = certified[q][channel]
                pop = bound(b["pairing"]) - r * bound(b["ordinary"])
                m = banks["signs"][q][channel]
                scale = iv.mpf(m["range"])
                means = {
                    arm: scale * (2 * iv.mpf(h) / n - 1)
                    for arm, h in m["plus_counts"].items()
                }
                raw = (means["pairing"] - means["equilibrium"]) + (
                    means["plus"] - means["minus"]
                ) * c / 2
                ordinary = (means["ordinary"] - means["equilibrium"]) + (
                    means["plus"] - means["minus"]
                ) * a / 2
                estimate = raw - r * ordinary
                # All 30 independent binary acquisitions on a single event.
                noise = scale * iv.sqrt(2 * iv.log(2 * 30 / iv.mpf(".01")) / n)
                width = (2 + abs(c - r * a)) * noise
                count = estimate + iv.mpf([-1, 1]) * width
                rows[channel] = {
                    "population": outward(pop),
                    "joint_99pct_counts": outward(count),
                }
            comparison[q] = rows
        comparisons[name] = comparison
    for meter in ("boundary", "free_rate"):
        for estimate in ("population", "joint_99pct_counts"):
            x = comparisons["fixed_full_entropy"]["-1"][meter][estimate]
            y = comparisons["fixed_full_entropy"]["1"][meter][estimate]
            assert x[1] < 0 < y[0] if meter == "boundary" else y[1] < 0 < x[0]
    return {
        "units": "outward intervals in trillionths of a native momentum-rate unit",
        "controls": comparisons,
    }


if __name__ == "__main__":
    print(json.dumps(evidence(), indent=2))
