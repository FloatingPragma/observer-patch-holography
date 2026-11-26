"""Compare the ordinary KG group forecast with archived, independently ranged wall packets."""

import json

import numpy as np

from ..packet_speed import measure as packet


def compare():
    meter = json.loads((packet.HERE.parent / "coupling_meter.json").read_text())[
        "results"
    ]
    h, native = (meter[field]["value"] for field in ("h", "J"))
    grid = 2 * np.pi * np.fft.fftfreq(4096)
    summary = {}
    with np.load(packet.HERE / "counts.npz", allow_pickle=False) as raw:
        for n in (29, 49):
            x = np.arange(n + 1)
            envelope = np.exp(-((x - n / 2) ** 2) / (4 * min(4.0, (n + 1) / 7) ** 2))
            for contrast in (8, 4):
                j = native * contrast / 8
                mass, c2 = 2 * abs(j - h), 4 * j * h
                continuum = c2 * grid / np.sqrt(mass**2 + c2 * grid**2)
                lattice = (
                    c2 * np.sin(grid) / np.sqrt(mass**2 + 2 * c2 * (1 - np.cos(grid)))
                )
                rows = []
                for idx in range(len(packet.MOMENTA)):
                    key = f"n{n}_p{idx}"
                    pos, neg = (packet.momentum(axis)[0] for axis in raw[key + "_wave"])
                    momentum = float((pos - neg) / 2)
                    power = (
                        abs(np.fft.fft(envelope * np.exp(1j * x * momentum), len(grid)))
                        ** 2
                    )
                    timed = f"{key}_j{contrast}"
                    _, observed = packet.slope(
                        *packet.clock(raw[timed + "_clock"]), raw[timed + "_wall"]
                    )
                    con = float(np.average(continuum, weights=power))
                    lat = float(np.average(lattice, weights=power))
                    rows.append(
                        dict(
                            k=momentum,
                            velocity=observed,
                            continuum=con,
                            lattice=lat,
                            continuum_miss=max(
                                0.0, con - observed[1], observed[0] - con
                            ),
                            lattice_miss=max(0.0, lat - observed[1], observed[0] - lat),
                        )
                    )

                def count(records, law):
                    return sum(r[law + "_miss"] == 0 for r in records)

                low = [r for r in rows if r["k"] <= 0.50]
                high = [r for r in rows if r["k"] > 0.50]
                first = rows[1]["velocity"][0] - rows[0]["velocity"][1]
                last = rows[-2]["velocity"][0] - rows[-1]["velocity"][1]
                summary[f"{n}/{contrast}"] = dict(
                    rows=rows,
                    first_speed_increase_at_least=first,
                    final_speed_decrease_at_least=last,
                    low_continuum=count(low, "continuum"),
                    high_continuum=count(high, "continuum"),
                    low_lattice=count(low, "lattice"),
                    high_lattice=count(high, "lattice"),
                )
    return summary


if __name__ == "__main__":
    print(json.dumps(compare(), indent=2, allow_nan=False))
