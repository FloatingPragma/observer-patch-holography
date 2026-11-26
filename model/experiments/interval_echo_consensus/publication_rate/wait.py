"""Read one continuously flying flag at its endpoint and reconstruct the bright ports."""

import argparse
import gzip
import json
from itertools import combinations
from pathlib import Path

import numpy as np

from . import apparatus as a
from . import measure as m

DELAYS = (0, 1 / 8, 1 / 4, 1 / 2, 1, 2, 3, 6)


def acquire():
    blocks = a.triangle_blocks(a.FINE)
    banks = []
    for appointment, delay in enumerate(DELAYS):
        flag = a.flag_probability(np.pi / 3 * delay, 1)
        hand = a.dimer(np.pi / 16 * delay)[:, 0]
        phases = [
            m.counts(abs(state) ** 2, 10, appointment, quad, seed=177032)
            for quad, state in enumerate((hand, a.dimer(np.pi / 2).conj().T @ hand))
        ]
        readings = []
        for k, (i, j) in enumerate(combinations(range(3), 2)):
            for phase in (0, 1):
                u = a.pair_word(i, j, phase * np.pi / 2)
                src = u[:, i]
                bank = a.appointment(
                    {("blank",): np.outer(src, src.conj())}, blocks, flag, 0
                )
                obs = m.read_bank(bank, u, [11, appointment, k, phase], seed=177032)
                readings.append(dict(pair=[i, j], phase=phase, **obs))
        banks.append(dict(delay=delay, hand=phases, readings=readings))
    return banks


def forecast(data, commission):
    calibrated = m.report(commission)
    one = calibrated["clock_phases"]["1"]
    angle = 2 * np.arcsin(np.sqrt(calibrated["flags"]["1.0"][0]))
    rate = calibrated["rates_per_radian"][0]
    result = []
    for row in data:
        cos, sin = (2 * x[0] / m.N - 1 for x in row["hand"])
        hand = np.arctan2(sin, cos)
        # Zero is a known rest preparation. A slightly negative statistical
        # estimate there is an origin read, not a negative flight.
        elapsed = max(0, hand) / one
        exp_forecast = (1 + np.exp(-rate * max(0, hand))) / 2
        flight_forecast = 1 - 0.5 * np.sin(angle * elapsed / 2) ** 2
        bright = [
            float(np.sum(r["terminal_ports"], axis=0)[r["pair"][0]] / m.N)
            for r in row["readings"]
        ]
        result.append(
            dict(
                delay=row["delay"],
                hand_radians=float(hand),
                continuous_semigroup=float(exp_forecast),
                single_flight=float(flight_forecast),
                smallest_bright=min(bright),
                largest_bright=max(bright),
            )
        )
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    path = parser.parse_args().output
    data = acquire()
    if path:
        path.write_bytes(gzip.compress(m.plain(data), mtime=0))
    elif gzip.decompress((m.HERE / "waits.json.gz").read_bytes()) != m.plain(data):
        raise ArithmeticError("Uninterrupted count replay differs")
    primary = json.loads(gzip.decompress((m.HERE / "records.json.gz").read_bytes()))
    print(json.dumps(forecast(data, primary), indent=2))
