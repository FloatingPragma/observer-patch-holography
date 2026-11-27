"""Acquire all frozen exchange routes and store actual simulated read cohorts."""

import argparse
import gzip
import json
from pathlib import Path

import numpy as np

from ..experiment import LAW
from ..measure import serialize
from .pump import (
    PAIRS,
    SPOKES,
    aperture,
    clocks,
    index_pair,
    route_certificate,
    run_case,
    source_geometry,
    stocks,
)


def fixed_guide():
    # Whole periods of the frozen pair and centre terms are identities.
    _, pairs, centre = source_geometry(2)
    static = (1 - pairs) @ aperture(0) + 16 * centre
    np.testing.assert_allclose(np.exp(-2j * np.pi * static), 1, atol=1e-12)
    settings = {}
    for step in (1 / 16, 1 / 32):
        energies, levels = clocks(2, 16, 0, step)
        initial = np.eye(6)[:, index_pair(PAIRS[0])]
        n = 65536
        binid = round(energies[3] * n * step / (2 * np.pi))
        q = (binid / n - energies * step / (2 * np.pi) + 0.5) % 1 - 0.5
        rec = np.exp(1j * np.pi * (n - 1) * q) * np.sinc(n * q) / np.sinc(q)
        psi = levels @ (rec * (levels.T @ initial))
        prepared = np.column_stack([initial, psi / np.linalg.norm(psi)])
        data = []
        for path in ((1, 3, 2, 1, 3, 2), (1, 3, 2, 2, 3, 1)):
            at = prepared.copy()
            frames = [(abs(at.T) ** 2).tolist()]
            for b in path:
                ix = []
                for w in stocks(2):
                    s = list(w)
                    s[0], s[b] = s[b], s[0]
                    ix.append(stocks(2).index(tuple(s)))
                at = at[ix]
                frames.append((abs(at.T) ** 2).tolist())
            data.append((at, frames))
        d = np.linalg.norm(data[0][0] - data[1][0], axis=0) ** 2 / 4
        settings[str(step)] = dict(
            dark=d.tolist(), forward_words=data[0][1], ref_words=data[1][1]
        )
    for t in (0, 0.343, 1.72):
        for n in (0, 1):
            _, pairs, centre = source_geometry(n)
            np.testing.assert_allclose(pairs @ aperture(t), 0, atol=1e-14)
    return settings


def run():
    turns = (4608, 18432, 73728)
    early = [
        run_case(g, dt, n, turns)
        for g in (8, -8)
        for dt in (1 / 8, 1 / 16)
        for n in (1536, 3072)
    ]
    held = [
        run_case(g, dt, n, (196608, 393216))
        for g in (16, -16)
        for dt, n in ((1 / 16, 49152), (1 / 32, 98304))
    ]
    guarded = [run_case(g, 1 / 16, 49152, (393216,), guarded=True) for g in (16, -16)]
    parked = [
        run_case(g, 1 / 16, 49152, (393216,), guarded="parked") for g in (16, -16)
    ]
    return dict(
        native=route_certificate(),
        fixed_angle=fixed_guide(),
        word_order=stocks(2),
        parked=parked,
        early=early,
        held=held,
        guarded=guarded,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    dest = parser.parse_args().output
    packet = serialize(run())
    record = Path(__file__).with_name("records.json.gz")
    if dest:
        dest.write_bytes(
            gzip.compress(
                (json.dumps(packet, separators=(",", ":")) + "\n").encode(), mtime=0
            )
        )
    else:
        previous = json.loads(gzip.decompress(record.read_bytes()))
        if previous != packet:
            raise ValueError(
                "The native exchange acquisition differs from its deposited record"
            )
        print("Replayed all exchange sources, apertures and detector outcomes.")
