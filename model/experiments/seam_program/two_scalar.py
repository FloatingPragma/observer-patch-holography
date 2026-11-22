"""Two protected mobile packets on the native scalar chart, with a matched contact control."""

import argparse
from hashlib import sha256
from itertools import combinations_with_replacement
import json
from pathlib import Path

import numpy as np

from model.candidate_laws.contrast_memory_spherical_v1.sectors import (
    PairSector,
    SeamEvents,
)
from .clocked_apparatus import clock_laws, patches
from .compiler import Compiler

HERE = Path(__file__).resolve().parent
LAUNCHES = ((0, 1), (0, 7))
ARMS = ("native", "quarter_contact_cancelled", "free_bosons")


def detectors(comp, sector):
    nstations = int(np.prod(comp.shape[:3]))
    station = np.full(sector.n, -1)
    station[comp.slots] = np.repeat(np.arange(nstations), len(comp.slots) // nstations)
    labels = list(combinations_with_replacement(range(nstations), 2))
    lookup = np.full((nstations + 1, nstations + 1), len(labels), dtype=int)
    for index, (a, b) in enumerate(labels):
        lookup[a, b] = lookup[b, a] = index
    return (
        labels,
        lookup[station[sector.aa], station[sector.bb]],
        lookup[station[comp.slots][:, None], station[comp.slots][None, :]],
    )


def initial(comp, sector):
    block = np.array(comp.slots).reshape(np.prod(comp.shape[:3]), -1)
    pairs = np.zeros((sector.P, 2 * len(LAUNCHES)), dtype=complex)
    singles = np.zeros((len(comp.slots), 2 * len(LAUNCHES)), dtype=complex)
    width = block.shape[1]
    for j, (a, b) in enumerate(LAUNCHES):
        for first in block[a]:
            for second in block[b]:
                pairs[sector.index(first, second), 2 * j : 2 * j + 2] = 1 / width
        singles[a * width : (a + 1) * width, 2 * j] = 1 / np.sqrt(width)
        singles[b * width : (b + 1) * width, 2 * j + 1] = 1 / np.sqrt(width)
    return pairs, singles


def joint(pairs, singles, pair_bins, boson_bins, nbins):
    rows = []
    for j in range(len(LAUNCHES)):
        a, b = singles[:, 2 * j : 2 * j + 2].T
        product = np.outer(a, b)
        symmetric = (product + product.T) / np.sqrt(2)
        distributions = [
            np.bincount(pair_bins, weights=abs(psi) ** 2, minlength=nbins)
            for psi in pairs[:, 2 * j : 2 * j + 2].T
        ]
        distributions.append(
            np.bincount(
                boson_bins.ravel(),
                weights=(abs(symmetric) ** 2).ravel(),
                minlength=nbins,
            )
        )
        rows.append(distributions)
    rows = np.array(rows)
    np.testing.assert_allclose(rows.sum(axis=-1), 1, rtol=0, atol=2e-11)
    return rows


def run():
    fixture = json.loads((HERE / "clocked_protocol.json").read_text())["fixture"]
    host, data, _, front, bank = patches(fixture["level"])
    weights, _ = clock_laws(host, front, bank, fixture)
    comp = Compiler(data, fixture["marker"], fixture["shape"])
    word = comp.build("nominal")
    step = comp.columns(word)[comp.slots]
    sector = PairSector(data)
    events = SeamEvents(sector, 8)
    labels, pair_bins, boson_bins = detectors(comp, sector)
    state, singles = initial(comp, sector)
    control_phase = np.exp(-1j * np.pi / 4)
    periods = []
    for period in range(fixture["alarm_after_periods"] + 1):
        periods.append(joint(state, singles, pair_bins, boson_bins, len(labels) + 1))
        if period == fixture["alarm_after_periods"]:
            break
        for edge, q in word:
            events.apply(state, edge, q * np.pi / 4)
            if q == 1:
                both = events.tables[edge][2]
                state[both, 1::2] *= control_phase
        singles = step @ singles
    periods = np.array(periods)
    observations = []
    for r in (1, 2, 3):
        events_at_stop = np.einsum("t,tjab->jab", weights[r], periods)
        timeout = np.full((*events_at_stop.shape[:2], 1), 1 - weights[r].sum())
        stopped = np.concatenate((events_at_stop, timeout), axis=-1)
        for j, launch in enumerate(LAUNCHES):
            native, changed, free = stopped[j]
            observations.append(
                dict(
                    launch=launch,
                    clock=r,
                    laws=stopped[j].round(12).tolist(),
                    tv_contact=float(sum(abs(native - changed)) / 2),
                    tv_free=float(sum(abs(native - free)) / 2),
                )
            )
    return dict(
        protocol_sha256=sha256(
            (HERE / "two_scalar_protocol.json").read_bytes()
        ).hexdigest(),
        source="conditional normalized two mobile packets with a stationary marker",
        labels=labels + ["outside chart"],
        arms=ARMS,
        launches=LAUNCHES,
        period_probabilities=periods.round(12).tolist(),
        stopped=observations,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.parse_args().output.write_text(json.dumps(run(), indent=2, allow_nan=False) + "\n")
