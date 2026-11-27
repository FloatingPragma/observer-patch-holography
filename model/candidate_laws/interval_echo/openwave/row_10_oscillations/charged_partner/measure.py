"""Detect the energy and total charge of one captured daughter at two stations."""

import argparse
from functools import cache
from itertools import product
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh, expm

from ..measure import (
    BASES,
    LAW,
    OWNERS,
    blocks,
    density,
    enrol_matrix,
    flight,
    geometry,
    receiver,
    source_weights,
)
from ..joint_receivers.measure import sample

HERE = Path(__file__).resolve().parent
SEED, STARTS, BENCH = 26112687, 2**52, 2**22
TIMES = (256, 768)
MODES = ("native", "position", "one_contact", "contrast_four", "pickup_parked")
N_CLOCK = 4096
RANGE = 32


@cache
def station(owner):
    ring = geometry.layout()[0]
    r = receiver(owner)
    i = ring.index(r[0])
    sites = tuple(ring[(i + j) % len(ring)] for j in range(5))
    assert sites[:3] == r and len(set(sites)) == 5
    blanks = {frozenset((a, b)) for a, b, flag in flight.edges(True) if flag is None}
    assert all(frozenset(edge) in blanks for edge in zip(sites, sites[1:]))
    block = blocks()["blank"] / 8
    sector = [LAW.SRC.index(w) for w in [(2, 1), (1, 2)]]
    pulse = expm(-1j * np.pi * block)
    np.testing.assert_allclose(
        abs(pulse[:, sector[0]]) ** 2, np.eye(len(block))[sector[1]], atol=1e-12
    )
    return sites


@cache
def hamiltonian(one_contact):
    """Host-subtracted two-contact energy on the one-vacancy trimer."""
    h = np.zeros((3, 3))
    block = blocks()["blank"] / 8
    host = LAW.SRC.index((1, 1))
    for left in range(1 if one_contact else 2):
        for v in range(3):
            pair = tuple(int(k != v) for k in (left, left + 1))
            col = LAW.SRC.index(pair)
            h[v, v] += block[col, col] - block[host, host]
            if min(pair) == 0:
                h[2 * left + 1 - v, v] += block[LAW.SRC.index(pair[::-1]), col]
    np.testing.assert_array_equal(h, h.T)
    return h


def clock(energies):
    # Geometric series for a rectangular 4096-path delay interferometer.
    distance = energies[:, None] / RANGE - np.arange(N_CLOCK)[None, :] / N_CLOCK
    p = (np.sinc(N_CLOCK * distance) / np.sinc(distance)) ** 2
    np.testing.assert_allclose(p.sum(axis=1), 1, atol=2e-12)
    return p


@cache
def instrument(one_contact):
    e, u = eigh(hamiltonian(one_contact))
    channels = clock(e)
    centers = np.r_[
        np.linalg.eigvalsh(hamiltonian(False)), np.linalg.eigvalsh(hamiltonian(True))
    ]
    bins = np.arange(N_CLOCK) * RANGE / N_CLOCK
    # Periodic clock output, so the zero standard straddles bins 4095 and 0.
    windows = abs((bins[:, None] - centers + RANGE / 2) % RANGE - RANGE / 2) <= 0.016
    assert np.all(windows.sum(axis=1) <= 1)
    selector = np.column_stack([windows, ~windows.any(axis=1)])
    return e, u, channels @ selector


def spectrum(rho, one_contact):
    _, u, r = instrument(one_contact)
    p = (u.T @ rho @ u).diagonal().real @ r
    np.testing.assert_allclose(p.sum(), rho.trace().real, rtol=1e-12, atol=1e-26)
    return p


def standards():
    cal = []
    for one_contact in (False, True):
        for stock in range(4):
            # Stock 3 is the independently filled trimer of zero relative energy.
            p = (
                np.eye(7)[3]
                if stock == 3
                else spectrum(np.diag(np.eye(3)[stock]), one_contact)
            )
            key = [one_contact, stock]
            cal.append(
                dict(key=key, counts=sample(["bench", *key], p, BENCH, seed=SEED))
            )
    return cal


def fields(strength):
    """Remote position coherence survives tracing the source on each accepted branch."""
    result = {o: [] for o in OWNERS}
    for mask in range(4):
        states, native = flight.chart((160, 167, mask), True)
        h = enrol_matrix(states, strength, False)
        if strength == 8:
            np.testing.assert_array_equal(h, native)
        values, vectors = eigh(h, driver="evd")
        field = (np.exp(-1j * np.outer(TIMES, values)) * vectors[0]) @ vectors.T
        index = {s: i for i, s in enumerate(states)}
        for owner in OWNERS:
            *rails, pick, donor = station(owner)
            ids = [index[v, pick, mask] for v in rails]
            ids_pick = [
                i for i, (v, e, m) in enumerate(states) if e == pick and v != donor
            ]
            assert len(ids) == 3 and set(ids) <= set(ids_pick)
            assert all(m == mask for v, e, m in states if v in rails and e == pick)
            result[owner].append(
                (field[:, ids], (abs(field[:, ids_pick]) ** 2).sum(axis=1))
            )
    return result


def acquire():
    standards_data = standards()
    stocks = source_weights()
    banks = []
    for strength in (8, 4):
        arrivals = fields(strength)
        for owner, (index, time), recipe in product(OWNERS, enumerate(TIMES), range(6)):
            weights = stocks[:, recipe]
            columns = np.array([item[0][index] for item in arrivals[owner]])
            rho = density(weights, columns)
            at_pick = float(weights @ [item[1][index] for item in arrivals[owner]])
            pair = float(rho.trace().real)
            supplied = float(weights.sum())
            for mode in MODES:
                if (mode == "contrast_four") != (strength == 4):
                    continue
                measured = np.diag(rho.diagonal()) if mode == "position" else rho
                transferred = mode != "pickup_parked"
                energy = (
                    spectrum(measured, mode == "one_contact")
                    if transferred
                    else np.zeros(7)
                )
                missed = [
                    at_pick - pair if transferred else 0.0,
                    0.0 if transferred else at_pick,
                    supplied - at_pick,
                    1 / 3 - supplied,
                    2 / 3,
                ]
                key = [mode, owner, time, BASES[recipe // 3], recipe % 3]
                banks.append(
                    dict(
                        key=key,
                        counts=sample(key, [*energy, *missed], STARTS, seed=SEED),
                    )
                )
    return dict(
        seed=SEED,
        starts=STARTS,
        calibration_starts=BENCH,
        stations={o: station(o) for o in OWNERS},
        energies={str(q): instrument(q)[0].tolist() for q in (False, True)},
        bins=[
            "three full-trimer windows",
            "three one-contact windows",
            "remaining accepted energy",
            "pickup transferred, vacancy outside trimer",
            "parked pickup failed verification",
            "pickup stock refused",
            "production miss",
            "emitter refused",
        ],
        standards=standards_data,
        banks=banks,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--acquire", type=Path, required=True)
    args = p.parse_args()
    args.acquire.write_text(
        json.dumps(acquire(), separators=(",", ":"), allow_nan=False) + "\n"
    )
