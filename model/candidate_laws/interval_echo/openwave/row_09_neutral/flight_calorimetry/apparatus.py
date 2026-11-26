"""Native charge positions, translations, local signed pointers and finite clocks."""

from functools import cache
from itertools import product

import numpy as np
from scipy.linalg import eigh, expm

from model.experiments.interval_echo_angular.charged_orientation.live_trimer.measure import (
    LAW,
    blocks,
)
from model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.cyclic_charge.neutral_ladder.release import (
    flight,
)
from model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.cyclic_charge.neutral_ladder.release.detectors import (
    SENSORS,
)

RING = flight.layout()[0]
EDGES = list(zip(RING, RING[1:] + RING[:1]))
DT = 2 * np.pi / 128
APPOINTMENT = 256


def loads(word, edge):
    vacancy, excess = word
    return tuple(0 if v == vacancy else 2 if v == excess else 1 for v in edge)


@cache
def sector(kind, contrast):
    if kind == "pair":
        words = [(a, b) for a in RING for b in RING if a != b]
    elif kind == "vacancy":
        words = [(v, -1) for v in RING]
    elif kind == "excess":
        words = [(-1, v) for v in RING]
    else:
        assert kind == "host"
        words = [(-1, -1)]
    ix = {w: i for i, w in enumerate(words)}
    native = blocks(contrast)["blank"] / 8
    swap = np.eye(9)[[LAW.SRC.index(w[::-1]) for w in LAW.SRC]]
    np.testing.assert_allclose(expm(-1j * np.pi * native), swap, atol=2e-13)
    h = np.zeros((len(words), len(words)))
    costs = []
    shift = []
    for j, w in enumerate(words):
        cost = 0
        for a, b in EDGES:
            word = loads(w, (a, b))
            c = LAW.SRC.index(word)
            cost += contrast * (word[0] - word[1]) ** 2
            h[j, j] += native[c, c]
            if word[0] != word[1]:
                swapped = tuple(b if v == a else a if v == b else v for v in w)
                h[ix[swapped], j] += native[LAW.SRC.index(word[::-1]), c]
        costs.append(cost)
        rolled = w
        # Swaps from last to first implement the cyclic permutation i -> i+1.
        for a, b in reversed(EDGES[:-1]):
            rolled = tuple(b if v == a else a if v == b else v for v in rolled)
        assert rolled == tuple(
            RING[(RING.index(v) + 1) % len(RING)] if v in RING else -1 for v in w
        )
        shift.append(ix[rolled])
    np.testing.assert_array_equal(h, h[np.ix_(shift, shift)])
    np.testing.assert_array_equal(h, h.T)
    return tuple(words), h, np.array(costs), np.array(shift)


@cache
def frame(kind, contrast, zero):
    words, h, _, shift = sector(kind, contrast)
    if zero:
        columns = []
        seen = set()
        for j in range(len(words)):
            if j in seen:
                continue
            orbit = [j]
            while shift[orbit[-1]] != j:
                orbit.append(shift[orbit[-1]])
            seen.update(orbit)
            column = np.zeros(len(words))
            column[orbit] = 1 / np.sqrt(len(orbit))
            columns.append(column)
        basis = np.column_stack(columns)
    else:
        basis = np.eye(len(words))
    block = basis.T @ h @ basis
    np.testing.assert_allclose(h @ basis, basis @ block, atol=3e-14)
    eigenvalues, vectors = eigh(block, driver="evd")
    return eigenvalues, basis @ vectors


def launch(kind, sign):
    pair = (160, 167) if sign == -1 else (167, 160)
    return {"pair": pair, "vacancy": (160, -1), "excess": (-1, 160), "host": (-1, -1)}[
        kind
    ]


@cache
def clock(kind, contrast, zero, sign, size):
    words, *_ = sector(kind, contrast)
    ev, v = frame(kind, contrast, zero)
    initial = v[words.index(launch(kind, sign)), :].conj()
    delay = np.arange(size) * DT
    fringes = np.exp(-1j * ev[:, None] * delay[None, :])
    # The clock changes the state before the signed meter. Its full complex
    # amplitudes enter each exit, including unresolved-level interference.
    before_read = np.exp(-1j * ev * APPOINTMENT) * initial
    wave = v @ (before_read[:, None] * np.fft.ifft(fringes, axis=1))
    probabilities = abs(wave) ** 2
    accepted = float(np.vdot(initial, initial).real)
    np.testing.assert_allclose(probabilities.sum(), accepted, atol=8e-13)
    readings = {}
    for name, edge in SENSORS.items():
        difference = np.array([np.subtract(*loads(w, edge)) for w in words])
        readings[name] = np.column_stack(
            [probabilities[difference == j].sum(axis=0) for j in range(-2, 3)]
        )
    return accepted, readings


def local_meter():
    echo = LAW.echo()
    for a, b in product(range(3), repeat=2):
        col = LAW.SRC.index((a, b))
        raw = echo[:, col].reshape(11, 9)
        labels = np.array([False, *[d == a - b for _, d in LAW.LABELS]])
        assert abs(np.sum(abs(raw[labels]) ** 2) - 1) < 2e-14
        labels = np.array([False, *[abs(d) == 2 for _, d in LAW.LABELS]])
        # Reading magnitude without orientation, then unwriting, implements
        # the strict signed-pair projector, including on coherent orientations.
        measured = (raw * labels[:, None]).ravel()
        np.testing.assert_allclose(
            echo @ measured, np.eye(99)[col] * int(abs(a - b) == 2), atol=2e-14
        )


@cache
def unfiltered(kind, sign, contrast=8):
    words, *_ = sector(kind, contrast)
    ev, v = frame(kind, contrast, False)
    p = (
        abs(
            v
            @ (
                np.exp(-1j * APPOINTMENT * ev)
                * v[words.index(launch(kind, sign))].conj()
            )
        )
        ** 2
    )
    return words, p
