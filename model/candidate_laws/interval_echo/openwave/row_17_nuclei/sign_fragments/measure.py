"""Spectral counts on joined signed stocks and independently occupied copies."""

import argparse
import json
from collections import defaultdict
from functools import cache
from itertools import product
from math import sqrt, log
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.row_14_profile.measure import (
    EDGES,
    native_geometry,
)

HERE = Path(__file__).resolve().parent
SHOTS = 2**20
STOCKS = ((1, 1), (2, 1), (2, 2), (4, 3))
STRENGTHS = (8, 4, 0.25)
PHASES = ("sharp", "even", "odd")
SIDES = ("joint", "excess", "vacancy")
SRC = candidate.native_law.SRC


@cache
def contact(k):
    exchange, contrast = candidate.native_law.local_term_parts()
    h = exchange + k * contrast
    np.testing.assert_allclose(h[9:, :9], 0, atol=3e-14)
    np.testing.assert_allclose(h[:9, :9].imag, 0, atol=3e-14)
    return h[:9, :9].real


def separate_indices():
    return [
        SRC.index(tuple(max(1, t) for t in w)) * 9
        + SRC.index(tuple(min(1, t) for t in w))
        for w in SRC
    ]


def local_certificate(k):
    h = contact(k)
    ids = separate_indices()
    separate = np.kron(h, np.eye(9)) + np.kron(np.eye(9), h)
    difference = h - separate[np.ix_(ids, ids)]
    pair = [SRC.index(w) for w in ((2, 0), (0, 2))]
    block = (2 * k - 0.5) * np.eye(2) - 0.5 * np.array([[0, 1], [1, 0]])
    expected = np.zeros((9, 9))
    expected[np.ix_(pair, pair)] = block
    np.testing.assert_allclose(difference, expected, atol=2e-13)
    return np.linalg.eigvalsh(block).tolist()


@cache
def meter(k):
    """A 256-delay phase meter resolves each quarter-unit local eigenspace."""
    values, vectors = np.linalg.eigh(contact(k))
    labels = np.rint(values * 4).astype(int)
    np.testing.assert_allclose(values * 4, labels, atol=2e-13)
    tick = 8 * np.pi / 256
    output = np.fft.ifft(np.exp(-1j * np.outer(np.arange(256), values) * tick), axis=0)
    np.testing.assert_allclose(abs(output[labels, np.arange(9)]), 1, atol=1e-13)
    return {
        int(q): vectors[:, labels == q] @ vectors[:, labels == q].T
        for q in np.unique(labels)
    }


@cache
def cohort(positive, negative, phase):
    occupied = [1] * 12
    occupied[0], occupied[1] = 2, 0
    occupied[2 : positive + 1] = [2] * (positive - 1)
    occupied[positive + 1 : positive + negative] = [0] * (negative - 1)
    original = tuple(occupied)
    if phase == "sharp":
        return [original], np.ones(1)
    swapped = (0, 2, *original[2:])
    ids = [SRC.index(w[:2]) for w in (original, swapped)]
    pre = expm(-1j * np.pi / 2 * contact(8))[ids, ids[0]]
    dwell = (3 if phase == "even" else 1) * np.pi / 64
    prepared = np.exp(-1j * 8 * np.array([4, 0]) * dwell) * pre
    ratio = 1 if phase == "even" else -1
    np.testing.assert_allclose(prepared / prepared[0], [1, ratio], atol=2e-13)
    return [original, swapped], prepared


@cache
def parity_effect():
    """Undo the native even-pair preparation and sharply read (2,0)."""
    ids = [SRC.index((2, 0)), SRC.index((0, 2))]
    # The blank exchange spectrum on this sector is 32,33. The forward
    # pi/2 contact is inverted by 3pi/2, up to a common phase.
    contact_inverse = expm(-1j * 3 * np.pi / 2 * contact(8))[np.ix_(ids, ids)]
    private_inverse = np.diag(np.exp(-1j * 8 * np.array([4, 0]) * np.pi / 64))
    pulse = contact_inverse @ private_inverse
    effect = np.outer(pulse[0].conj(), pulse[0])
    ideal = np.full((2, 2), 0.5)
    np.testing.assert_allclose(effect, ideal, atol=5e-13, rtol=0)
    return ideal


def local_density(words, amplitudes, side, sites):
    traced = defaultdict(lambda: np.zeros(3 ** len(sites), complex))
    for word, amplitude in zip(words, amplitudes, strict=True):
        if side == "joint":
            state, partner = word, ()
        else:
            positive = tuple(max(t, 1) for t in word)
            negative = tuple(min(t, 1) for t in word)
            state, partner = (
                (positive, negative) if side == "excess" else (negative, positive)
            )
        key = tuple(state[i] for i in range(12) if i not in sites) + partner
        symbol = tuple(state[i] for i in sites)
        idx = symbol[0] if len(sites) == 1 else SRC.index(symbol)
        traced[key][idx] += amplitude
    return sum(np.outer(v, v.conj()) for v in traced.values())


def sample(rng, probabilities):
    p = np.array(probabilities)
    if p.min() < -1e-12 or abs(p.sum() - 1) > 1e-11:
        raise ValueError("incomplete local energy instrument")
    p = np.clip(p, 0, None)
    return rng.multinomial(SHOTS, p / p.sum()).tolist()


def energy_read(rng, k, words, amplitudes, side, sites):
    rho = local_density(words, amplitudes, side, sites)
    if len(sites) == 1:
        labels = [int(4 * k * (n * n - 1)) for n in range(3)]
        probabilities = rho.diagonal().real
    else:
        projectors = meter(k)
        labels = list(projectors)
        probabilities = [np.trace(p @ rho).real for p in projectors.values()]
    return [labels, sample(rng, probabilities)]


def acquire():
    rng = np.random.default_rng(170126)
    heat = []
    mask, comparison = native_geometry()
    sites = [tuple(edge) for edge in EDGES] + [(i,) for i in np.flatnonzero(mask)]
    for (plus, minus), k, phase in product(STOCKS, STRENGTHS, PHASES):
        words, amps = cohort(plus, minus, phase)
        banks = {
            side: [energy_read(rng, k, words, amps, side, term) for term in sites]
            for side in SIDES
        }
        occupancy = sample(rng, abs(amps) ** 2)
        ia, ib = SRC.index((2, 0)), SRC.index((0, 2))
        rho = local_density(words, amps, "joint", (0, 1))
        twice_coherence = 2 * rho[ia, ib].real
        sharp_return = np.trace(parity_effect() @ rho[np.ix_([ia, ib], [ia, ib])]).real
        even = (1 + twice_coherence) / 2
        np.testing.assert_allclose(sharp_return, even, atol=5e-13, rtol=0)
        parity = sample(rng, [even, (1 - twice_coherence) / 2])
        heat.append(
            dict(
                plus=plus,
                minus=minus,
                k=k,
                phase=phase,
                energy=banks,
                sharp_word_counts=occupancy,
                parity=parity,
            )
        )
    return dict(
        shots=SHOTS,
        seed=170126,
        second_owner=comparison,
        read_ports=[[int(i) for i in row] for row in sites],
        banks=heat,
    )


def decode(data):
    mask, _ = native_geometry()
    required_ports = [list(map(int, row)) for row in EDGES] + [
        [int(i)] for i in np.flatnonzero(mask)
    ]
    enrolled = {
        (p, m, k, phase) for (p, m), k, phase in product(STOCKS, STRENGTHS, PHASES)
    }
    seen = {(r["plus"], r["minus"], r["k"], r["phase"]) for r in data["banks"]}
    if (
        data["shots"] != SHOTS
        or data["read_ports"] != required_ports
        or len(data["banks"]) != len(enrolled)
        or seen != enrolled
    ):
        raise ValueError("wrong exposure or number of sources")
    # For 36 energy and 36 independent predictor intervals, simultaneous 99%.
    logtail = log(2 * 72 / 0.01)
    out = []
    for bank in data["banks"]:

        def check(counts):
            if (
                not counts
                or any(type(c) is not int for c in counts)
                or sum(counts) != SHOTS
                or min(counts) < 0
            ):
                raise ValueError("malformed bank counts")

        k, phase = bank["k"], bank["phase"]
        means, ranges_sq = [], 0
        for side in SIDES:
            if len(bank["energy"][side]) != 39:
                raise ValueError("missing local energy appointment")
            total = 0
            for term, (labels, counts) in zip(
                required_ports, bank["energy"][side], strict=True
            ):
                expected = (
                    list(meter(k))
                    if len(term) == 2
                    else [int(4 * k * (n * n - 1)) for n in range(3)]
                )
                if labels != expected or len(counts) != len(labels):
                    raise ValueError("invalid terminal energy labels")
                check(counts)
                total += np.dot(labels, counts) / (4 * SHOTS)
                ranges_sq += ((max(labels) - min(labels)) / 4) ** 2
            means.append(float(total))
        words, _ = cohort(bank["plus"], bank["minus"], phase)
        adjacent = [sum(abs(w[a] - w[b]) == 2 for a, b in EDGES) for w in words]
        if len(bank["sharp_word_counts"]) != len(words) or len(bank["parity"]) != 2:
            raise ValueError("invalid parity or stock record")
        check(bank["sharp_word_counts"])
        check(bank["parity"])
        contact_count = np.dot(adjacent, bank["sharp_word_counts"]) / SHOTS
        parity = (bank["parity"][0] - bank["parity"][1]) / SHOTS
        # On other occupied opposite seams the expected swap is zero:
        # their swap leads outside the two enrolled source words.
        pred = (2 * k - 0.5) * contact_count - 0.5 * parity
        pred_range_sq = ((2 * k - 0.5) * (max(adjacent) - min(adjacent))) ** 2 + 1
        out.append(
            dict(
                stock=[bank["plus"], bank["minus"]],
                phase=phase,
                k=k,
                measured=means[0] - means[1] - means[2],
                means=means,
                bound=sqrt(logtail * ranges_sq / (2 * SHOTS)),
                forecast=pred,
                forecast_bound=sqrt(logtail * pred_range_sq / (2 * SHOTS)),
            )
        )
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquire", type=Path)
    a = parser.parse_args()
    if a.acquire:
        a.acquire.write_text(json.dumps(acquire(), separators=(",", ":")) + "\n")
    else:
        print(
            json.dumps(decode(json.loads((HERE / "counts.json").read_text())), indent=2)
        )
