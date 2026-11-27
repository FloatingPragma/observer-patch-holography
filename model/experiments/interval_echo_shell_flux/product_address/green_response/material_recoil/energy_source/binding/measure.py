"""Acquire local energy marks, a parked Ramsey beat and two distal material records."""

import argparse
from functools import cache
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from .. import instrument
from ...common_receiver import measure as guide

HERE = Path(__file__).parent
N = 2**31
WIDTH = np.sqrt(np.log(2 * 64 / 0.005) / (2 * N)) + 2**-32
MATERIALS = ((1, 1), (1, 2), (2, 2))
STOCKS = {"near": [1, 0], "far": [0, 1], "even": [1, 1], "odd": [1, -1]}
DWELLS = (0, 0.13)
SHUTTER = 1 / 32
SCALE = 512


def probability(x):
    return round(float(np.clip(x, 0, 1)) * 2**32) / 2**32


@cache
def source(resident, mobile):
    board = instrument.receiver.original.rails.layout()[0]
    anchor, near, _, far = instrument.PORTS
    assert near in board.neighbours[anchor] and far not in board.neighbours[anchor]
    assert far in board.neighbours[near] and board.degree[near] == board.degree[far]
    k = instrument.receiver.original.rails.candidate.KAPPA
    diag, hopping = instrument.receiver.original.terms(mobile, k)
    costs = []
    for address in (near, far):
        load = {anchor: resident, address: mobile}
        edges = {tuple(sorted((i, int(j)))) for i in load for j in board.neighbours[i]}
        value = sum(diag[load.get(i, 0), load.get(j, 0)] for i, j in edges)
        costs.append(
            value + sum((6 - board.degree[i]) * k * q**2 for i, q in load.items())
        )
    h = np.diag(costs)
    h[0, 1] = h[1, 0] = hopping[0, 1]
    odd = instrument.instruments(mobile)[1][0][:2, :2]
    # A permutation sends (near,far) to the first two odd-reader addresses
    # and parks the resident away from the reader's phase-parking operation.
    destination = instrument.PORTS[:2]
    park = next(v for v in board.neighbours[destination[1]] if board.degree[v] == 6)
    clear = {*destination, park}
    assert any(
        r not in clear and all(r not in board.neighbours[v] for v in clear)
        for r in range(board.n)
    )
    return h, odd


def readings(pair, stock, dwell):
    h, odd = source(*pair)
    initial = np.array(STOCKS[stock], dtype=complex)
    psi = expm(-1j * dwell * (h - h.trace() * np.eye(2) / 2)) @ (
        initial / np.linalg.norm(initial)
    )
    return abs(psi[0]) ** 2, (psi.conj() @ odd @ psi).real


def beat(pair):
    h, _ = source(*pair)
    _, link = instrument.receiver.original.terms(
        pair[1], instrument.receiver.original.rails.candidate.KAPPA
    )
    split = expm(-1j * np.pi / 2 * (link - link[0, 0] * np.eye(2)))
    arm = split.conj().T @ np.diag(np.exp(-1j * np.diag(h) * SHUTTER)) @ split[:, 0]
    return abs(arm[0]) ** 2


def acquire():
    rng = np.random.default_rng(261127)
    banks, distal, beats = [], [], []
    for pair in MATERIALS:
        h, _ = source(*pair)
        diag = h.diagonal() - 0.5
        beats.append(int(rng.binomial(N, probability(beat(pair)))))
        for stock in STOCKS:
            for dwell in DWELLS:
                near, odd = readings(pair, stock, dwell)
                parity_starts = int(rng.binomial(N, 1 / SCALE))
                near_starts = int(rng.binomial(N - parity_starts, probability(near)))
                nside = [near_starts, N - parity_starts - near_starts]
                clicks = [
                    int(rng.binomial(n, probability(d / (SCALE - 1))))
                    for n, d in zip(nside, diag, strict=True)
                ]
                parity_clicks = int(rng.binomial(parity_starts, probability(odd)))
                banks.append(
                    dict(
                        pair=pair,
                        stock=stock,
                        dwell=dwell,
                        parity_starts=parity_starts,
                        parity_clicks=parity_clicks,
                        sharp_starts=nside,
                        sharp_clicks=clicks,
                    )
                )
                if dwell == DWELLS[1] and stock in ("near", "far"):
                    releases = sum(clicks) + parity_clicks
                    visits, ticks = guide.parent.visits(rng, releases)
                    for probe in (1, 2):
                        setting = guide.Setting(
                            "matched", (6, 5, 5, 0), probe, (2, 2, 2)
                        )
                        arms = guide.count_panel(rng, visits, setting)
                        distal.append(
                            dict(
                                pair=pair,
                                stock=stock,
                                probe=probe,
                                releases=releases,
                                last_drain=ticks,
                                arms=arms,
                            )
                        )
    return dict(seed=261127, starts=N, beats=beats, banks=banks, distal=distal)


def enclosure(lo, hi):
    return [int(np.floor(lo * 1e6)) - 1, int(np.ceil(hi * 1e6)) + 1]


def audit(data):
    assert data["starts"] == N and len(data["beats"]) == len(MATERIALS)
    assert {(tuple(r["pair"]), r["stock"], r["dwell"]) for r in data["banks"]} == {
        (p, s, t) for p in MATERIALS for s in STOCKS for t in DWELLS
    }
    expected_distal = {
        (p, s, q) for p in MATERIALS for s in ("near", "far") for q in (1, 2)
    }
    assert len(data["distal"]) == len(expected_distal)
    assert {
        (tuple(r["pair"]), r["stock"], r["probe"]) for r in data["distal"]
    } == expected_distal
    u = guide.visits_before_drain()
    outcomes = []
    for i, pair in enumerate(MATERIALS):
        h, _ = source(*pair)
        cal = data["beats"][i] / N

        def gap(p):
            return np.arccos(np.clip(2 * p - 1, -1, 1)) / SHUTTER

        lower, upper = gap(cal + WIDTH), gap(cal - WIDTH)
        predicted_gap = h[1, 1] - h[0, 0]
        assert lower < predicted_gap < upper
        matched = [r for r in data["banks"] if tuple(r["pair"]) == pair]
        assert len(matched) == len(STOCKS) * len(DWELLS)
        by_stock = {}
        for row in matched:
            stock, dwell = row["stock"], row["dwell"]
            initial = np.array(STOCKS[stock])
            expect = initial @ h @ initial / (initial @ initial)
            sharp = sum(row["sharp_clicks"])
            total = sharp + row["parity_clicks"]
            assert sum(row["sharp_starts"]) + row["parity_starts"] == N
            assert 0 <= row["parity_clicks"] <= row["parity_starts"]
            assert all(
                0 <= c <= n
                for c, n in zip(row["sharp_clicks"], row["sharp_starts"], strict=True)
            )
            near, _ = readings(pair, stock, dwell)
            assert abs(expect / SCALE - total / N) < WIDTH
            assert (
                abs(
                    (np.array([near, 1 - near]) @ (np.diag(h) - 0.5)) / SCALE
                    - sharp / N
                )
                < WIDTH
            )
            if dwell == DWELLS[-1]:
                by_stock[stock] = total
        source_gap = (by_stock["far"] - by_stock["near"]) / N * SCALE
        source_interval = enclosure(
            source_gap - 2 * WIDTH * SCALE, source_gap + 2 * WIDTH * SCALE
        )
        probes = []
        for probe in (1, 2):
            setting = guide.Setting("matched", (6, 5, 5, 0), probe, (2, 2, 2))
            reward = guide.signed_readout(setting)
            per_release = np.sum(u * reward)
            assert abs(per_release + 1 / 15) < 1e-12
            samples = {}
            half_width = 0
            for row in data["distal"]:
                if tuple(row["pair"]) != pair or row["probe"] != probe:
                    continue
                assert row["releases"] == by_stock[row["stock"]]
                assert 0 < row["last_drain"] <= guide.PLAN["maximum_appointments"]
                assert len(row["arms"]) == 2 and all(
                    len(arm["starts"]) == len(arm["right"]) == 8 for arm in row["arms"]
                )
                walked = sum(sum(arm["starts"]) for arm in row["arms"])
                right = sum(sum(arm["right"]) for arm in row["arms"])
                assert all(
                    0 <= v <= n
                    for arm in row["arms"]
                    for v, n in zip(arm["right"], arm["starts"], strict=True)
                )
                samples[row["stock"]] = (right - walked / 2) / N
                half_width += abs(per_release) * WIDTH + row[
                    "releases"
                ] / N * guide.bound(reward, u, 12, row["releases"], 0.005)
            delta = samples["near"] - samples["far"]
            assert delta - half_width > 0
            assert (
                delta - half_width
                <= -per_release * predicted_gap / SCALE
                <= delta + half_width
            )
            probes.append(enclosure(delta - half_width, delta + half_width))
        outcomes.append(
            dict(
                pair=pair,
                ramsey_gap=enclosure(lower, upper),
                click_gap=source_interval,
                recoil_near_minus_far=probes,
                forecast_from_ramsey=enclosure(
                    lower / (15 * SCALE), upper / (15 * SCALE)
                ),
            )
        )
    return outcomes


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--acquire", action="store_true")
    if p.parse_args().acquire:
        data = acquire()
        (HERE / "records.json").write_text(json.dumps(data, indent=2) + "\n")
        (HERE / "results.json").write_text(json.dumps(audit(data), indent=2) + "\n")
    else:
        assert json.dumps(
            audit(json.loads((HERE / "records.json").read_text()))
        ) == json.dumps(json.loads((HERE / "results.json").read_text()))
        print(
            "Verified the three binding clocks, all source banks and twelve distal material banks."
        )
