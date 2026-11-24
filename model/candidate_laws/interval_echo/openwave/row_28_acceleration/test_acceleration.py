"""Independent native and observational checks on the expansion hypothesis."""

import json
from collections import deque

import numpy as np
import pytest
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from model.experiments.interval_echo_consensus.asynchronous import SharpService
from .experiment import HERE


def test_independent_sharp_controller_preserves_the_deposited_addresses():
    board = federation(0)
    service = SharpService(range(board.n), board.seams)
    rng = np.random.default_rng(280210)
    loads = (2 * rng.integers(0, 2, board.n, dtype=np.int8)).tolist()
    collars = [0] * len(board.seams)
    with np.load(HERE / "records.npz") as bank:
        times = bank["epoch_280210"]
        for epoch in range(1, int(times.max()) + 1):
            draws = rng.random(len(collars))
            chosen = [
                i
                for i in range(len(collars))
                if draws[i] < 0.25 and service.enabled(i, loads, collars)
            ]
            assert chosen == np.flatnonzero(times == epoch).tolist()
            for edge in chosen:
                service.commit(edge, loads, collars)
        assert collars == bank["sign_280210"].tolist()
        assert not any(service.ready(i, loads, collars) for i in range(len(collars)))


def test_native_pulse_restores_every_addressed_seam_and_material():
    law = candidate.native_law
    kinetic, contrast = law.local_term_parts()
    pulse = expm(-1j * np.pi * (kinetic + 8 * contrast))
    for input_word in ((0, 2), (2, 0), None):
        if input_word is None:
            mem = np.eye(1 + law.N_LAB)[0]
        else:
            column = law.echo()[:, law.SRC.index(input_word)].reshape(-1, law.N_SRC)
            mem = column[:, law.SRC.index((1, 1))]
        for loads in ((0, 1), (0, 2), (2, 0), (1, 0)):
            start = np.kron(mem, np.eye(law.N_SRC)[law.SRC.index(loads)])
            arrival = pulse @ start
            detected = loads if input_word else loads[::-1]
            ptr = arrival.reshape(-1, law.N_SRC)[:, law.SRC.index(detected)]
            assert np.vdot(ptr, ptr).real == pytest.approx(1)
            np.testing.assert_allclose(pulse.conj().T @ arrival, start, atol=3e-13)


def test_room_for_competing_scale_and_period_functions():
    primary = json.loads((HERE / "results.json").read_text())
    scales = [
        s
        for p in primary["panels"]
        if p["attempt_probability"] == 0.25
        for s in p["standard"]
        if s
    ]
    assert len(scales) == 67
    assert sum(s["initial_acceleration"] > 0 for s in scales) == 27
    assert all(x["rejected"] for s in scales for x in s["held"])
    assert all(
        s["scales"][-1] == s["scales"][-2] or s["scales"][-1] >= s["scales"][-2]
        for s in scales
    )
    signals = json.loads((HERE / "signals.json").read_text())["cases"]
    tests = [m["pairs"][0] for p in signals if not p["dilute"] for m in p["metrics"]]
    assert sum(r["missed_beyond_one_round"] for r in tests) == 22016
    for p in signals:
        for m in p["metrics"]:
            last = m["pairs"][-1]
            assert last["both"] == last["tested"] == last["period_histogram"][4]


def test_event_graph_replays_the_same_first_stamp_and_diverging_second_stamps():
    level, seed = 1, 280222
    board = federation(level)
    lookup = {tuple(sorted(e)): i for i, e in enumerate(board.seams)}
    with np.load(HERE / "records.npz") as bank:
        commit = bank[f"epoch_{seed}"]
    # Dynamic programming over actual spacetime events, allowing waits;
    # do not assume that first-arrived paths dominate later paths.
    for emit, expected in ((0, [2, 2]), (4, [6, 59])):
        sources = {265}
        first = {}
        for t in range(emit, 60):
            for port in (267, 286):
                if port in sources:
                    first.setdefault(port, t)
            nxt = set(sources)
            for u in sources:
                for v in board.neighbours[u]:
                    close = commit[lookup[tuple(sorted((u, v)))]]
                    if close == 0 or t < close:
                        nxt.add(v)
            sources = nxt
        assert [first[v] for v in (267, 286)] == expected
    assert all(len(board.neighbours[v]) >= 5 for v in (265, 267, 286))


def test_boundary_and_capacity_witness_without_component_labels():
    board = federation(1)
    with np.load(HERE / "records.npz") as bank:
        broken = bank["epoch_280222"] > 0
    live = [[] for _ in range(board.n)]
    for (u, v), cut in zip(board.seams, broken, strict=True):
        if not cut:
            live[u].append(v)
            live[v].append(u)
    for origin, expected in ((0, 41), (265, 72)):
        q = deque([origin])
        distance = {origin: 0}
        while q:
            u = q.popleft()
            for v in live[u]:
                if v not in distance:
                    distance[v] = 1 + distance[u]
                    q.append(v)
        assert (len(distance), max(distance.values())) == (468, expected)
        boundary = [
            (a, b) for a, b in board.seams if (a in distance) != (b in distance)
        ]
        assert len(boundary) == 25


def test_echo_uses_only_origin_stamps_and_agrees_with_separate_return_propagation():
    from .signals import adjacency, arrival

    (example,) = [
        row
        for row in json.loads((HERE / "radar.json").read_text())
        if row["seed"] == 280222 and not row["dilute"]
    ]
    board = federation(1)
    edges = np.asarray(board.seams)
    with np.load(HERE / "records.npz") as bank:
        closing = bank["epoch_280222"]
    graphs = [
        adjacency(board.n, edges[(closing == 0) | (closing > t)])
        for t in range(1 + int(closing.max()))
    ]
    for emit, expected in zip(example["emissions"], example["returns"], strict=True):
        outward = arrival(graphs, [265], emit)[0]
        observed = [
            int(arrival(graphs, [target], int(outward[target]))[0, 265])
            for target in (267, 286)
        ]
        assert observed == expected
    assert example["returns"][:2] == [[4, 4], [8, 114]]
    assert example["transmitting_control"][:2] == [[4, 4], [8, 8]]


def test_same_source_bit_crosses_both_orientations_of_every_calibrated_cut():
    from .cut_signal import SHOTS

    data = json.loads((HERE / "cut_signal.json").read_text())
    native = [r for r in data["records"] if r[2:5] == [None, 8, 2]]
    exchange = [r for r in data["records"] if r[2:5] == [None, 0, 2]]
    assert {(r[0], r[1]) for r in native} == {
        (int(k), s) for k in data["endpoints"] for s in ("left", "right")
    }
    assert min((r[6] - r[5]) / SHOTS for r in native) > 0.90
    assert max(abs(r[6] - r[5]) / SHOTS for r in exchange) < 0.10
    # Reversal of |+2> is orthogonal to it; its even and odd components
    # accumulate phases 32t and 33t whenever the special 11 word is absent.
    t = np.pi / 32
    expected = 0.5 * (1 - 0.5 * (np.cos(32 * t) + np.cos(33 * t)))
    p = sum(r[6] for r in native) / (SHOTS * len(native))
    assert abs(expected - p) < 0.001


def test_two_receiver_code_needs_coherence_and_both_calibrated_durations():
    import sympy as sp

    # Each preserved sign parity has a different null factor in the product
    # of quiet amplitudes. This checks exact zero independent of simulation tolerances.
    for energy in (32, 33):
        assert sp.prod(sp.cos(sp.pi * energy / (2 * d)) for d in (32, 33)) == 0
    data = json.loads((HERE / "paired_probe.json").read_text())
    for r in data["records"]:
        assert sum(r["joint_counts"]) == data["shots"]
        if r["bit"] == 0:
            assert r["joint_counts"][0] == data["shots"]
        elif r["duration_scale"] == 1 and r["kappa"] == 8:
            if r["coin"] is None:
                assert r["joint_counts"][0] == 0 and r["p_quiet"] == 0
            else:
                assert r["p_quiet"] == 0.5
        elif r["bit"] == 1:
            assert r["p_quiet"] > 0
