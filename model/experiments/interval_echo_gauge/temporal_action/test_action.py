import gzip
import json
from math import log, sqrt

import numpy as np
import pytest
from scipy.linalg import expm

from . import action
from .material import geometry, reader, source
from .measure import HERE, SHOTS, constraints, run
from .rails import GRID, Meter, fringe, hardware, read_reference, significant


def deposited():
    return json.loads(gzip.decompress((HERE / "records.json.gz").read_bytes()))


def test_remote_source_census_and_signed_cycle():
    d, curl, _, pin, face, hole = geometry()
    assert [pin[x] for x in face] == [0, 1, 5] and pin[hole] == 7
    reader_setup = reader()
    assert reader_setup == deposited()["receiver"]
    assert read_reference() == 17 and hardware() == (12, 13)
    assert all(
        17 not in path and 12 not in path and 13 not in path
        for path in reader_setup["paths"]
    )
    idle, ji, _ = source("idle")
    closed, jc, report = source("loop")
    reverse, jr, _ = source("reverse")
    np.testing.assert_array_equal(closed, idle)
    np.testing.assert_array_equal(reverse, idle)
    np.testing.assert_array_equal(jc[0], -2 * curl[0])
    np.testing.assert_array_equal(jr, -jc)
    assert not ji.any() and len(report["contacts"][0]) == 3
    for case in action.CASES:
        rho, j, _ = source(case)
        np.testing.assert_array_equal(2 * np.diff(rho, axis=0), -j @ d)


def test_phase_decoder_and_native_data_pulses():
    meter = Meter()
    values = [(3, -17, 12), (33, 27), (-2, 9, 11), (1 << 35,)]
    assert meter.read(values) == (1 << 35) + 81
    assert meter.read([(-1, 1 << 35), (7,), (-127, 101)]) == -(1 << 35) + 7 - 127 * 101
    assert meter.read([]) == 0
    with pytest.raises(OverflowError):
        meter.read([(GRID // 2,)])
    # Direct propagation of the two amplitudes, with independent recombiners.
    x = np.array([[0, 1], [1, 0]])
    gate = 32.5 * np.eye(2) - x / 2
    private = np.diag([0.0, 32])
    quarter = expm(-1j * np.pi * private / 64)
    analyzer = quarter @ expm(-1j * np.pi * gate / 2) @ quarter
    for value in (-123, 0, 73):
        for altered in (False, True):
            state = analyzer @ np.array([1, 0])
            for bit in significant(value):
                angle = np.sign(value) * np.pi / 2 * 2**bit / 512
                t = (-angle / 32) % (2 * np.pi / 32)
                state = expm(-1j * t * private * (0.5 if altered else 1)) @ state
            direct = [
                abs((analyzer @ (phase * state))[0]) ** 2
                for phase in (np.array([1, 1]), np.array([1, -1j]))
            ]
            np.testing.assert_allclose(fringe(value, 512, altered), direct, atol=2e-13)


def test_independent_held_front_and_window_failures():
    archive = deposited()
    rows = archive["settings"]
    c = geometry()[1]
    loop, idle = (np.array(rows["native_" + case]["B"]) for case in ("loop", "idle"))
    delta = loop - idle
    normal = c @ c.T
    np.testing.assert_array_equal(delta[:2], 0)
    np.testing.assert_array_equal(delta[2], -(action.L // 2) * normal[:, 0])
    np.testing.assert_array_equal(
        4 * delta[3], (8 * np.eye(20, dtype=int) - normal) @ delta[2]
    )
    assert delta[3, 3] == action.L // 8
    for case in action.CASES:
        h = {
            k: np.array(v)
            for k, v in rows["native_" + case].items()
            if k in {"A", "E", "B", "J", "rho", "phi"}
        }
        assert not any(constraints(h).values())
        stationary = rows["native_" + case]["stationarity"]
        assert action.independent(h) == stationary["action"]
        for key, shape in [("A", (4, 30)), ("phi", (5, 12))]:
            assert np.shape(stationary["differences"][key]) == shape
            assert not np.any(stationary["differences"][key])
    expected = {
        "off_shell": (7, 4),
        "no_kinetic": (65, 10),
        "no_magnetic": (65, 0),
        "reversed_current": (2, 0),
        "no_port_charge": (0, 10),
    }
    for key, failures in expected.items():
        measured = archive["challenges"][key]["nonzero"]
        assert (measured["A"], measured["phi"]) == failures
    assert archive["gauge"]["packet"]["action"] == 227472


def test_seeded_binomial_counts_and_complete_replay():
    archive = deposited()
    assert run(archive["sampling_probabilities"]) == archive
    assert len(archive["fringes"]) == 720
    bound = sqrt(log(2 * 1441 / 0.01) / (2 * SHOTS))
    metadata = archive["provenance"]
    assert (
        metadata["seed"] == 662901
        and metadata["hoeffding"]["number_of_frequencies"] == 1441
    )
    assert metadata["hoeffding"]["half_width"] == bound
    rng = np.random.Generator(np.random.PCG64(metadata["seed"]))
    crossing_p = archive["source_records"]["open_swap"]["half_probabilities"][1]
    assert rng.binomial(SHOTS, crossing_p) == archive["half_source_original_starts"][1]
    for row, p in zip(
        archive["fringes"], archive["sampling_probabilities"], strict=True
    ):
        np.testing.assert_array_equal(rng.binomial(SHOTS, p), row[-1])
        assert np.max(abs(np.array(row[-1]) / SHOTS - p)) < bound
    quiet, swap = archive["half_source_original_starts"]
    assert quiet + swap == SHOTS and abs(swap / SHOTS - 0.5) < bound
