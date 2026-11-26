import gzip
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from . import pump, statistics_control as ctl


def test_antisymmetric_hops_are_the_exterior_square_with_native_diagonal():
    words = pump.stocks(2)
    holes = [tuple(i for i, v in enumerate(w) if v == 0) for w in words]
    # Independent construction on two labeled vacancies followed by antisymmetrization.
    one = ctl.hopping(1, False)
    one -= np.diag(one.diagonal())
    two = np.kron(one, np.eye(4)) + np.kron(np.eye(4), one)
    sectors = []
    for parity in (1, -1):
        embedding = np.zeros((16, 6))
        for j, (a, b) in enumerate(holes):
            embedding[4 * a + b, j] = 1 / np.sqrt(2)
            embedding[4 * b + a, j] = parity / np.sqrt(2)
        expected = embedding.T @ two @ embedding
        h = ctl.hopping(2, parity == -1)
        np.testing.assert_allclose(h - np.diag(h.diagonal()), expected, atol=3e-15)
        sectors.append(h)
    np.testing.assert_array_equal(sectors[0].diagonal(), sectors[1].diagonal())
    ports = pump.route_certificate()["ports"]
    # Sorting the Fock basis by physical port number only rephases the six words.
    gauges = np.array([1 if ports[a] < ports[b] else -1 for a, b in holes])
    physical = sectors[0].copy()
    for i, j in zip(*np.nonzero(np.triu(physical, 1)), strict=True):
        left, right = (set(ports[k] for k in holes[v]) for v in (i, j))
        a, b = sorted(left ^ right)
        spectator, = left & right
        physical[i, j] *= -1 if a < spectator < b else 1
        physical[j, i] = physical[i, j]
    np.testing.assert_allclose(gauges[:, None] * sectors[1] * gauges, physical)
    np.testing.assert_array_equal(ctl.hopping(1, False), ctl.hopping(1, True))
    assert ctl.direct_route(False) == 0 and ctl.direct_route(True) == 1


def test_both_arms_against_explicit_microsteps():
    n = 24
    dt = 1 / 16
    duration = 3
    source = np.eye(6)[:, pump.index_pair(pump.PAIRS[0])]
    _, pairs, centre = pump.source_geometry(2)
    for detuning, anti in itertools.product((16, -16), (False, True)):
        flight = expm(-1j * dt * ctl.hopping(2, anti))
        pulses = []
        for j in range(n):
            d = detuning * centre - pairs @ pump.aperture((j + 0.5) * 2 * np.pi / n)
            kick = np.diag(np.exp(-0.5j * dt * d))
            pulses.append(kick @ flight @ kick)
        def arm(sequence):
            u = np.eye(6)
            for j in sequence:
                for _ in range(2):
                    u = pulses[j] @ u
            return u
        f = arm(range(n))
        r = arm([*range(n // 2), *reversed(range(n // 2))])
        row = ctl.run_case(detuning, dt, n, duration, anti)
        assert abs(row["dark_once"] - ctl.read_exits(f, r, source)) < 1e-12
        assert abs(row["dark_twice"] - ctl.read_exits(f @ f, r @ r, source)) < 1e-12


def test_full_duration_controls_against_held_native_receipt():
    folder = Path(__file__).parent
    with gzip.open(folder / "records.json.gz", "rt") as handle:
        original = json.load(handle)["held"]
    capture = json.loads((folder / "statistics_control.json").read_text())
    live = ctl.run()
    assert list(live["ports"]) == capture["ports"]
    assert list(live["initial_stock"]) == capture["initial_stock"]
    assert live["direct_sharp_dark"] == capture["direct_sharp_dark"]
    for a, b in zip(live["pumps"], capture["pumps"], strict=True):
        for field in ("detuning", "microstep", "appointments", "duration", "statistics"):
            assert a[field] == b[field]
        for field in ("dark_once", "dark_twice", "pinned_probabilities", "crossing_probabilities"):
            np.testing.assert_allclose(a[field], b[field], rtol=0, atol=2e-8)
    for case in original:
        native, fermion = [
            x for x in capture["pumps"]
            if x["microstep"] == case["microstep"] and x["detuning"] == case["detuning"]
        ]
        assert native["statistics"] == "native" and fermion["statistics"] == "antisymmetric"
        held = next(r for r in case["records"] if r["duration"] == native["duration"])
        np.testing.assert_allclose(
            [native["dark_once"], native["dark_twice"], min(native["pinned_probabilities"]),
             min(native["crossing_probabilities"])],
            [held["once"]["dark"][0], held["twice"]["dark"][0], held["min_pinned"][0],
             held["min_crossing_shared_two_pairs"][0]], atol=2e-8, rtol=0,
        )
        assert (fermion["dark_once"] > 0.98) == (native["detuning"] < 0)
        assert max(native["dark_twice"], fermion["dark_twice"]) < 0.02
        assert min(fermion["pinned_probabilities"]) > 0.98
