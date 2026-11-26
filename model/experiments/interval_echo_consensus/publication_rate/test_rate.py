import gzip
from itertools import product
import json

import numpy as np

from . import apparatus as a
from . import measure as m
from . import wait


def test_protected_instrument_on_all_matrix_units():
    for grouped in (False, True):
        blocks = a.repair_blocks(grouped)
        keys = [abs(u - v) if grouped else u - v for u, v in a.LAW.SRC]
        for i, j in product(range(9), repeat=2):
            unit = np.zeros((9, 9))
            unit[i, j] = 1
            out = sum(k @ unit @ k.conj().T for k in blocks.values())
            np.testing.assert_allclose(out, unit * (keys[i] == keys[j]), atol=3e-15)
        sharp = sum(k.conj().T @ k for k in blocks.values())
        np.testing.assert_allclose(sharp, np.eye(9), atol=3e-15)


def test_counts_and_interference_controls():
    data = m.acquire()
    assert m.plain(data) == gzip.decompress((m.HERE / "records.json.gz").read_bytes())
    for row in data["sharp"] + data["panels"]:
        assert sum(map(sum, row["terminal_ports"])) == m.N
    for row in data["sharp"]:
        assert np.sum(row["terminal_ports"], axis=0)[row["address"]] == m.N
    assert len(data["panels"]) == 300
    report = m.report(data)
    for mode, error in report["max_held_probability_error"].items():
        assert error > 0.3 if mode == "latched" else error < 0.006
    for n, hand in report["clock_phases"].items():
        assert abs(hand - int(n) * report["clock_phases"]["1"]) < 0.008
    assert data["raw_sharp_and_mixed"] == [[0, m.N, 0]] * 4
    assert data["free_returns"] == [[m.N, 0], [0, m.N], [m.N, 0]]


def test_independent_clock_forecasts_on_uninterrupted_flight():
    data = wait.acquire()
    assert m.plain(data) == gzip.decompress((m.HERE / "waits.json.gz").read_bytes())
    calibration = json.loads(gzip.decompress((m.HERE / "records.json.gz").read_bytes()))
    rows = wait.forecast(data, calibration)
    for r in rows:
        assert (
            max(
                abs(r[k] - r["single_flight"])
                for k in ("smallest_bright", "largest_bright")
            )
            < 0.002
        )
    assert rows[-1]["smallest_bright"] - rows[-1]["continuous_semigroup"] > 0.41


def test_raw_echo_fails_unitality_on_the_full_algebra():
    ket = a.LAW.echo()[:, : a.LAW.N_SRC].reshape(-1, a.LAW.N_SRC, a.LAW.N_SRC)
    raw_identity = sum(k @ k.conj().T for k in ket)
    row = a.LAW.SRC.index((1, 1))
    assert abs(raw_identity[row, row] - 3) < 1e-13
    for actual, expected in zip(
        raw_identity.diagonal(), [1, 1, 0, 1, 3, 1, 0, 1, 1], strict=True
    ):
        assert abs(actual - expected) < 1e-13
