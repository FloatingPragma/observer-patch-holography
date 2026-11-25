import json
from itertools import permutations, product

import numpy as np
import pytest
from scipy.linalg import expm

from . import centre, local_bases, measure


def epsilon():
    tensor = np.zeros((3, 3, 3))
    for permutation in permutations(range(3)):
        inv = sum(
            x > y for k, x in enumerate(permutation) for y in permutation[k + 1 :]
        )
        tensor[permutation] = (-1) ** inv / np.sqrt(6)
    return tensor


def test_hodge_response_in_every_complex_basis():
    rng = np.random.default_rng(392)
    for _ in range(6):
        raw = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
        u, _ = np.linalg.qr(raw)
        v, _ = np.linalg.qr(raw.T)
        data = measure.detector(epsilon(), v, u, u)
        expected = np.zeros((3, 3, 3))
        for a, b, c in measure.TUPLES:
            if b != c:
                expected[a, b, c] = abs((v @ u.conj().T)[a, 3 - b - c]) ** 2 / 6
        np.testing.assert_allclose(data, expected, atol=2e-15, rtol=0)
        wedge = np.sqrt(3) * epsilon().reshape(3, 9).T
        np.testing.assert_allclose(
            wedge.T @ np.kron(u, u) @ wedge,
            np.linalg.det(u) * u.conj(),
            rtol=0,
            atol=2e-15,
        )


def test_native_gates_and_physical_routes():
    doc = json.loads((measure.HERE / "counts.json").read_bytes())
    board = measure.carrier.federation(0)
    route_docs = doc["routes"]
    labels = np.arange(board.n)
    for record in route_docs:
        for i, j in record["pulses"]:
            assert j in board.neighbours[i]
            labels[[i, j]] = labels[[j, i]]
    for record in route_docs:
        np.testing.assert_array_equal(labels[record["triangle"]], record["source"])
        for i, j in zip(record["triangle"], record["triangle"][1:]):
            assert j in board.neighbours[i] and board.degree[i] == 5
    for record in reversed(route_docs):
        for i, j in reversed(record["pulses"]):
            labels[[i, j]] = labels[[j, i]]
    np.testing.assert_array_equal(labels, np.arange(board.n))
    h = measure.native_law.seam_term(measure.KAPPA)
    swap = np.eye(9)[
        [measure.native_law.SRC.index((b, a)) for a, b in measure.native_law.SRC]
    ]
    np.testing.assert_allclose(expm(-1j * np.pi * h), swap, rtol=0, atol=4e-14)
    for word in ("A", "B", "AB", "BA", "BAB"):
        u = measure.native_word(word)
        np.testing.assert_allclose(
            measure.native_word(word, operation="inverse"),
            u.conj().T,
            rtol=0,
            atol=5e-14,
        )
        np.testing.assert_allclose(
            measure.native_word(word, operation="transpose"), u.T, rtol=0, atol=5e-14
        )


def test_unconditioned_arrays_and_all_held_decisions():
    doc = json.loads((measure.HERE / "counts.json").read_bytes())
    for k, bank in enumerate(doc["banks"]):
        for mode, data in bank["records"].items():
            for rule, cal in bank["calibration"].items():
                redo = measure.score(data["counts"], cal)
                assert redo == data["scores"][rule]
                verdict = redo["covered"]
                if mode == "zero_contrast":
                    assert verdict is None
                elif mode == "native":
                    assert verdict is (rule == "inverse" or k == 3)
                else:
                    assert verdict is False
        broken = bank["records"]["native"]["counts"].copy()
        broken[0] -= 1
        with pytest.raises(ValueError, match="Incomplete"):
            measure.score(broken, bank["calibration"]["inverse"])


def test_split_centre_and_sharp_clock_columns():
    doc = json.loads((measure.HERE / "centre.json").read_bytes())
    board = measure.carrier.federation(0)
    clock, (x, y) = centre.clock_instruments()
    np.testing.assert_allclose(
        np.outer(clock, clock.conj()), np.ones((2, 2)) / 2, atol=3e-14, rtol=0
    )
    for instrument, signed_readout in ((x, [[0, 1], [1, 0]]), (y, [[0, -1j], [1j, 0]])):
        np.testing.assert_allclose(
            instrument.conj().T @ np.diag([1, -1]) @ instrument,
            signed_readout,
            atol=1e-13,
            rtol=0,
        )
    for route in doc["routes"]:
        labels = np.arange(board.n)
        for i, j in route["swaps"]:
            assert j in board.neighbours[i]
            labels[[i, j]] = labels[[j, i]]
        p, q = route["ancilla"], route["target"]
        assert labels[p] == doc["clock"][0] and labels[q] == q
        assert q in board.neighbours[p] and board.degree[p] == board.degree[q] == 5
    for theta in centre.ANGLES:
        for weight in centre.WEIGHTS:
            np.testing.assert_allclose(
                centre.endpoint_word(theta * weight),
                [1, 1, 1, np.exp(1j * theta * weight)],
                atol=3e-14,
                rtol=0,
            )
    for row in doc["rows"]:
        assert row["fits"]
        for r in row["counts"]:
            assert sum(r) == measure.STARTS
    for ref in doc["single_load_standards"]:
        z = np.exp(1j * ref["theta"] * centre.WEIGHTS[ref["colour"]])
        for bank, truth in zip(ref["counts"], (z.real, z.imag), strict=True):
            assert sum(bank) == measure.CAL
            assert abs((bank[0] - bank[1]) / measure.CAL - truth) < np.sqrt(
                2 * np.log(2 * 128 / 0.01) / measure.CAL
            )
    z = np.exp(2j * np.pi / 3)
    group = [
        exps
        for exps in product(range(3), repeat=3)
        if abs(np.prod([z**k for k in exps]) - 1) < 1e-12
    ]
    assert len(group) == 9 and (1, 1, 1) in group and (1, 0, 0) not in group


def test_four_native_bases_certify_singlet_from_local_counts():
    swap = np.eye(9).reshape(3, 3, 9).transpose(1, 0, 2).reshape(9, 9)
    collisions = np.zeros((9, 9), dtype=complex)
    for frame in local_bases.bases():
        compiled = local_bases.compile_basis(frame)
        for direction in compiled.conj():
            pair = np.kron(direction, direction)
            collisions += np.outer(pair, pair.conj())
    np.testing.assert_allclose(collisions, np.eye(9) + swap, atol=2e-13, rtol=0)

    basis = np.eye(27).reshape(3, 3, 3, 27)
    total = 3 * np.eye(27) + sum(
        basis.swapaxes(a, b).reshape(27, 27) for a, b in ((0, 1), (0, 2), (1, 2))
    )
    np.testing.assert_allclose(
        np.linalg.eigvalsh(total), [0] + [3] * 16 + [6] * 10, atol=2e-14, rtol=0
    )
    np.testing.assert_allclose(total @ epsilon().ravel(), 0, atol=1e-15, rtol=0)
    doc = json.loads((measure.HERE / "local.json").read_bytes())
    for case, experiment in doc["cases"].items():
        rows = experiment["bases"]
        assert len(rows) == 4
        # Inverse-beta endpoints depend on the numerical library. Compare each
        # interval in event units, to one ten-millionth of a count. The four
        # integer inventories and the physical fidelity thresholds are exact.
        event_budget = 1e-7
        yielded = local_bases.launch_yield(rows)
        recorded = experiment["launch_yield"].copy()
        np.testing.assert_allclose(
            yielded.pop("success_per_million_interval"),
            recorded.pop("success_per_million_interval"),
            rtol=0, atol=1e6 * event_budget / measure.STARTS,
        )
        assert yielded == recorded
        intervals = [local_bases.pair_intervals(row["counts"]) for row in rows]
        for row, actual in zip(rows, intervals, strict=True):
            enrolled = sum(row["counts"][2:])
            np.testing.assert_allclose(
                actual, row["equal_pair_intervals"], rtol=0, atol=event_budget / enrolled
            )
        # The reported limits must be precisely the prescribed transformation
        # of the archived endpoints, independently of inverse-CDF roundoff.
        assert local_bases.fidelity([r["equal_pair_intervals"] for r in rows]) == experiment["fidelity_interval"]
        limits = local_bases.fidelity(intervals)
        if case in ("idle", "common_AB"):
            assert limits[0] > 0.9997
            assert experiment["launch_yield"]["success_per_million_interval"][0] > 33
        else:
            assert limits[1] < 0.9


def test_archived_acquisitions_replay_and_reject_changed_probabilities():
    from model.candidate_laws.interval_echo.openwave.electron_field.replay import (
        compare,
    )

    for module, name in ((measure, "counts.json"), (centre, "centre.json"), (local_bases, "local.json")):
        archive = json.loads((measure.HERE / name).read_bytes())
        data = module.run(archive["acquisition"]["launch_probabilities"])
        compare(json.loads(json.dumps(data)), archive, name)
    draw, _ = measure.sampler(1, [[0.5, 0.5]])
    with pytest.raises(ValueError, match="archived acquisition"):
        draw(measure.CAL, [0.5 + 1e-9, 0.5 - 1e-9])
