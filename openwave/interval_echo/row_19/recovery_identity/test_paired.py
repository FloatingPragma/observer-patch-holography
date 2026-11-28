from copy import deepcopy
import json

import numpy as np
import pytest
from scipy.sparse.linalg import expm_multiply

from openwave.interval_echo.row_19.recovery_identity import paired_arrivals as arrivals
from openwave.interval_echo.row_19.recovery_identity import paired_relay as pair
from openwave.interval_echo.row_19.recovery_identity import paired_survey as radar
from openwave.interval_echo.row_19.recovery_identity import relay


def test_visible_and_hidden_scores_and_trace_bound():
    d = json.loads(pair.RECEIPT.read_text())
    report = pair.decode(d)
    assert report == d["measurement"]
    lo, hi = report["modular_contrast"]
    assert lo > d["source"]["expected_charge"] - 1e-8 and lo > pair.H**2 and hi > lo
    for i, arm in enumerate(("vacancy", "excess")):
        limits = report["receiver_distances"][arm]["source_distance"]
        assert limits[0] <= d["population_trace_distances"][i] <= limits[1] < 2.3e-9
    p = np.array([3 * np.eye(4) / 16, np.eye(4) / 16])
    delta = pair.H * np.einsum(
        "i,ibjk->bjk", d["source"]["direction"], pair.meter.head.joint_basis()
    )
    vals = np.linalg.eigvalsh(p + delta)
    # P is scalar within each parity block and commutes with R.
    direct = np.sum(
        (vals - p.diagonal(axis1=1, axis2=2))
        * np.log(vals / p.diagonal(axis1=1, axis2=2))
    )
    assert lo < direct < hi
    c = np.asarray(
        json.loads(pair.meter.RECEIPT.read_text())["gates"]["0.0"][
            "pauli_coefficients"
        ]["original"]
    )
    assert abs(c @ np.asarray(d["source"]["direction"])) > 0.01
    for bad in ("readers", "source_counts"):
        damaged = deepcopy(d)
        damaged[bad][0][0][0] = [] if bad == "readers" else -1
        with pytest.raises(ValueError):
            pair.decode(damaged)


def test_radar_check_and_direct_terminal_round_trip():
    d = json.loads(radar.RECEIPT.read_text())
    data = radar.decode(d)
    assert data == d["measurements"]
    assert sum(r["source_compatible"] is False for r in data["surveys"].values()) == 1
    for r in data["surveys"].values():
        assert r["blank_control"][0] < 0 < r["blank_control"][1]
    delivery = json.loads(pair.RECEIPT.read_text())
    paths = delivery["paths"]
    board = relay.FORMATION.federation(1)
    sites = [*paths[0], paths[1][0]]
    h = relay.RANGING.generator(board, (*paths[0][:2], *paths[1][:2]), 1, 8, False)
    state = pair.unpack(delivery["installed"][1])[0]

    # Earliest round trip has one propagation in each direction. Build both
    # from the original two-head solver and project the *returned* load.
    def one(first, second):
        cols = np.zeros((board.n * 4, 4))
        cols[4 * first : 4 * first + 4] = np.eye(4)
        w = expm_multiply(-0.5j * h, cols, traceA=0)
        return w[4 * second : 4 * second + 4]

    a = one(sites[1], sites[0]) @ one(sites[0], sites[1])
    expected = np.trace(a @ state @ a.conj().T).real
    n = d["banks"]["1,1"][0][0][0]
    assert abs(expected - n / radar.SHOTS) < radar.precision.radius(n, radar.L) + 2e-13


def test_arrival_flag_and_control():
    d = json.loads(arrivals.RECEIPT.read_text())
    m = arrivals.decode(d)
    for q in ("-1", "1"):
        assert m[q]["control"]["event_lower_bound"] == 0
        assert d["distances"][q]["complete"] >= d["distances"][q]["forgotten"]
        assert d["distances"][q]["complete"] >= m[q]["source"]["event_lower_bound"]
    bad = deepcopy(d)
    bad["banks"]["1"][0][0][0][0] += 1
    with pytest.raises(ValueError, match="Incomplete"):
        arrivals.decode(bad)


def test_arrival_contrasts_from_independent_generator_and_polynomial():
    from openwave.interval_echo.row_19.recovery_identity.bridge_dynamics import flight

    d = json.loads(arrivals.RECEIPT.read_text())
    original = json.loads(pair.RECEIPT.read_text())
    board = relay.FORMATION.federation(1)
    route = relay.FORMATION.region(board)
    paths = original["paths"]
    owners = np.arange(board.n) // 12
    first, second = (x[0] // 12 for x in paths)
    masks = (owners == first, owners == second, (owners != first) & (owners != second))
    direction = -pair.H * np.einsum(
        "i,ibjk->bjk", original["source"]["direction"], pair.meter.head.joint_basis()
    )
    for q in (-1, 1):
        differences = np.zeros((3, 4, 4), complex)
        for parity in range(2):
            u = flight(board, route, paths, q, parity)
            for k, mask in enumerate(masks):
                v = u[mask].reshape(-1, 4, 4)
                differences[k] += np.einsum(
                    "pca,ab,pdb->cd", v, direction[parity], v.conj(), optimize=True
                )
        target = d["distances"][str(q)]
        tol = 3e-13
        np.testing.assert_allclose(
            abs(np.linalg.eigvalsh(differences)).sum() / 2,
            target["complete"],
            atol=tol,
            rtol=0,
        )
        np.testing.assert_allclose(
            abs(np.linalg.eigvalsh(differences.sum(axis=0))).sum() / 2,
            target["forgotten"],
            atol=tol,
            rtol=0,
        )
        measures = d["measurements"][str(q)]["source"]["flags"]
        for diff, (lo, hi) in zip(differences, measures, strict=True):
            assert lo - tol <= diff.trace().real <= hi + tol


def test_flagged_access_and_observational_budgets():
    from itertools import product
    from openwave.interval_echo.row_19.recovery_identity import access, retention_budget

    receipt = json.loads(access.RECEIPT.read_text())
    for p in receipt["panels"].values():
        rows = np.array(p["boards"], dtype=int)
        n = int(rows[:, -1].sum())
        assert n == receipt["boards_per_setting"]
        for k, s in enumerate(p["settings"]):
            assert int(rows[:, -1] @ rows[:, k + 1]) == s["failed"]
            assert int(rows[:, -1] @ rows[:, 0]) == s["cuts"]
            assert (
                s["observed_score"][0]
                <= s["predicted_conditional_score"]
                <= s["observed_score"][1]
            )
            assert (
                abs(s["observed_failure"] - s["population_failure"])
                < 4 * s["board_standard_error"]
            )
    # Enumerate all four witness flags at a rational occupancy, conditioning
    # on the two distinct ends. Each weight is counted before any routing.
    p = 0.25
    choices = [
        (bits, p ** sum(bits) * (1 - p) ** (4 - sum(bits)))
        for bits in product((0, 1), repeat=4)
    ]
    assert sum(prob for bits, prob in choices if not any(bits)) == (1 - p) ** 4
    gain = retention_budget.run()["galaxy_min_failed_fraction_gain"]
    assert access.homogeneous_ceiling(receipt["witness_counts"][0], 0.16, gain) is None
    max_p = access.homogeneous_ceiling(receipt["witness_counts"][1], 0.16, gain)
    assert 0.06209 < max_p < 0.06210
