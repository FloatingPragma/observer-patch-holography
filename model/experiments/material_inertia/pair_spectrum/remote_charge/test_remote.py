import json
from math import log, sqrt

import networkx as nx
import numpy as np
from model.executable_model.carrier import federation
from . import measure
from .geometry import guides


def test_full_three_site_charge_source_and_causal_remote_observable():
    x = np.array([[0, 1], [1, 0]])
    z = np.diag([1, -1])
    eye = np.eye(2)

    def kron(a, b, c):
        return np.kron(a, np.kron(b, c))

    u = measure.unitary_layers("native")
    obs = kron(eye, eye, z)
    np.testing.assert_allclose(u.conj().T @ obs @ u, -kron(eye, z, x), atol=2e-13)
    u2 = u @ u
    np.testing.assert_allclose(u2.conj().T @ obs @ u2, kron(z, x, x), atol=2e-13)
    for name, targets in measure.TARGETS.items():
        q = sum(n - 1 for n in targets)
        for setting in ("sign", "parity"):
            p = measure.calibration(targets, setting, "native")[:, 0]
            final = abs(q) if setting == "parity" else (1 - q) / 2
            np.testing.assert_allclose(p, [0.5, 0.5, final], atol=2e-13)
            np.testing.assert_allclose(
                measure.calibration(targets, setting, "cut_bridge")[:, 0],
                [0.5, 0.5, 0.5],
                atol=2e-13,
            )

    # Changing the order of the compensated target or inserting ordinary
    # background ports changes only target-only phases in the remote effects.
    for setting in ("parity", "sign"):
        host = measure.calibration((1,), setting, "native")
        for neutral in ((0, 2), (2, 0), (1, 0, 1, 2), (2, 1, 0, 1)):
            np.testing.assert_allclose(
                measure.calibration(neutral, setting, "native"), host, atol=2e-13
            )


def test_routed_inputs_have_local_contacts_and_never_enter_receiver_owner():
    board = federation(0)
    net = nx.Graph(board.seams)
    options = guides()
    census = json.loads(measure.SOURCE.read_bytes())
    record = json.loads((measure.HERE / "results.json").read_bytes())
    assert len(census["profiles"]) == len(record["fixtures"]) == 60
    for profile in record["fixtures"]:
        used = set()
        for site in profile["sites"]:
            path = site["guide"]
            assert all(net.has_edge(a, b) for a, b in zip(path, path[1:]))
            assert all(
                board.degree[a] == board.degree[b]
                for a, b in zip(path[::2], path[1::2])
            )
            assert net.has_edge(path[1], site["contact"]) and board.degree[path[1]] == 5
            inside = {p // 12 for p in path}
            assert len(inside) == 2 and not inside & (used | set(profile["owners"]))
            used.update(inside)
            routes = next(
                t for _, g, j, t in options if tuple(path) == g and j == site["contact"]
            )
            for owner in profile["owners"]:
                for p in range(12 * owner, 12 * owner + 12):
                    route = routes[p]
                    assert route[0] == site["contact"] and route[-1] == p
                    assert all(
                        z // 12 != path[-1] // 12 and z not in path for z in route
                    )


def test_original_denominators_complete_and_remote_controls_fail():
    r = json.loads((measure.HERE / "results.json").read_text())
    draws = r["calibration_shots_per_setting"]
    margin = sqrt(log(2 * 120 / 0.01) / (2 * draws))
    for panel in r["calibrations"]:
        for obs in panel["observations"]:
            assert sum(obs["terminal_counts"]) == draws
            assert (
                abs(obs["terminal_counts"][0] / draws - obs["probability_word20"])
                <= margin
            )
        p = panel["observations"][-1]["probability_word20"]
        if panel["guide"] == "exchange_only_seams":
            assert abs(p - 0.5) < 0.0003
        if panel["source"] == "excess" and panel["quadrature"] == "parity":
            assert ("native" == panel["guide"]) == (p > 0.99)
    for launch in r["equivalent_source_banks"]:
        for group in launch["panels"].values():
            assert len(group) == 60
            for row in group:
                assert (
                    sum(row["signed_source_counts"]) + row["initial_filter_refusals"]
                    == r["initial_pair_shots_per_setting"]
                )
        assert (
            max(
                row["parity_dark_probability_interval"][1]
                + launch["native_per_event_error"]
                for row in launch["panels"]["native"]
            )
            < 0.001
        )
