import hashlib
import json

import networkx as nx
import numpy as np
import pytest

from model.candidate_laws.interval_echo.candidate import KAPPA
from model.executable_model import cell, federation
from model.experiments.interval_echo_gauge.singlet_flight import measure, orientation

P, WORD, ANTIPODES, CAL = measure.colour_meter()


def test_centre_and_energy_identify_same_triplet_on_all_ports():
    energy, modes = np.linalg.eigh(cell().laplacian().toarray() / 2)
    triplet = modes[:, np.isclose(energy, (5 - np.sqrt(5)) / 2)]
    np.testing.assert_allclose(P, triplet @ triplet.T, atol=1e-12)
    np.testing.assert_allclose(P.diagonal(), 0.25)
    expected = np.eye(12) + (np.exp(2j * np.pi / 3) - 1) * P
    assert np.linalg.norm(WORD - expected, 2) < 8e-5
    np.testing.assert_allclose(WORD.conj().T @ WORD, np.eye(12), atol=5e-12)


@pytest.mark.parametrize("level", range(3))
def test_all_gateways_and_exact_native_tangent_equation(level):
    board = federation(level)
    by_owner = [[] for _ in range(board.cells)]
    for a, b in board.external:
        by_owner[a // 12].append((a, b))
        by_owner[b // 12].append((b, a))
    # Verify the incidence witness directly in the actual numerical band basis.
    # Glue endpoint coefficients with +/- orientation for an empty background.
    link = nx.Graph()
    constrained = set()
    for ends in by_owner:
        ports = [a % 12 for a, b in ends]
        assert len({b // 12 for a, b in ends}) == len(ports) == 3
        columns = P[:, ports]
        _, sv, vt = np.linalg.svd(columns, full_matrices=False)
        pair_ids = [tuple(sorted((a, b))) for a, b in ends]
        if sv[-1] > 1e-8:
            constrained.update(pair_ids)
        else:
            assert sv[-2] > 1e-8
            kernel = vt[-1]
            zero, left, right = np.argsort(abs(kernel))
            assert abs(kernel[zero]) < 1e-12
            np.testing.assert_allclose(kernel[left], kernel[right], atol=1e-12)
            constrained.add(pair_ids[zero])
            link.add_edge(pair_ids[left], pair_ids[right])
    assert all(not c.isdisjoint(constrained) for c in nx.connected_components(link))
    witness = measure.gateway_constraints(board, ANTIPODES)
    assert len(constrained) == witness["locally_forced_bridges"]
    assert not witness["empty_background_undecided_components"]
    rng = np.random.default_rng(127 + level)
    values = rng.normal(size=(board.cells, 12)) + 1j * rng.normal(
        size=(board.cells, 12)
    )
    psi = (values @ (np.eye(12) - P)).ravel()
    for arm, sign in (("one_load", 0), ("vacancy", -1), ("excess", 1)):
        h_psi = (measure.operator(board, arm, KAPPA) @ psi).reshape(-1, 12)
        alpha = 1 - 4 * KAPPA * sign
        for i, ends in enumerate(by_owner):
            tangent = sum(P[:, a % 12] * (alpha * psi[a] - psi[b]) / 2 for a, b in ends)
            np.testing.assert_allclose(P @ h_psi[i], tangent, atol=2e-13)


def test_flight_receipt_and_centre_probabilities():
    archived = json.loads((measure.HERE / "results.json").read_text())
    fresh = measure.run()
    assert fresh["protocol_sha256"] == archived["protocol_sha256"]
    for old, new in zip(archived["boards"], fresh["boards"], strict=True):
        assert old["constraints"] == new["constraints"]
        for arm in old["arms"]:
            for past, now in zip(old["arms"][arm], new["arms"][arm], strict=True):
                for key in (
                    "outside_probability",
                    "outside_colour_probability",
                    "centre_dark_probability",
                ):
                    assert past[key] == pytest.approx(now[key], abs=1e-11)
                for detector, target in (
                    ("energy_read", "outside_colour_probability"),
                    ("centre_read", "centre_dark_probability"),
                ):
                    assert sum(past[detector]) == archived["shots"]
                    assert (
                        abs(past[detector][0] / archived["shots"] - now[target])
                        < archived["simultaneous_99_percent_probability_radius"]
                    )
                assert (
                    abs(
                        now["centre_dark_probability"]
                        - 0.75 * now["outside_colour_probability"]
                    )
                    < CAL["compiled_centre_error"]
                )
        for direct, changed in zip(
            old["arms"]["one_load"], old["arms"]["exchange_only"], strict=True
        ):
            assert direct["outside_probability"] == pytest.approx(
                changed["outside_probability"], abs=1e-12
            )


def test_joint_filters_preserve_trace_and_reference():
    old = json.loads((measure.HERE / "orientation.json").read_text())
    fresh = orientation.run()
    assert (
        old["protocol_sha256"]
        == hashlib.sha256(
            (measure.HERE / "ORIENTATION_PROTOCOL.md").read_bytes()
        ).hexdigest()
    )
    for past, now in zip(old["reads"], fresh["reads"], strict=True):
        for key in ("launch_acceptance", "accepted", "present", "exact_negative_pt"):
            assert past[key] == pytest.approx(now[key], abs=1e-11)
        assert (
            0 <= now["accepted"] <= now["present"] <= now["launch_acceptance"] + 1e-12
        )
        for prep in past["counts"]:
            for setting in prep:
                assert sum(setting) == old["shots_per_setting"]
        if now["launch"] == now["receive"] == "singlet":
            assert past["measured_negative_pt"] > old["negativity_threshold"]
    for start in range(0, len(old["reads"]), 3):
        whole, singlet, colour = old["reads"][start : start + 3]
        assert whole["present"] == pytest.approx(singlet["present"] + colour["present"])
