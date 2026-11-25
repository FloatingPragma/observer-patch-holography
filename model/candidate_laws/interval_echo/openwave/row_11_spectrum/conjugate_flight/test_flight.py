import json
from itertools import combinations

import numpy as np
import pytest
from scipy.linalg import expm

from model.executable_model import cell, federation
from . import measure as ex


def test_native_background_writes_and_cell_rank():
    h, _ = ex.seam_data(8)
    reverse = np.eye(9)[[ex.LAW.SRC.index((2 - a, 2 - b)) for a, b in ex.LAW.SRC]]
    np.testing.assert_allclose(reverse @ h @ reverse, h)
    a = cell().adjacency().toarray()
    for gateways in combinations(range(12), 3):
        bulk = [j for j in range(12) if j not in gateways]
        assert np.linalg.matrix_rank(a[np.ix_(bulk, gateways)]) == 3
    for level in range(3):
        b = federation(level)
        private = 6 - b.degree
        np.testing.assert_array_equal(
            8 * (b.degree + 3 * private) - 8 * (b.degree - private), 32 * private
        )
        assert all(
            sum(b.degree[12 * c : 12 * c + 12] == 6) == 3 for c in range(b.cells)
        )


def test_shutter_and_signed_readout_against_observer_records():
    data = json.loads((ex.HERE / "results.json").read_text())
    board = ex.federation(0)
    ring = data["ring"]
    n = len(ring)
    assert len(set(ring)) == 15 and len(data["owners"]) == 5
    live = set(board.seams) - set(map(tuple, data["cuts"]))
    assert {e for e in live if set(e) & set(ring)} == {
        ex.edge(ring[i - 1], ring[i]) for i in range(n)
    }
    assert all(
        p in ring and board.degree[p] == 5 for p in data["preparation"]["dark_private"]
    )
    assert all(
        b in board.neighbours[a]
        for path in data["preparation"]["phase_routes"]
        for a, b in zip(path, path[1:])
    )
    archive = ex.measure()
    for case, rerun in zip(data["cases"], archive["cases"], strict=True):
        gen = (
            2 * case["kappa"] * np.eye(n)
            + np.eye(n)
            - (np.roll(np.eye(n), 1, axis=0) + np.roll(np.eye(n), -1, axis=0)) / 2
        )
        if not case["private_marks"]:
            gen += np.diag(
                case["kappa"] * ((1 + case["sign"]) ** 2 - 1) * (6 - board.degree[ring])
            )
        for panel in case["sharp_launch"]:
            expected = abs(expm(-1j * panel["time"] * gen)[:, 1]) ** 2
            np.testing.assert_allclose(
                np.array(panel["position"]) / data["shots_per_bank"],
                expected,
                atol=0.0012,
            )
        for r, recomputed in zip(case["banks"], rerun["banks"], strict=True):
            assert len(r["clock_and_position"]) == len(ex.TIMES)
            for panel, again in zip(
                r["clock_and_position"], recomputed["clock_and_position"], strict=True
            ):
                assert sum(panel["position"]) == data["shots_per_bank"]
                np.testing.assert_allclose(
                    np.array(panel["position"]) / data["shots_per_bank"],
                    np.array(again["position"]) / data["shots_per_bank"],
                    atol=0.0012,
                )
            assert r["measured_total_charge"] == case["sign"]
            if case["private_marks"]:
                E = 2 * case["kappa"]
                delta = 1 - np.cos(2 * np.pi * r["mode"] / n)
                for key, x in [
                    ("gap", delta),
                    ("lower_line", E),
                    ("upper_line", E + delta),
                ]:
                    lo, hi = r["simultaneous_intervals"][key]
                    assert lo <= x <= hi
                assert r["inferred_frequency"] == pytest.approx(delta, abs=0.001)
                assert (
                    case["sharp_launch"][-1]["outside_launch_owner"] > 0.69 * ex.SHOTS
                )
            else:
                assert r["inferred_lower_and_upper"] is None
                assert (
                    case["sharp_launch"][-1]["outside_launch_owner"] < 0.01 * ex.SHOTS
                )


def test_conjugation_defect_is_invariant_to_energy_zero_and_recoil_slopes():
    board, graph, ring, cuts, ref, donor = ex.layout()
    n = len(ring)
    H = {}
    for dark in (False, True):
        for sign in (-1, 1):
            setting = ex.prepare_seams(board, graph, ring, cuts, ref, donor, sign, dark)
            H[dark, sign] = ex.ring_hamiltonian(board, ring, setting, 8, sign)[0]
    private = 6 - board.degree[ring]
    np.testing.assert_allclose(
        H[False, 1] - H[False, -1], 32 * np.diag(private), atol=1e-12
    )
    np.testing.assert_allclose(H[True, 1], H[True, -1], atol=1e-12)
    psi = np.ones(n) / np.sqrt(n)
    w = np.mean(private)
    for shift in (0, 3, 32 * w):
        assert np.linalg.norm(
            (H[False, 1] - H[False, -1] - shift * np.eye(n)) @ psi
        ) ** 2 == pytest.approx((shift - 32 * w) ** 2 + 32**2 * w * (1 - w))
    u, v = [2 * (1 - np.cos(2 * np.pi * k / n)) for k in (1, 2)]
    for E in (16.0, 8.0, 0.0):
        assert ((E + v / 2) ** 2 - E**2) / v - (
            (E + u / 2) ** 2 - E**2
        ) / u == pytest.approx((v - u) / 4)


def test_closed_route_rest_kicks_and_three_port_energy_reader():
    from . import rest_kick
    from .pulses import native_word

    saved = json.loads((ex.HERE / "rest_kick.json").read_text())
    pulses = json.loads((ex.HERE / "pulses.json").read_text())
    n = len(saved["ring"])
    for case in saved["cases"]:
        sign = case["sign"]
        # A fresh terminal reference retains its complex amplitudes. Check the
        # read as a linear operator on all three addresses, including a spectator.
        matrices = np.column_stack(
            [
                native_word(v, pulses["bond_minus"], 8, sign, inverse=True)
                for v in np.eye(3, dtype=complex)
            ]
        )
        effect = np.outer(matrices[1].conj(), matrices[1])
        anti = np.array([1, -1, 0]) / np.sqrt(2)
        np.testing.assert_allclose(effect, np.outer(anti, anti), atol=5e-12)
        for r in case["records"]:
            k = 2 * np.pi * r["mode"] / n
            expected = 16 + 1 - np.cos(k)
            assert (
                r["simultaneous_energy_interval"][0]
                <= expected
                <= r["simultaneous_energy_interval"][1]
            )
            assert r["signed_current"] == pytest.approx(np.sin(k), abs=0.0004)
            for key in (
                "position_counts",
                "current_left_right_other",
                "bond_minus_and_other",
            ):
                assert sum(r[key]) == rest_kick.SHOTS
            assert r["word_norm_error"] < 3e-11 and r["free_eigen_residual"] < 3e-11
        assert (
            case["records"][1]["mass_shell_speed_squared_interval"][1]
            < case["records"][2]["mass_shell_speed_squared_interval"][0]
        )
        assert case["changed_law_uniform_overlap"] < 0.02
    actual = rest_kick.run()
    for observed, fresh in zip(saved["cases"], actual["cases"], strict=True):
        np.testing.assert_allclose(
            [r["energy_from_native_bond"] for r in observed["records"]],
            [r["energy_from_native_bond"] for r in fresh["records"]],
            atol=0.0004,
        )
