import itertools
import json

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from model.executable_model import carrier
from ..frames import EPS
from ..open_release.flight_check import neighbourhood
from .measure import HERE, STARTS, TARGET, apparatus, interval_bank


def records():
    return json.loads((HERE / "results.json").read_bytes())


def test_arm_permutation_and_unique_alternating_read():
    sec, gates, arms, sign, n = apparatus()
    assert set(arms[:, 1]) == gates
    h = sec.hamiltonian(8, False) + np.diag(n / 2)
    vals, vec = np.linalg.eigh(h)
    chosen = abs(vals - TARGET) < 1e-7
    assert chosen.sum() == 1
    psi = vec[:, chosen][:, 0]
    assert np.linalg.norm(sign @ psi - 6 * psi) < 2e-10
    for omit in range(4):
        indices = [
            sec.index[tuple(sorted(word))]
            for word in itertools.product(*np.delete(arms, omit, axis=1))
        ]
        epsilon = np.rint(np.sqrt(6) * EPS).astype(int).ravel()
        # Every column of the exact restricted sign projector is alternating.
        restricted = sign[indices]
        np.testing.assert_array_equal(
            6 * restricted, np.outer(epsilon, epsilon @ restricted)
        )
        assert np.linalg.norm(psi[indices]) > 0.009


def test_principal_compression_and_integer_degrees():
    sec, _, _, _, n = apparatus()
    for level in range(3):
        patch, missing = neighbourhood(level)
        home = np.flatnonzero(np.max(patch.tri, axis=1) < 12)
        h = patch.hamiltonian(8) + sparse.diags(missing / 2)
        np.testing.assert_array_equal(
            h[home][:, home].toarray(), sec.hamiltonian(8) + np.diag(n / 2)
        )
        just_out = np.flatnonzero(np.max(patch.tri, axis=1) >= 12)
        b = h[just_out][:, home]
        np.testing.assert_array_equal((b.T @ b).toarray(), np.diag(n / 4))
        survey = records()["open_boards"][str(20 * 4**level)]["walk_survey"]
        car = carrier.federation(level)
        v = (car.degree == 6).astype(int)
        for _ in range(survey["length"]):
            v = np.array([sum(v[q] for q in row) for row in car.neighbours])
        assert v[[1, 3, 10]].tolist() == survey["gate_walk_counts"]


def test_independent_free_propagator_and_paired_colour_controls():
    data = records()
    sec, _, arms, _, n = apparatus()
    val, vec = np.linalg.eigh(sec.hamiltonian(8, False) + np.diag(n / 2))
    source = vec[:, np.argmin(abs(val - TARGET))]
    patch, missing = neighbourhood(0)
    h = patch.hamiltonian(8) + sparse.diags(missing / 2)
    initial = np.zeros(patch.dim, complex)
    home = np.flatnonzero(np.max(patch.tri, axis=1) < 12)
    initial[home] = source
    final = expm_multiply(-1j * h, initial)
    groups = np.delete(arms, 3, axis=1)
    indices = [patch.rank(*sorted(w)) for w in itertools.product(*groups)]
    projected = final[indices].reshape(3, 3, 3)
    row = data["open_boards"]["20"]["native"]["1.0"]["colours"][3]
    assert abs(abs(np.vdot(EPS, projected)) ** 2 - row["singlet"]) < 1e-12
    assert data["all_time_connected_lower_given_source"]["20"] > 0.98
    phys = data["material_half_turn"]["fringes"]
    assert phys["half_turn_all"]["dark_given_triplets"] < 1e-6
    assert phys["half_turn_one"]["dark_given_triplets"] > 0.01
    assert data["controls"]["symmetric_ground"][3]["nonsinglet"] > 0.9


def test_all_attempts_and_adverse_propagation_enclosure():
    data = records()
    for board in data["open_boards"].values():
        for law in ("native", "exchange"):
            for panel in board[law].values():
                bank = panel["full_board_count_bank"]
                assert (
                    sum(bank["certain"].values())
                    + sum(x["count"] for x in bank["unresolved"])
                    == STARTS
                )
                assert bank["certain"]["singlet"] > 10000
    for key in ("guarded_counts", "filter_stage_counts"):
        groups = data[key] if key == "guarded_counts" else [data[key]]
        for tally in groups:
            assert sum(tally) == STARTS
    wide = interval_bank(
        np.random.default_rng(7), dict(nonsinglet=0.01, triplet=0.05), 0.4, (0.08, 0.09)
    )
    assert sum(x["count"] for x in wide["unresolved"]) > 1_000_000


def test_material_number_phase_is_distinct_from_colour_centre():
    data = json.loads((HERE / "number.json").read_bytes())
    for law in ("native", "half_strength"):
        for source in (
            "vacuum",
            "standard",
            "constituent_0",
            "constituent_1",
            "constituent_2",
        ):
            result = data["controls"][law + "/" + source]
            assert [k for k, v in result["covered"].items() if v] == [
                "empty" if source == "vacuum" else "unit"
            ]
            assert all(sum(bank["counts"]) == STARTS for bank in result["rows"])
            pulse = result["rows"][-2]
            assert pulse["turn"] == 1 and pulse["quadrature"] == "X"
            expected_dark = 0 if law == "native" or source == "vacuum" else 1
            assert pulse["counts"][2 + expected_dark] > 0
            assert pulse["counts"][3 - expected_dark] == 0
    assert abs(data["earlier_colour_centre_dark"] - 0.75) < 1e-8
