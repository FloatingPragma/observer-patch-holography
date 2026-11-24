"""Independent linear solves and checks on the exposed port records."""

import json

import numpy as np
from scipy.linalg import solve

from model.executable_model import carrier
from . import experiment as e


def read(name):
    return json.loads((e.HERE / name).read_text())["studies"]


def test_full_board_resolvent_from_native_terms():
    # Assemble each blank contact and the private degree reference separately.
    board = carrier.federation(0)
    off = np.array([[0.5, -0.5], [-0.5, 0.5]])
    seam = carrier.EMITTER_SEAM[0]
    h = np.diag(8 * (6 - board.degree.astype(float)))
    for a, b in board.seams:
        if (a, b) != tuple(seam):
            h[np.ix_([a, b], [a, b])] += off + 8 * np.eye(2)
    port = read("results.json")[0]["ports"]
    z = 1.3819660112500993 + 0.25j
    exact = solve((z + 48) * np.eye(board.n) - h, np.eye(board.n)[:, port])[port]
    with np.load(e.HERE / "results.npz") as raw:
        observed = e.measure(raw["level0_written"], z)[0]
    assert np.max(abs(observed - exact)) < 2 * e.allowance(z)
    radius = read("results.json")[0]["core"]
    counts = np.asarray(radius["counts_by_distance_then_failed_tag"])
    assert (counts.sum(axis=1) == 2**23).all()
    empirical = np.average(
        (0.5 + np.arange(counts.shape[1] - 1)) ** 2, weights=counts[:, :-1].sum(0)
    )
    np.testing.assert_allclose(empirical, radius["radius_links"] ** 2)


def test_same_held_input_discriminates_the_two_interactions():
    for s in read("results.json"):
        fine = s["two_resolutions"][1]
        assert fine["matched_error"] < fine["matched_joint_bound"]
        assert fine["coarse_in_reduced_space_error"] > 100 * fine["coarse_joint_bound"]
        assert fine["bare_in_reduced_space_error"] > fine["bare_joint_bound"]
        for family in ("written", "cut_only"):
            for trial in s[family]:
                for shell in trial["shells"]:
                    assert (
                        shell["constant_interaction_difference"]
                        < shell["constant_interaction_joint_bound"]
                    )
                    assert (
                        shell["held_prediction_error"] < shell["held_joint_error_bound"]
                    )
                    if family == "cut_only":
                        a = shell["vertex_parity"]
                        assert (
                            abs(complex(a["real"][0][0], a["imag"][0][0]))
                            < shell["vertex_operator_error_bound"]
                        )


def test_first_crossing_survives_radius_and_counting_errors():
    cores = {s["cells"]: s["core"] for s in read("results.json")}
    for s in read("cutoffs.json"):
        lo, hi = s["crossing_indices"]
        rows = s["shells"]
        assert rows[lo]["threshold"] == "below" and rows[hi]["threshold"] == "above"
        assert all(
            r["ratio_interval"][0] <= 1.1 and r["ratio_interval"][1] >= 0.9
            for r in rows[:hi]
        )
        assert rows[hi]["ratio_interval"][0] > 1.1
        assert rows[hi]["k"] * cores[s["cells"]]["radius_upper_simultaneous_99"] < 0.5
        post = s["held_checks"][1]
        assert post["matched"]["absolute_error"] < post["matched"]["joint_bound"]
        assert (
            post["frozen_pre_crossing"]["absolute_error"]
            > post["frozen_pre_crossing"]["joint_bound"]
        )
