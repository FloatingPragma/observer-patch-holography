from copy import deepcopy
import json
from unittest.mock import patch

import numpy as np
import pytest

from .. import experiment as ex
from . import spatial


@pytest.fixture(scope="module")
def exposed_bank(tmp_path_factory):
    # The first declared pilot is already public. A fresh draw is never made.
    board = json.loads((ex.HERE / "holdout.json").read_text())["boards"][0]
    directory = tmp_path_factory.mktemp("spatial_pilot")

    def keep(stream, key, counts):
        if key.startswith("clock_") or "_scan" in key or key.endswith("_sum"):
            np.save(directory / (key + ".npy"), np.asarray(counts, dtype="<u8"))

    with patch.object(ex.CountStream, "keep", keep):
        ex.spectroscopy(
            board["dimers"], np.random.default_rng([board["seed"], board["level"], 47])
        )
    return board, directory


def test_independent_prediction_and_identical_count_foil(exposed_bank):
    board, directory = exposed_bank
    observed = spatial.score(board, directory)
    assert observed["status"] == "pass"
    assert (
        observed["candidate_q"]
        == observed["even_trained_q"]
        <= observed["quadratic_cutoff"]
    )
    assert observed["full_recurrence_diagnostic_q"] != observed["candidate_q"]
    assert np.allclose(
        observed["omega_from_pair_counts"],
        [r["omega"] for r in board["observer"]["modes"]],
        atol=1e-8,
        rtol=0,
    )
    launch_k = (
        np.pi
        * np.arange(1, 1 + len(board["observer"]["modes"]))
        / (board["dimers"] + 2)
    )
    foil = spatial.score(board, directory, launch_k)
    assert foil["records_complete"] and foil["status"] == "physical_fail"
    assert foil["launch_label_rejected"] and not foil["profile_pass"]
    assert foil["candidate_q"] == observed["launch_label_q"]
    assert foil["even_trained_k"] == observed["even_trained_k"]
    legacy_foil = deepcopy(board)
    for row, label in zip(legacy_foil["observer"]["modes"], launch_k, strict=True):
        row["k"] = float(label)
    unchanged = spatial.score(legacy_foil, directory)
    assert (
        unchanged["status"] == "pass"
        and unchanged["candidate_q"] == observed["candidate_q"]
    )
    assert unchanged["full_recurrence_diagnostic_q"] == observed["launch_label_q"]


def test_missing_invalid_and_physically_quiet_are_distinct(exposed_bank, tmp_path):
    board, directory = exposed_bank
    missing = spatial.score(board, tmp_path)
    assert missing["status"] == "missing_records" and not missing["records_complete"]
    for item in directory.glob("*.npy"):
        (tmp_path / item.name).symlink_to(item)
    scan = tmp_path / "p1_2_scan0.npy"
    scan.unlink()
    np.save(scan, np.array([-1], dtype=int))
    invalid = spatial.score(board, tmp_path)
    assert invalid["status"] == "invalid_records" and not invalid["records_complete"]
    quiet = np.full(
        (spatial.SAMPLES, board["dimers"] + 1), spatial.SIGN_SHOTS // 2, dtype="<u8"
    )
    np.save(scan, quiet)
    complete = spatial.score(board, tmp_path)
    assert complete["status"] == "physical_fail" and complete["records_complete"]
    assert not complete["profile_pass"]
    (tmp_path / "p2_5_scan1.npy").unlink()
    incomplete = spatial.score(board, tmp_path)
    assert (
        incomplete["status"] == "missing_records" and not incomplete["records_complete"]
    )
