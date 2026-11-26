"""Check the deposited raw histories through the native observation program."""

import json
from pathlib import Path

import numpy as np

from .observations import run


def test_replay_and_physical_records(tmp_path):
    here = Path(__file__).resolve().parent
    reference = json.loads((here / "observations.json").read_text())
    actual = run(output=tmp_path)
    native_p = actual.pop("native_half_transfer_probability_02_20")
    assert np.isclose(
        native_p,
        reference.pop("native_half_transfer_probability_02_20"),
        rtol=0,
        atol=1e-12,
    )
    assert actual == reference
    with (
        np.load(tmp_path / "patch_records.npz") as new,
        np.load(here / "patch_records.npz") as archived,
    ):
        assert set(new.files) == set(archived.files)
        for key in new.files:
            np.testing.assert_array_equal(new[key], archived[key])
        for board in actual["boards"]:
            prefix = f"L{board['level']}_"
            history = new[prefix + "singleton_history"]
            tag, mobile, primary, hand = history.transpose(2, 0, 1)
            transfers = new[prefix + "transfers_per_service"]
            np.testing.assert_array_equal((np.diff(hand, axis=0) % 5), primary[1:])
            np.testing.assert_array_equal(
                np.count_nonzero(np.diff(mobile, axis=0), axis=1), 2 * transfers
            )
            assert np.all(mobile.sum(axis=1) == board["initial_doubles"])
            assert np.all(tag == tag[0])
            assert np.any(primary[1:] == 0) and np.any(primary[1:] == 1)
            patches = new[prefix + "triple_support"]
            np.testing.assert_array_equal(
                history[-1, patches], new[prefix + "triple_checkpoint"]
            )


def test_contact_receipts_refer_to_each_owners_actual_last_transfer():
    from .contact_history import ContactHistory, interval_summary

    log = ContactHistory(4)
    start = np.array([0, 2])
    end = np.array([1, 3])
    positive = np.array([True, False, True, True])
    log.record(start, end, np.array([True, True]), positive)
    log.record(start, end, np.array([False, False]), positive)
    log.record(np.array([0]), np.array([3]), np.array([True]), positive)
    np.testing.assert_array_equal(log.array(), [
        [0, 1, -1, -1, 1, 0], [2, 3, -1, -1, 1, 1], [0, 3, 0, 1, 3, 3],
    ])
    summary = interval_summary(log.array(), minimum=2, maximum=3)
    assert summary["source_height"] == 2
    assert summary["flat_diamond_expected_fraction"] == dict.fromkeys((1, 2, 3, 4), 1.0)
    assert summary["intervals"] == 2
    assert summary["mean_ordering_fraction"] == 1


def test_independent_owner_service_changes_only_attending_hands_and_pairs(tmp_path):
    actual = run((0,), services=128, output=tmp_path, provenance=True,
                 shuffled_schedule=True, owner_service_probability=0.7)
    board = actual["boards"][0]
    with np.load(tmp_path / "patch_records.npz") as archive:
        active = archive["L0_owner_service"]
        records = archive["L0_singleton_history"]
        mobile, primary, hands = records[:, :, 1:].transpose(2, 0, 1)
        marks = active & primary[1:].astype(bool)
        np.testing.assert_array_equal(np.diff(hands, axis=0) % 5, marks)
        assert board["positive_marks"] == marks.sum()
        assert board["quiet_marks"] == active.sum() - marks.sum()
        assert np.any(primary[1:][~active])  # A held positive rail need not post a fresh mark.
        endpoints = np.diff(mobile, axis=0) != 0
        assert endpoints.any() and not np.any(endpoints & ~active)
        np.testing.assert_array_equal(endpoints.sum(axis=1), 2 * archive["L0_transfers_per_service"])
        lo, hi = board["owner_opportunity_range"]
        assert 0 < lo < hi < 128
        assert [lo, hi] == [active.sum(axis=0).min(), active.sum(axis=0).max()]
        assert len(archive["L0_contact_history"]) == board["transfers"] > 0
