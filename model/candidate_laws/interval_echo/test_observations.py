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
