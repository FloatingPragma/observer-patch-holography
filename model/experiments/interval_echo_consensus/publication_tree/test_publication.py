import json

import numpy as np

from .measure import HERE, SHOTS, run


def test_all_owners_and_the_reference_agree_on_the_same_emission():
    result = json.loads(json.dumps(run()))
    assert result == json.loads((HERE / "results.json").read_text())
    for board in result["boards"]:
        expected = np.diag([0, 1, 1, 2])
        sent = board["panels"]["broadcast"]
        np.testing.assert_array_equal(sent["exact_source_last"], expected)
        assert (
            sent["delivered_everywhere"] == SHOTS
            and sent["starts_with_conflicting_readings"] == 0
        )
        assert len(board["edges"]) == board["owners"] - 1
        assert sent["readers_without_a_word"] == []
        for ref, counts in enumerate(sent["private_reference_source_last"]):
            a = np.asarray(counts)
            assert np.trace(a) == a.sum() and a[2 - ref].sum() == 0


def test_marginals_cannot_certify_a_shared_record():
    result = json.loads((HERE / "results.json").read_text())
    for board in result["boards"]:
        cut = board["panels"]["cut"]
        assert 0 < len(cut["readers_without_a_word"]) < board["owners"]
        shuffled = np.asarray(board["panels"]["independent"]["exact_source_last"])
        np.testing.assert_array_equal(shuffled.sum(0), [0, 4, 4, 8])
        np.testing.assert_array_equal(shuffled.sum(1), shuffled.sum(0))
        assert np.trace(shuffled) * 8 == shuffled.sum() * 3
    for source in (0, 1):
        assert result["private_phase"][f"{source}/unread"]["restore_probability"] == 1
        assert (
            result["private_phase"][f"{source}/locally_read"]["restore_probability"] == 0.5
        )
