import gzip
import json
from copy import deepcopy

import numpy as np
import pytest

from . import apparatus as a
from . import measure as m


def saved():
    return json.loads(gzip.decompress(m.RECEIPT.read_bytes()))


def test_replay_and_complete_support_classification():
    archive = saved()
    assert m.acquire() == archive
    summary = m.summarize(archive)
    assert summary["banks"] == 160 and summary["starts"] == 10485760
    assert summary["supports_by_number_of_sharp_labels"] == {1: 393, 2: 111, 3: 7}
    partitions = dict(archive["source_partitions"])
    even = [j for j, x in enumerate(a.INPUTS) if sum(x) % 2 == 0]
    mask = sum(1 << j for j in even)
    assert partitions[mask] == [1, (1 << 2) + (1 << 4) + (1 << 6), 1 << 8]
    assert partitions[511] == [511]


def test_continuation_and_complementary_coherences():
    reads, eps = m.decode(saved())
    for a0, b in a.INPUTS:
        s = a0 + b
        for k in (0, 8):
            for t in (1, 2):
                assert m.mismatch(reads[f"wait/{k}/0/{t}/{a0}{b}"]) == (s % 2) * m.SHOTS
    assert m.mismatch(reads["wait/8/1/2/11"]) == m.SHOTS
    assert 0 < m.mismatch(reads["wait/8/1/2/02"]) < 0.001 * m.SHOTS
    assert m.mismatch(reads["wait/0/1/2/02"]) > 0.64 * m.SHOTS
    for turns in (1, 3):
        assert reads[f"phase/{turns}/0"][1, 2] == m.SHOTS
        assert abs(reads[f"phase/{turns}/1"][1, 2] / m.SHOTS - 0.5) < eps
    assert reads["20/parked/02"][1, 2] == m.SHOTS
    assert reads["20/unchanged_loads/02"][1, 2] == m.SHOTS
    targets = [0, 0.5, 1]
    assert all(
        abs(reads[f"clock/{t}"][2, 1] / m.SHOTS - p) <= eps
        for t, p in enumerate(targets)
    )


def test_joint_outcomes_and_missing_acquisition_have_separate_errors():
    data = saved()
    partial = deepcopy(data)
    partial["banks"].pop("20/repair/02")
    with pytest.raises(ValueError, match="Missing or unregistered"):
        m.decode(partial)
    damaged = deepcopy(data)
    # Same single-endpoint distributions on an odd pair, different event correlations.
    counts = np.array(damaged["banks"]["20/repair/01"]["counts"]).reshape(4, 4)
    step = min(counts[0, 1], counts[1, 0])
    counts[0, 1] -= step
    counts[1, 0] -= step
    counts[0, 0] += step
    counts[1, 1] += step
    damaged["banks"]["20/repair/01"]["counts"] = counts.ravel().tolist()
    with pytest.raises(ValueError, match="Native repair forecast failed"):
        m.decode(damaged)
