import json
from pathlib import Path

import numpy as np

from . import rotation_domain as domain


def test_received_mask_on_the_raw_native_terms_and_all_spatial_turns():
    path = Path(__file__).with_name("rotation_domain.json")
    record = domain.run()
    assert record == json.loads(path.read_text())
    # Independently project each complete numerical native term between all
    # distinct head sectors, including blank and written zero contrast.
    labels = np.repeat([-1, *(abs(d) for _, d in domain.LAW.LABELS)], 9)
    for h in domain.LAW.local_term_parts():
        for left in (-1, 0, 1, 2):
            for right in (-1, 0, 1, 2):
                if left != right:
                    np.testing.assert_allclose(
                        h[np.ix_(labels == left, labels == right)], 0, atol=1e-13
                    )
    receiver = record["receiver"]
    live = {tuple(e) for e in receiver["absolute_one_edges"]}
    assert live == {(157, 161), (158, 166), (159, 165)}
    for row in receiver["proper_rotations"]:
        r = row["permutation"]
        moved = {tuple(sorted(156 + r[x - 156] for x in e)) for e in live}
        assert (moved == live) == (row["order"] == 1)
    assert receiver["half_turn_overlap_histogram"] == {"0": 10, "1": 3, "2": 2}
