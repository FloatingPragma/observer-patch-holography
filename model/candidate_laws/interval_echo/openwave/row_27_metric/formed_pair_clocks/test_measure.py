import json

import numpy as np

from . import measure as m


def test_all_nine_live_clock_columns_and_joint_encounter():
    m.native_clocks()
    board = m.federation(0)
    origin, detected, sites = m.layout(board)
    x, y, adjacent, _, h, _ = m.pair.generators(board)
    idx = m.pair.ordered_index(*detected, board.n)
    assert not adjacent[m.pair.ordered_index(*origin, board.n)] and adjacent[idx]
    ket = np.zeros(len(x))
    ket[idx] = 1
    for name, seam in sites.items():
        base, occupied = m.enrolled_clock(board, x, y, h, seam)
        assert not occupied[idx]
        accessible = np.linalg.norm((base @ ket)[occupied]) ** 2
        assert accessible == (0.25 if name == "first" else 0.0)


def test_acquired_harmonics_reject_shared_random_time():
    record = json.loads((m.HERE / "results.json").read_text())
    panels = record["panels"]
    for row in panels:
        assert (
            row["nonheralds"] + sum(row["accepted_bins"])
            == record["starts_per_setting"]
        )
        assert (
            row["herald_interval"][0]
            <= record["formation_probability"]
            <= row["herald_interval"][1]
        )
    for t in m.DWELLS:
        one, even, odd = [
            r for r in panels if r["mount"] == "first" and r["dwell"] == t
        ]
        for line in (even, odd):
            assert (
                line["wrong_clock_interval"][0]
                > line["clock_line"] ** 2 * one["wrong_clock_interval"][1]
            )
            assert line["pair_contact_interval"][0] > 0.98
    slow, fast = record["unconditioned"]
    assert sum(slow["bins"]) == sum(fast["bins"]) == record["starts_per_setting"]
    assert fast["wrong_clock_interval"][0] > 32**2 * slow["wrong_clock_interval"][1]
    assert record["source_energy"]["incoming"]["mean_above_uniform"][0] > 0
