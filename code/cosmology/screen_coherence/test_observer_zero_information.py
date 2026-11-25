import json
from pathlib import Path

import numpy as np

from observer_zero_information import (
    channel,
    clock_read,
    frame_choices,
    local_gram,
    midpoint_accelerations,
    send_word,
)


def test_written_information_and_changed_readout():
    word = [1, 0, 1, 1, 0, 0] * 4
    for sign in (-2, 2):
        for direction in ("left", "right"):
            got, success = send_word(word, sign, direction)
            assert got == word and abs(success - 1) < 1e-11
            for row in channel(sign, direction):
                x, r = row["challenge"], row["remote_choice"]
                assert abs(row["joint"][x][x ^ r] - 1) < 1e-12
            # A prior sharp coin read leaves an exactly ambiguous quiet outcome
            # for an outbound zero with a remotely chosen one.
            row = channel(sign, direction, coin=0)[1]
            np.testing.assert_allclose(row["joint"], [[0.5, 0.5], [0, 0]], atol=1e-12)
            assert channel(sign, direction, kappa=0)[1]["joint"][0][1] < 0.01


def test_intrinsic_local_frame_and_clock():
    g = local_gram()
    np.testing.assert_allclose(g @ g, 4 * g, atol=1e-12)
    assert np.linalg.matrix_rank(g, tol=1e-9) == 3
    for i in range(12):
        maps = frame_choices(g, i, (i + 1) % 12)
        assert (
            maps["maps_given_one_direction"]
            == len(maps["possible_second_near_directions"])
            == 5
        )
    assert clock_read([0, 1, 0], 1, 2) == 2.5


def test_radar_midpoint_with_exact_synthetic_motion():
    bank = dict(
        origin_emission_ticks=[0, 4, 8],
        reflectors=[[1, 0], [2, 0], [3, 0]],
        echo_waits=[[4, 4, 4], [6, 4, 6], [8, 6, 6]],
    )
    rows = midpoint_accelerations(bank)
    assert [r["sign"] for r in rows] == [0, 1, -1]
    assert rows[0]["middle_radar_time"] == 7
    # Recover the complete deposited classification from the primary echo bank.
    folder = Path(__file__).parent
    bank = json.loads((folder / "observer_zero_capacity.json").read_text())
    report = json.loads((folder / "observer_zero_information.json").read_text())
    assert midpoint_accelerations(bank) == report["midpoint_radar_accelerations"]
