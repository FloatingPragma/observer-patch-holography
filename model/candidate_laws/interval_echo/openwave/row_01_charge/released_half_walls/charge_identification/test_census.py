import json

import numpy as np
import pytest

from . import census


def test_nonadditive_normalization_and_alternating_backgrounds_are_distinct_checks():
    for width in census.WIDTHS:
        test = census.calibrations(14, width)["diagnostics"]
        assert (
            test["normalized_two_end"]
            == test["redistributed_local"]
            == {"failed_total_stocks": 0, "alternating_outputs_twice": [0, 0]}
        )
        assert test["raw_two_end"]["failed_total_stocks"] == 2**width
        assert test["alternating_bulk_response"]["failed_total_stocks"] == 0
        assert set(test["alternating_bulk_response"]["alternating_outputs_twice"]) == {
            -1,
            1,
        }
    with pytest.raises(ValueError, match="padding"):
        census.calibrations(8, 6)


def test_six_site_redistribution_changes_local_readings_and_preserves_each_wall_total():
    n, width = 14, 6
    with np.load(census.DATA, allow_pickle=False) as archive:
        words = archive["words_14"]
        gate, _ = census.masks(words, n, 1)
        codes = words[gate]
    assert len(codes) > 1
    for code in codes:
        bits = [(int(code) >> j) & 1 for j in range(n)]
        walls = [i for i in range(n) if bits[i] == bits[(i + 1) % n]]
        read = census.window_codes(bits, width)
        primitive = np.array([sum(census.digits(w, width)[:2]) - 1 for w in read])
        nonlinear = np.array([int(w % 4 == 3) - int((w >> 1) % 4 == 3) for w in read])
        assert np.any(nonlinear)
        for wall in walls:
            patch = [(wall - j) % n for j in range(width - 1)]
            assert (
                primitive[patch].sum() == (primitive + 3 * nonlinear)[patch].sum() == 1
            )


def test_every_deposited_original_start_cohort_replays():
    assert json.loads(json.dumps(census.decode())) == json.loads(
        (census.HERE / "counts.json").read_text()
    )
