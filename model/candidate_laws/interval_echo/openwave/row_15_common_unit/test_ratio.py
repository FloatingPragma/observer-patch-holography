import gzip
import json

import numpy as np

from ..electron_field import replay
from . import measure


def test_calibration_bound_is_sharp_for_both_possible_orderings():
    for ratio in (0.6, 1.0, measure.TARGET, 1.2, 2.0):
        r = ratio / measure.TARGET
        a = 2 / (1 + r)
        best = measure.scale_score([ratio, ratio])[
            "best_common_unit_relative_error_at_least"
        ]
        assert np.isclose(max(abs(a - 1), abs(r * a - 1)), best)
        for factor in np.linspace(0.5, 1.5, 51):
            assert max(abs(factor * a - 1), abs(factor * r * a - 1)) + 1e-14 >= best


def test_zero_kick_matches_direct_double_commutator():
    import numpy as np
    from ..row_14_profile import measure as source

    words, swaps = source.sector(0, 1)
    coordinates = (words != 0) @ measure.impulse.vertices()
    rng = np.random.default_rng(889)
    amp = rng.normal(size=len(words))
    amp /= np.linalg.norm(amp)
    rho = np.outer(amp, amp)
    h, _, _ = source.operators(words, swaps, np.ones(12, dtype=bool))
    direct = 0.0
    for x in coordinates.T:
        delta = x[:, None] - x
        direct -= np.sum(h * rho * delta**2) / 3
    a, b = source.EDGES[0]
    ell2 = ((measure.impulse.vertices()[a] - measure.impulse.vertices()[b]) ** 2).sum()
    bank = measure.zero_kick(rng, words, swaps, rho, 1.0, 1, ell2)
    lo, hi = measure.interval(bank)
    assert lo < direct < hi


def test_complete_records_and_previous_held_impulse_comparison():
    receipt = json.loads(
        gzip.decompress((measure.HERE / "counts.json.gz").read_bytes())
    )
    earlier = json.loads(
        gzip.decompress((measure.impulse.HERE / "counts.json.gz").read_bytes())
    )
    for tag, panel in receipt["materials"].items():
        source = panel["sources"]
        # Interval ends are binary64 dot products of the stored counts with their levels;
        # the platform BLAS can differ in the last bit.
        comparison = json.loads(json.dumps(measure.compare(source)))
        replay.compare(comparison, panel["comparison"], f"row_15/{tag}")
        for row in source.values():
            for what in ("energy", "curvature"):
                bank = row[what]
                assert sum(bank["counts"]) + bank["refused"] == measure.STARTS
        if tag == "native":
            for name in ("positive", "neutral"):
                lo, hi = measure.interval(source[name]["curvature"])
                for v in earlier["panels"][name + "/native"]["recoil"].values():
                    assert lo <= v[2] + v[3] and v[2] - v[3] <= hi


def test_all_impulse_shares_exclude_half_contrast_and_bracket_native_fitting():
    from . import momentum

    splits = momentum.possible_splits()
    assert splits["half_contrast"] == []
    ((lo, hi),) = splits["native"]
    assert 0.8085 < lo / momentum.PARTS < hi / momentum.PARTS < 0.8089
    assert (
        measure.scale_score([measure.TARGET - 0.001, measure.TARGET + 0.001])[
            "best_common_unit_relative_error_at_least"
        ]
        == 0
    )
