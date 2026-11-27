import json
from itertools import permutations

from .measure import HERE, bands, decode, strongest


def test_lossy_fixed_detector_and_moving_control():
    # Effects diagonal in the detection basis; all three outcomes can be missed.
    effects = [(1, 3, 6), (8, 3, 1), (1, 2, 3)]
    stocks = (300, 900, 1500)
    arrival = (1, 3, 2)
    groups = []
    for order in permutations(range(3)):
        clicks = [
            [stocks[i] * arrival[j] * count for count in effects[j]]
            for i, j in enumerate(order)
        ]
        groups.append([bands(row) for row in clicks])
    assert all(strongest(groups[0], other)["gap"] <= 0 for other in groups)
    before = [bands([20000, 30000, 50000])] * 3
    after = [bands([10000, 10000, 80000])] * 3
    assert strongest(before, after)["gap"] > 0.27
    assert strongest([bands([0, 0, 0])] * 3, after)["gap"] <= 0


def test_registered_counts_and_source_mixture_exclusion():
    fit = decode(json.loads((HERE / "counts.json").read_text()))
    assert fit["resolved"] == {"native": 11, "which_head": 19, "contrast_four": 1}
    trial = next(
        r
        for r in fit["comparisons"]
        if (r["mode"], r["owner"], r["panel"], r["times"])
        == ("native", 17, 1, [320, 640])
    )
    assert trial["gap"] > 0.0567
    assert trial["weights"] == (0, 0, 1)
