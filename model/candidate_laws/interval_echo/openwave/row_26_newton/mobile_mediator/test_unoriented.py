"""Check the unbiased path coefficient directly from the clamped native graph."""

from fractions import Fraction
from itertools import product
import json

import numpy as np
from scipy import sparse

from . import matched_energy, measure, unoriented


def test_unoriented_coupling_at_equal_energy():
    board = measure.federation(0)
    raw = json.loads((measure.HERE / "unoriented.json").read_text())
    ordered = [p for p in measure.sources(board) if p[0] >= 2]
    original = json.loads((measure.HERE / "matched_energy.json").read_text())["panels"]
    adjacency = board.adjacency().tocsr()
    for panel, (distance, _, portals) in zip(raw["panels"], ordered, strict=True):
        assert panel["distance"] == distance
        a = matched_energy.regular_pair(board, portals[0])
        b = matched_energy.regular_pair(board, portals[2])
        for pair in (a, b):
            assert [
                matched_energy.source_energy(board, pair, swap)
                for swap in (False, True)
            ] == [116.5] * 2
        allowed = np.setdiff1d(np.arange(board.n), (*a, *b))
        matrix = sparse.csr_matrix(adjacency[allowed][:, allowed], dtype=int)
        ma, mb = [
            np.asarray(adjacency[allowed][:, pins].sum(axis=1)).ravel().astype(int)
            for pins in (list(a), list(b))
        ]
        support = np.flatnonzero(ma)
        paths = np.eye(len(allowed), dtype=int)[:, support]
        for _ in range(len(allowed)):
            weights = mb[:, None] * ma[support]
            if np.any(paths * weights):
                break
            paths = matrix @ paths
        else:
            raise AssertionError("no path reaches the second source")
        exact = Fraction(9, 4) * int(np.sum(weights * paths**2))
        assert exact == Fraction(panel["mean_path_weight"]) > 0
        arms = panel["orientations"]
        assert {(v["q"], v["reverse_a"], v["reverse_b"]) for v in arms} == set(
            product((1, -1), (False, True), (False, True))
        )
        assert sum(Fraction(v["path_weight"]) for v in arms) / 8 == exact
        for v in arms:
            if not v["reverse_a"]:
                ref = next(
                    x
                    for x in original
                    if (x["distance"], x["mobile_charge"], x["reversed_receiver"])
                    == (distance, v["q"], v["reverse_b"])
                )
                lo, hi = ref["heat_trace_beta_eighth"]
                assert max(lo, v["heat"][0]) <= min(hi, v["heat"][1])
        assert 0 < panel["mean_heat_trace"][0] <= panel["mean_heat_trace"][1]
    source = [
        dict(distance=p["distance"], group=p["orientations"]) for p in raw["panels"]
    ]
    assert unoriented.summarize(source) == raw
    gains = [entry["gain"] for entry in raw["inverse_radius_calibrations"]]
    assert all(
        current[0] > 10**12 * following[1] > 0
        for current, following in zip(gains, gains[1:])
    )
