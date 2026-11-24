"""Independent station-log and time-dependent route checks."""

import json
from collections import Counter
import numpy as np
from .experiment import HERE, layout
from .capacity import path_records


def test_all_enrolled_station_markers_and_settled_growth():
    record = json.loads((HERE / "volume.json").read_text())
    assert len(record["panels"]) == 8 and record["selected_blank_tetra"] == "1,4,6"
    assert all(a != b for a, b in zip(record["clock_bits"], record["clock_bits"][1:]))
    fixed = 0
    for p in record["panels"]:
        assert len(p["return_delays"]) == 10 and all(
            len(r) == 240 for r in p["return_delays"]
        )
        assert p["counts"] == [sum(e >= 0 for e in r) for r in p["return_delays"]]
        assert all(b <= a for a, b in zip(p["counts"], p["counts"][1:]))
        assert all(len(t["determinants"]) == 165 for t in p["tetrahedra"])
        fixed += all(t["determinants"]["1,4,6"] == 5632 for t in p["tetrahedra"])
    assert fixed == 6
    assert (
        Counter(record["panels"][0]["tetrahedra"][0]["determinants"].values())[
            "non_euclidean"
        ]
        == 135
    )
    v = np.array(record["panels"][0]["growing_aperture"]["marker_counts"])
    assert v.tolist() == [1, 5, 12, 18, 23, 24, 24]
    assert np.diff(v, 2)[0] > 0 and np.diff(np.cbrt(v), 2)[0] < 0


def test_independent_routes_and_exact_volume():
    record = json.loads((HERE / "volume.json").read_text())
    p = record["panels"][0]
    n, edges, *_ = layout(0)
    with np.load(HERE / "records.npz") as bank:
        closed = bank["epoch_280210"]
    assert len(closed) == len(edges) and n == 240
    for k in (0, 1, 7):
        t = record["emissions"][k]
        # Re-run individual point-to-point legs with the independent route search.
        for dest in (1, 4, 6, 144, 178):
            (out,) = path_records(0, closed, 0, t, [dest])
            delay = -1
            if out:
                (back,) = path_records(0, closed, dest, out["arrival"], [0])
                if back and back["arrival"] <= t + 12:
                    delay = back["arrival"] - t
            assert delay == p["return_delays"][k][dest]
        h = p["tetrahedra"][k]
        points = [0, 1, 4, 6]
        rad = np.array([h["echo"][i - 1] for i in points[1:]]) / 2
        g = np.diag(rad**2)
        for i, a in enumerate(points[1:]):
            for j, b in enumerate(points[1 : i + 1]):
                outward = h["tours"][a - 1][b - 1]
                reverse = h["tours"][b - 1][a - 1]
                assert outward == reverse
                leg = outward - rad[i] - rad[j]
                g[i, j] = g[j, i] = (rad[i] ** 2 + rad[j] ** 2 - leg**2) / 2
        np.testing.assert_allclose(np.linalg.det(g), 11, atol=1e-13)
