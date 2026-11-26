"""Independent source, shape and blinding checks for the matched annulus diagnostic."""

from fractions import Fraction as F
import json

from . import measure as m


def test_velocity_columns_are_blind_during_pair_choice():
    t1, t2 = [
        (m.DATA / n).read_text().splitlines() for n in ("table1.dat", "table2.dat")
    ]
    masked = [r[:26] + "?" * 12 + r[38:] for r in t2]
    names, rings = m.photometric_sources(t1, t2)
    blinded, ring_blind = m.photometric_sources(t1, masked)
    assert names == blinded
    for cut in m.CUTS.values():
        assert {n: m.select(rings[n], cut) for n in names} == {
            n: m.select(ring_blind[n], cut) for n in names
        }


def test_shape_identity_and_common_distance_and_inclination():
    for a, b in ((1, 2), (3, 2), (4, 4)):
        t = F(a * a, b * b)
        for u in (F(0), F(3, 7), F(11)):
            pred = (t + u * F(a, b)) / (1 + u)
            assert m.ordering(pred, pred, t) == "overlap"
            faint = (t + (u + 2) * F(a, b)) / (1 + u)
            assert m.ordering(faint, faint, t) != "smaller_faint"
    vf, vb, rf, rb = F(149), F(152), F(13), F(5)
    for d, i in ((F(2), F(3, 4)), (F(7, 11), F(16, 13))):
        ratio = ((i * vf) ** 2 / (d * rf)) / ((i * vb) ** 2 / (d * rb))
        assert ratio == (vf / vb) ** 2 * rb / rf


def test_full_receipt_and_two_opposite_shapes():
    d = m.run()
    assert d == json.loads((m.HERE / "results.json").read_text())
    deep = d["cuts"]["deep"]
    assert deep["count_by_number_of_speed_errors"]["2"] == {
        "larger_faint": 2,
        "overlap": 60,
        "smaller_faint": 8,
    }
    for c in d["cuts"].values():
        assert len({x["galaxy"] for x in c["pairs"]}) == c["selected"]
        assert c["eligible_quality_inclination_without_pair"] + c["selected"] == 153
    by = {r["galaxy"]: r for r in deep["pairs"]}
    assert by["D631-7"]["speed_error_boxes"]["2"]["observed_range"][0] > 1.44
    assert by["UGC05764"]["speed_error_boxes"]["2"]["observed_range"][1] < 0.79
