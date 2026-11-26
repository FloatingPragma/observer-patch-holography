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


def test_mass_rectangle_and_force_tolerance():
    from . import sensitivity as s

    result = s.run()
    assert result == json.loads((m.HERE / "sensitivity.json").read_text())
    stellar = {
        r["galaxy"] for r in result if r["remains_below_through_stellar_rectangle"]
    }
    matched = {
        r["galaxy"]
        for r in result
        if r["remains_below_for_any_ten_percent_force_match"]
    }
    assert len(stellar) == len(matched) == 5
    assert len(stellar | matched) == 7
    robust_eligible = [
        r["galaxy"]
        for r in result
        if r["always_satisfies_deep_cut"]
        and r["always_satisfies_pair_force_match"]
        and r["galaxy"] in stellar
    ]
    assert robust_eligible == ["UGC06667"]
    middle = next(r for r in result if r["galaxy"] == "UGC07125")
    assert [0.6, 0.56] in middle["eligible_corners_reaching_shape_bound"]
    source = (s.HERE.parent / "data/table2.dat").read_text().splitlines()
    for row in result:
        f, h = [
            s.measured(source[row[k] - 1])
            for k in ("faint_table_line", "bright_table_line")
        ]
        high = h[0] * (f[1] + 2 * f[2]) ** 2 / (f[0] * (h[1] - 2 * h[2]) ** 2)
        t = s.force(f[3], F(".5"), F(".7")) / s.force(h[3], F(".5"), F(".7"))
        limit = max(high, high * high)
        e = (t - limit) / (t + limit)
        assert row["minimum_independent_symmetric_fraction"] == float(e)
        assert t * (1 - e) == limit * (1 + e)
        for d in (F(".4"), F(".49"), F(".6")):
            for b in (F(".56"), F(".73"), F(".84")):
                trial = float(s.force(f[3], d, b) / s.force(h[3], d, b))
                lo, hi = row["stellar_rectangle_baryonic_ratio"]
                assert lo <= trial <= hi
