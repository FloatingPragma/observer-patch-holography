"""Independent geometry, compatible alphabet and entropy differential checks."""
from collections import Counter
from decimal import Decimal, localcontext
from fractions import Fraction as F
from itertools import product
import json
import math

import numpy as np
import pytest

from . import measure as m


def archived():
    return json.loads((m.HERE / "results.json").read_text())


@pytest.mark.parametrize("field", ["area", "entropy_slope_over_ell4", "swapped_balance_over_ell4",
                                   "worldlines", "gain"])
def test_replay_rejects_changed_survey_readings(tmp_path, field):
    original = archived()
    changed = archived()
    row = changed["rows"][0]
    if field == "area":
        row["plus"]["area"] = "100000"
    elif field == "worldlines":
        row[field][0][0] += 1e-8
    elif field == "gain":
        changed[field] = 2
    else:
        row[field] = 987.0
    path = tmp_path / "results.json"
    path.write_text(json.dumps(changed))
    with pytest.raises(ValueError, match="Different"):
        m.verify_saved(original, path)


def test_shared_worldlines_and_nested_restrictions():
    saved = archived()
    rows = saved["rows"]
    whole = m.memberships(rows)
    assert len(rows) == 21 and len(saved["frames"]) == 10
    for i, (r, internal) in enumerate(zip(rows, whole, strict=True)):
        assert internal == r["whole_bond_indices"]
        assert i not in internal
        assert r["relative_radial_clearance"] > .11
        for j, s in enumerate(rows):
            if r["shape"] == s["shape"] and s["inverse_ell"] > r["inverse_ell"]:
                assert j in internal and set(whole[j]) <= set(internal)
    for axis in "xyz":
        assert saved["frames"]["+"+axis] == saved["frames"]["-"+axis] == axis


def test_frame_and_direct_pullback_volume_and_area(monkeypatch):
    repo = m.HERE.parents[3]
    monkeypatch.syspath_prepend(str(repo/"code"))
    from matched_einstein_limit import frame
    from matched_einstein_limit.source_implication.full_geometry.station_balls.timelike_areas.measure import Slice

    look = dict(((r["inverse_ell"], r["shape"]), r) for r in archived()["rows"])
    for label, native in zip(m.FRAMES, frame.FRAME, strict=True):
        gamma, *spatial = map(float, native)
        spatial = np.array(spatial)
        g, L = m.lift(m.FRAMES[label])
        assert gamma == g
        np.testing.assert_allclose(L@L, np.eye(3)+np.outer(spatial, spatial), atol=4e-15)
        for d in m.SCALES:
            row = look[d, m.FRAMES[label]]
            for sign, key in ((1, "plus"), (-1, "minus")):
                direct = Slice(gamma, sign, 8)
                direct.lift = L
                direct.dt_squared = np.outer(spatial, spatial)
                r = float(row[key]["radius"])
                assert math.isclose(direct.volume(r), 4*math.pi/(3*d**3), rel_tol=4e-13)
                assert math.isclose(direct.area(r), float(row[key]["area"]), rel_tol=4e-13)


def test_floor_replay_and_held_geometry():
    saved = archived()["rows"]
    replay = m.measure()["rows"]
    for old, new in zip(saved, replay, strict=True):
        d = old["inverse_ell"]
        assert old["M"] == d**5
        assert old["plus_d"] == new["plus_d"] and old["minus_d"] == new["minus_d"]
        assert abs(old["balance_over_ell4"]-new["balance_over_ell4"]) < 1e-9
        assert new["M_times_max_area_order_change"] < new["min_word_floor_margin"]/100
        assert new["max_relative_volume_residual"] < 1e-16
        N = old["legal_pairs"]
        with localcontext() as ctx:
            ctx.prec = 50
            for key in ("plus", "minus"):
                count = old[key+"_d"]
                assert N % count == 0 and N >= count**2
                assert int(old["M"] * Decimal(old[key]["area"]).exp()) == count
        assert abs(old["balance_over_ell4"]) < 4/d
        assert -0.694 < old["biased_reference_residual_nats"] < -0.692
        if old["shape"] in ("rest", "z", "xy"):
            assert abs(old["area_slope_over_ell4"]) < 1e-10
        else:
            assert abs(old["area_slope_over_ell4"]) > 8
            assert abs(old["swapped_balance_over_ell4"]) > 17
            ratio = old["homogeneous_line_area_slope_over_ell4"]/old["area_slope_over_ell4"]
            assert -.668 < ratio < -.666


def test_entropy_derivative_for_a_global_uneven_microstate_tangent():
    # Small exact enumerator with a partially read bond and a wholly internal bond.
    # The signed tangent changes conditional microscopic frequencies as well as g.
    hist, tangent, dg = Counter(), Counter(), Counter()
    count = {0: 2, 1: 3}
    N, internal_N = 12, 6
    for tag, rank, inner in product(range(2), range(N), range(internal_N)):
        left, right = m.word_from_rank(rank, count[tag], N)
        inside = m.word_from_rank(inner, 2, internal_N)
        record = tag, left, inside
        hist[record] += F(1, 2*N*internal_N)
        delta = F(2*tag-1, 2*N*internal_N) + F(2*right-(N//count[tag]-1), 5*N*internal_N)
        tangent[record] += delta
        dg[tag] += delta
    assert sum(tangent.values()) == 0 and dg == {0: -F(1, 2), 1: F(1, 2)}
    derivative = -sum(float(tangent[k])*math.log(float(p)) for k, p in hist.items())
    assert math.isclose(derivative, -float(dg[1])*math.log(3/2), abs_tol=1e-13)
    for tag in count:
        pairs = [m.word_from_rank(i, count[tag], N) for i in range(N)]
        assert len(set(pairs)) == N
        assert set(Counter(i for i, _ in pairs).values()) == {count[tag]}
        assert set(Counter(j for _, j in pairs).values()) == {count[tag]}
