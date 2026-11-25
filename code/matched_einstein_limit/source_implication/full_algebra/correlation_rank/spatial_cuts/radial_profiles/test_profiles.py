import json

import numpy as np
import pytest
from mpmath import iv

from .. import measure as cuts
from . import effects, held, measure, proofs


def test_native_partition_against_resolved_contacts():
    iv.dps = 50
    raw = json.loads((effects.HERE / "effects.json").read_text())
    grid = raw["shells"]
    patch = cuts.native.chosen(cuts.native.Atlas(1))
    atlas = cuts.native.Atlas(1)
    assert grid == [
        sorted(u for u in patch["cells"] if atlas.dist[patch["center"]][u] == d)
        for d in range(6)
    ]
    regional, shells, cells, node_rows, _, _ = measure.readings()
    older = json.loads((cuts.HERE / "effects.json").read_text())
    edges = [r["channel"] for r in older["cuts"]]
    readings = cuts.interval_rows([r["row"] for r in older["cuts"]])
    prior = json.loads((cuts.survey.HERE / "translation.json").read_text())
    assert all(0 in x for x in (regional - cuts.interval_rows(prior["rows"])).flat)
    for sign, q in enumerate((-1, 1)):
        for i, u in enumerate(cells):
            total = sum(
                (
                    row * ((a == u) - (b == u))
                    for (charge, a, b), row in zip(edges, readings)
                    if charge == q
                )
            )
            assert all(0 in x for x in total - node_rows[sign, i])
        # Adversarial sign control on the largest single shell.
        assert any(0 not in x for x in shells[sign, 1] + shells[sign, 1])


def test_common_shapes_and_two_independent_readers():
    iv.dps = 50
    j, s, cells, nodes, a, v = measure.readings()
    cert = proofs.projective(j, s, a, v)
    assert sum(c["certified_boxes"] for c in cert) >= 2
    assert min(c["minimum_entropy_slope_piconats"] for c in cert) > 10**6
    assert max(c["maximum_defect_ppm"] for c in cert) <= 500000
    assert proofs.independent(j, s)["rank"] == 18
    assert proofs.outer_jump(j, s, a, v)["outer_price_jump_lower"] >= 10**21
    for signs in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
        r = proofs.signed_cone(j, s, a, v, signs)
        assert r is not None and r["minimum_entropy_loss_piconats"] > 5 * 10**7
    angular = proofs.angular_profile(j, cells, nodes)
    assert angular["rank"] == 18 and angular["raw_price_l1_upper"] < 3.46e22
    equal = np.vstack((j, v, s.sum(axis=0)[1:]))
    ent, spread = measure.angular_spread(equal, a, cells, nodes)
    assert ent.a > 4.389e-5 and spread >= 10**18
    with pytest.raises((np.linalg.LinAlgError, ValueError)):
        proofs.independent(j, np.array([s[0], s[0]]))


def test_held_equal_volume_pair():
    iv.dps = 40
    j, s, *_ = measure.readings()
    reader = json.loads((effects.HERE / "results.json").read_text())[
        "common_cell_profile"
    ]
    report = held.run(j, s, reader["regional_and_cell_prices"])
    assert abs(report["arm_volumes"][0] - report["arm_volumes"][1]) < 2e-13
    assert 0.000246 < report["spectral_entropy"] < 0.000249
    assert abs(report["spectral_entropy"] - report["ordinary_calibrated_area"]) < 1e-7
    assert abs(report["common_cell_flux_slope"] - report["spectral_entropy"]) < 2e-5
    assert report["max_relative_equal_mix_drift"] < 1e-6
    assert report["max_relative_regional_drift"] < 1e-8
