import json

import numpy as np

from ..apparatus import fixture
from . import measure as m
from .certificate import certify


def test_full_graph_bound_and_coherent_label_control():
    res = certify()
    assert res["twice_spectral_interval"] == [16, 54]
    assert res["head_matrix_units"] == 64**2
    assert res["common_commutant_dimension"] == 1
    # Independent local degree bound before any matrix diagonalization:
    # blank edges give [8,9] per occupied endpoint; a single live magnitude-
    # one seam gives [-8,0]. Each active port has at most one such incidence.
    f = fixture()
    ds = [
        (sum(p in edge for edge in f["blank"]), sum(p in edge for edge in f["marked"]))
        for p in f["ports"]
    ]
    assert min(8 * d - 8 * w for d, w in ds) == 8
    assert max(9 * d for d, _ in ds) == 27
    assert max(w for _, w in ds) == 1
    assert abs(m.positive_controls()["koide"] - 2 / 3) < 1e-12
    assert res["filled_seam_energy"] == {"blank": 0, "silent": 0, "unit_even": 8}


def test_fresh_energy_tag_and_disturbed_return_banks():
    saved = json.loads((m.HERE / "results.json").read_text())
    run = m.run()
    p = run["observer_panel"]
    old = saved["observer_panel"]
    assert p["reached_positive_groups"] == 949
    np.testing.assert_allclose(
        [i["center"] for i in p["centers"]],
        [i["center"] for i in old["centers"]],
        atol=1e-8,
        rtol=0,
    )
    assert min(i["component_weight"] for i in p["centers"]) > 0.9999
    assert (
        min(
            abs(i["before_location_read"] - i["after_location_read"])
            for i in p["centers"]
        )
        > 0.42
    )
    true = np.asarray(p["ideal_position_laws"])
    np.testing.assert_allclose(true.sum(axis=1), 1, atol=2e-11)
    maximum = float(
        np.mean(np.sum(true * np.eye(3)[np.argmax(true, axis=0)].T, axis=1))
    )
    assert abs(maximum - p["ideal_balanced_location_ceiling"]) < 1e-10
    a, b = p["balanced_accuracy_simultaneous99"]
    assert a < maximum < b and a > 0.51 and b < 0.53
    assert all(
        sum(row["assigned_counts"]) + row["preparation_refusals"]
        == p["position_starts"]
        for row in p["classification"]
    )
    for i, row in enumerate(p["second_filter_return"]):
        for j, d in enumerate(row):
            if i != j:
                assert (
                    d["location_return_probability"] > 2e-12
                    and d["idle_return_probability"] < 1e-16
                )
                assert d["location_hits"] > 30 and d["idle_hits"] == 0
    assert run["external_targets"]["lab_koide_interval"][0] > 3 / 8
