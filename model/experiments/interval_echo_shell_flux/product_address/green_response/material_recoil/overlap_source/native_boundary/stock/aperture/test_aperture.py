from fractions import Fraction
import json

import numpy as np
from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    measure,
    flight,
    observation,
)


def test_exact_vacuum_and_all_cumulative_formation_costs():
    rows = measure.run()
    for board in rows["boards"]:
        path = board["geodesic"]["apertures"]
        assert all(p["radius"] == p["inside"] for p in path)
        assert all(p["closed_formation_excess"] == 8.5 * p["inside"] for p in path[:-1])
        assert (
            path[-1]["evacuated_excess"] - path[-1]["closed_formation_excess"]
            == board["extraction_work"]
            == 12.75
        )
        assert board["whole_deficit_interval"] == [10.25, 16.25]
    sam = rows["empty_board"]["sampled"]
    for r in sam:
        assert sum(r["counts"]) == 2**18
        assert abs(r["sampled_modular_increment"] - r["exact_modular_increment"]) < 0.01


def test_distance_bounds_and_full_contact_sensitivity():
    result = json.loads((measure.HERE / "flight.json").read_text())
    board = measure.trail.measure.observations.federation(2)
    path = json.loads((measure.trail.HERE / "trail.json").read_text())["boards"][-1][
        "path"
    ]
    witness = result["federations"][-1]["far"]
    r = flight.propagation_bound(board, path, witness["port"], (4,))
    exact = Fraction(r["bounds"][0]["integer_bound"])
    assert r["distance"] == 35 and 0 < exact < Fraction(1, 10**23)
    rows = result["one_cell"]["flights"]
    assert max(row["trace_distance"] for row in rows if row["time"] == 1) > 0.89


def test_physical_monitoring_and_reporting_conventions():
    r = observation.run()
    assert (
        r["sign_branches_verified"] == 8
        and r["max_ordinary_density_difference"] < 1e-10
    )
    np.testing.assert_allclose([v["energy"] for v in r["physical_observations"]], 8.5)
    flag = r["thinning"][3]
    assert (
        flag["flagged_increment"]
        < r["thinning"][-1]["flagged_increment"]
        < flag["missing_as_blank_increment"]
    )
