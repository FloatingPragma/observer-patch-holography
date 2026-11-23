"""External transfer invariants and custody of the archived data."""

import importlib.util
import math
from pathlib import Path


def module(name):
    path = Path(__file__).parent / (name + ".py")
    spec = importlib.util.spec_from_file_location("udg_" + name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def test_sample_inverse_units_and_endpoint_conservatism():
    measure = module("measure")
    archive = measure.load(measure.HERE / "results.json")
    assert measure.report() == archive
    assert len(archive["six_galaxies"]) == 6
    assert archive["misses"] == [6, 6]
    for row in archive["six_galaxies"]:
        assert row["bounds"][1]["minimum_gap_kms"] < row["bounds"][0]["minimum_gap_kms"]
        inverse_from_ratios = archive["acceleration_scale"] / row["ratio"] ** 4
        assert math.isclose(row["effective_a0_m_s2"], inverse_from_ratios)
    source_names = (measure.HERE.parent / "data" / "table1.dat").read_text()
    assert all(str(row["agc"]) not in source_names for row in archive["six_galaxies"])


def test_resolved_disk_prediction_and_inversion():
    measure = module("measure")
    report = measure.report()
    assert report["followup"]["max_abs_baryonic_two_curve_discrepancy_kms"] < 0.06
    for rows in report["resolved_rings"].values():
        assert len(rows) == 5 and all(r["selected"] for r in rows)
        for r in rows:
            radius = r["radius_kpc"] * measure.KPC
            vb = r["baryonic_kms"] * 1000
            delta_v2 = r["predicted_kms"] ** 2 * 1e6 - vb**2
            assert math.isclose(
                delta_v2**2, vb**2 * radius * report["acceleration_scale"]
            )
            observed_scale = (r["observed_kms"] ** 2 * 1e6 - vb**2) ** 2 / (
                radius * vb**2
            )
            assert math.isclose(r["signed_inverse_a0_m_s2"], observed_scale)
    final = report["resolved_rings"]["figure8"][-1]
    assert 26.93 < final["observed_kms"] < 26.95
    assert 62.90 < final["predicted_kms"] < 62.92
    assert 18.1 < final["scale_floor_residual"] < 18.3
