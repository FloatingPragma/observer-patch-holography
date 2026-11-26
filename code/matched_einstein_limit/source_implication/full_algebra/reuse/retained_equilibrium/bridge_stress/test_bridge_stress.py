import json

from .measure import HERE, report


def test_material_and_entropy_holds():
    fresh = report()
    assert fresh == json.loads((HERE / "results.json").read_text())
    for bank in ("source_certificate", "simulated_counts"):
        material = fresh[bank]
        assert material["determinant"][0] > 27000
        assert (
            material["vacancy"]["work_over_current"][0]
            > material["excess"]["work_over_current"][1]
        )
    independent = fresh["count_only_material_calibration"]
    assert independent["vacancy_current_at_ratio_endpoints"][0][1] < 0
    assert independent["vacancy_current_at_ratio_endpoints"][1][0] > 0
    assert independent["held_excess_current"][1] < -6000
    constraint = fresh["fixed_vacancy_current_and_volume"]
    assert (
        constraint["excess_current"][1] < -10000 and constraint["entropy"][0] > 115000
    )
    assert fresh["liminf_area_loss_per_step_squared_lower"] >= 680
