import json
from pathlib import Path

from . import measure


def test_certified_sector_separation_and_finite_receiver():
    receipt = measure.run()
    saved = json.loads(Path(measure.__file__).with_name("results.json").read_text())
    assert receipt["joint_lie_dimension"] == 23
    assert receipt["central_weights"] == [1,2]
    assert receipt["simply_connected_semisimple_times_circle_kernel_order"] == 12
    single, pair = receipt["sectors"]
    assert [s["port_zero_traceless_squared_HS"] for s in receipt["sectors"]] == [
        s["port_zero_traceless_squared_HS"] for s in saved["sectors"]]
    word = single["common_turn"]["2048"]
    assert word["rigorous_operator_error_upper"] < .0121
    assert word["0"]["total_variation"] > 1-1e-11
    assert abs(word["16"]["ideal_total_variation"]-.344050982153)<1e-10
    assert word["zero_gain_operator_deviation"] < 1e-5
    assert sum(word["counts"]["pulse"][1:3]) == 0
    assert sum(word["counts"]["sham"][1:3]) == sum(word["counts"]["occupation_gain_zero"][1:3]) == 8192
    pair_word = pair["common_turn"]["2048"]
    assert pair_word["0"]["total_variation"] < 1e-5
    assert pair_word["rigorous_operator_error_upper"] < .0024
    assert all(len(counts)==31 and sum(counts)==8192 for counts in pair_word["counts"].values())
