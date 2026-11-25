from fractions import Fraction

from .flight_confluence import run


def test_protected_mark_fork_and_uniform_guard():
    d = run()
    assert d["installed_winner"] == [0, 1]
    assert d["isolated_pi_pulse"]["other_winner"] == [1, 82]
    assert d["isolated_pi_pulse"]["output_address_total_variation"] == 1
    assert d["raw_contact_pi_pulse_error"] < 8e-14
    assert d["exchange_off_other_address_probability"] == 0
    lower = float(Fraction(d["all_sphere_short_time_lower"]))
    board = d["native_panels"][0]
    assert board["source_columns"] == 240
    assert board["times"]["1/4"]["all_tested_minimum"][0] > lower
    assert board["times"]["1"]["all_tested_minimum"][0] > 0.61
    for row in d["complete_first_address_counts"]:
        assert sum(row["pre_sweep_occupied_each_seam"]) == d["shots_per_schedule"]
        error = abs(
            row["other_address_count"]
            - d["shots_per_schedule"] * row["other_address_probability"]
        )
        assert (
            error
            < 5
            * (
                d["shots_per_schedule"]
                * row["other_address_probability"]
                * (1 - row["other_address_probability"])
            )
            ** 0.5
        )
    guard = d["random_guard"]
    assert (
        guard["exact_one_offer_probability"] == "1/630"
        and guard["expected_offers"] == 630
    )
    for row in guard["simulated_first_offers"]:
        assert row["hit"] + row["quiet"] == d["shots_per_schedule"]
        mean = d["shots_per_schedule"] / 630
        assert abs(row["hit"] - mean) < 5 * mean**0.5
