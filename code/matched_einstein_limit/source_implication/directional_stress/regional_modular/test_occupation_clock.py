import json

from . import occupation_clock as clock


def test_clocks_built_from_the_two_running_loads():
    fresh = clock.run()
    saved = json.loads(clock.Path(clock.__file__).with_suffix(".json").read_text())
    # The material inverse is a proposed tangent. Certify the observable
    # consequences of both independently rounded proposals.
    for record in (fresh, saved):
        assert record["held_probability_and_material_rows"] == 26
        dwell = record["dwell_clock"]
        assert clock.parse_interval(dwell["clock_speed"]).a > 0
        assert clock.parse_interval(dwell["clock_acceleration"]).a > 0
        for result in dwell["source_controls"].values():
            derivative = clock.parse_interval(
                result["odd_focusing_derivative_in_dwell_time"]
            )
            assert derivative.b < -700000
        controls = record["fixed_occupation_and_rate_controls"]
        assert [a["target_head_coordinate"] for a in controls] == [9, 4, 20, 8]
        for control in controls:
            assert (
                clock.parse_interval(control["augmented_right_inverse_defect_upper"]).b
                < 0.01
            )
        h0, h1 = [
            clock.parse_interval(a["odd_expansion_derivative"]) for a in controls[2:]
        ]
        g0, g1 = [
            clock.parse_interval(a["odd_focusing_derivative"]) for a in controls[2:]
        ]
        det = h0 * g1 - h1 * g0
        assert det.a > 47 and det.b < 70
        claimed = clock.parse_interval(record["nonzero_two_control_minors"]["2,3"])
        assert claimed.a > 47 and claimed.b < 70
