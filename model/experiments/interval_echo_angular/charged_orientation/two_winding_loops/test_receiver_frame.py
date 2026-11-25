from model.experiments.interval_echo_angular.charged_orientation.two_winding_loops.receiver_frame import (
    inspect,
)


def test_all_receive_frames_and_both_oriented_head_words():
    record = inspect()
    assert record["turns"] == 60
    assert record["fixed_port_turns"] == 5
    assert record["fixed_port_and_face_turns"] == 1
    assert record["proper_half_turns"] == 15 and record["fixed_port_half_turns"] == 0
    assert record["frame_orbit_sizes"] == [60, 60]
    assert (
        record["perpendicular_pairs"] == 15
        and record["compared_frame_orderings"] == 1800
    )
