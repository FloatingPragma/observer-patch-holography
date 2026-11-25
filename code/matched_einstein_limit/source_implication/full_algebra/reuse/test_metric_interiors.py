from .metric_interiors import verify


def test_same_survey_and_proper_volume_allow_opposite_local_null_responses():
    result = verify()
    assert result["null_ricci_at_center"] == "-source"
    assert result["relative_small_ball_quadratic_coefficient"] == "source/10"
    assert result["proper_volume_offset"] == "0"
