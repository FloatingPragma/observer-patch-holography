import numpy as np

from .feedback import run


def test_incidence_current_controls_the_same_cage_link():
    result = run()
    assert result["written_cage_receiver_gain"] > 0.1
    for record in result["records"]:
        joint = np.asarray(record["feedback"])
        rejected = np.asarray(record["positive_band_rejections_by_current"])
        assert record["total_variation"] > 1e-4
        np.testing.assert_allclose(joint.sum() + rejected.sum(), 1, atol=2e-10)
