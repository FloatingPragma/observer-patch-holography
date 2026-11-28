import json
import numpy as np
from cosmology.rar_deep_regime.compensation.trail_profile import measure


def test_discrete_envelope_and_common_mass_scaling():
    boards = json.loads((measure.native.HERE / "results.json").read_text())["boards"]
    middle = boards[1]
    for scale in (1.0, 0.1, 10.0):
        results = measure.envelope(
            62.71 * scale, 291.07 * scale, middle, "geodesic", 1.0
        )
        assert 0.8549 < results["all_free_times"]["speed_ratio_ceiling"] < 0.8551
        assert 0.8052 < results["at_release"]["speed_ratio_ceiling"] < 0.8054


def test_no_reported_significance_from_independent_errors():
    result = json.loads((measure.HERE / "results.json").read_text())
    chosen = next(
        c
        for c in result["cases"]
        if c["primary"]
        and c["owners"] == 80
        and c["layout"] == "geodesic"
        and c["inner_baryon_fraction"] == 1.0
    )
    assert not np.any(chosen["comparison"]["above_bonferroni_four_one_sided_95"])
