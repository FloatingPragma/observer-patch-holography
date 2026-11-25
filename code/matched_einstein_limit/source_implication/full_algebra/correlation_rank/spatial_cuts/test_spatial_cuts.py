import json

import numpy as np
import pytest

from . import held, measure


def test_material_injectivity_and_aperture_control():
    measure.iv.dps = 36
    data = json.loads((measure.HERE / "effects.json").read_text())
    report = measure.certify(data)
    assert report["rank"] == 18 and report["aperture"]["kernel_dimension"] == 7
    assert report["condition_lower_integer"] > 40_000_000
    assert report["modular_error_upper"] / data["scale"] < 1.59e-5
    # A duplicate replaces an independent oriented reading. It must lose rank.
    cuts = {tuple(row["channel"]): row for row in data["cuts"]}
    cuts[measure.CHANNELS[1]]["row"] = cuts[measure.CHANNELS[0]]["row"]
    with pytest.raises((ValueError, np.linalg.LinAlgError), match="rank|[Ss]ingular"):
        measure.certify(data)


def test_held_spectrum_against_separate_flight():
    result = json.loads((measure.HERE / "results.json").read_text())
    measured = held.run(result["spectrum_weights"])
    assert (
        abs(measured["material_entropy_slope"] - measured["spectrum_entropy_slope"])
        < 5e-8
    )
    assert abs(measured["spectrum_entropy_slope"] - 0.00159405) < 1e-7


def test_complete_reader_changes_both_marginals_and_parity():
    from . import full_source

    measure.iv.dps = 36
    source = json.loads((measure.HERE / "local_effects.json").read_text())
    report = full_source.certify(source["rows"])
    assert (
        report["trace_one_rank"] == 31
        and report["survey_volume_tangent_dimension"] == 30
    )
    assert report["full_modular_error_upper"] / measure.survey.SCALE < 3e-5
    arm = full_source.held_volume_arms(report)
    assert abs(arm["combined"] - arm["direct_entropy"]) < 5e-8
    assert abs(arm["ordinary_calibrated_area"] - arm["direct_entropy"]) < 5e-8
    assert min(abs(arm["marginal_and_parity_slopes"][k]) for k in [0, 11, 12]) > 0.05
