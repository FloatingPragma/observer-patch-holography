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
