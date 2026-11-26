import json
from copy import deepcopy

import pytest

from . import measure


def counts():
    return json.loads((measure.HERE / "records.json").read_text())


def test_replay_and_signed_coverage():
    archive = counts()
    assert measure.acquire() == archive
    contrasts, radius = measure.decode(archive)
    assert 0.016 < radius < 0.017
    for (menu, _, _), mean in contrasts.items():
        assert abs(mean - measure.SIGNED_TARGET[menu]) < radius
        if menu in {"pair", "reversed_polarity"}:
            assert abs(mean) > radius
    # Both probe frames and the independent dark doublet output matter.
    for row in archive["banks"]:
        if row["menu"] == "pair":
            other = 1 - row["level"]
            assert 0.39 < row["counts"][other] / measure.SHOTS < 0.41


def test_decoder_exposes_missing_or_relabeled_data_and_cut_hits():
    archive = counts()
    missing = deepcopy(archive)
    missing["banks"].pop()
    with pytest.raises(ValueError, match="Missing receiver bank"):
        measure.decode(missing)
    wrong_labels = deepcopy(archive)
    wrong_labels["detector_words"][:2] = reversed(wrong_labels["detector_words"][:2])
    with pytest.raises(ValueError, match="calibration"):
        measure.decode(wrong_labels)
    failed = deepcopy(archive)
    failed["banks"][-1]["counts"] = [1, 0, 0, measure.SHOTS - 1]
    with pytest.raises(ValueError, match="Cut transfer"):
        measure.decode(failed)
