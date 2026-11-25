import json
from pathlib import Path
import numpy as np
import pytest

from . import complex_effects as ce


def test_source_effects_and_all_complete_counts():
    data = ce.report()
    saved = json.loads(Path(ce.__file__).with_suffix(".json").read_text())
    assert data["simulated_counts"] == saved["simulated_counts"]
    assert data["effect_operator_error"] < 2e-12
    assert data["commutator_norm"] > 0.01
    assert (
        abs(sum(data["disjoint_probabilities"][:2]) - data["disjoint_probabilities"][2])
        < 1e-14
    )


def test_rank_defects_and_incompatible_pair():
    identity, zero = np.eye(4), np.zeros((4, 4))
    sharp = np.diag([0, 1, 0, 1])
    for a, b in (
        (zero, identity),
        (identity, zero),
        (sharp, identity - sharp),
        (zero, zero),
    ):
        v = ce.three_effects(a, b)
        columns, _ = ce.compile_isometry(v.astype(complex), ce.native_block())
        np.testing.assert_allclose(columns, v, atol=1e-12)
    with pytest.raises(ValueError):
        ce.three_effects(identity, sharp)
