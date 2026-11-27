import json

import numpy as np
import pytest

from ..joint_receivers.blind import decode
from .column_sums import HERE, SOURCE, analyze, decision, read_banks, evaluate


def test_held_columns_agree_with_separate_decoder():
    fresh = analyze()
    published = json.loads((HERE / "columns.json").read_text())
    assert fresh["input_sha256"] == published["input_sha256"]
    independent = decode(json.loads(SOURCE.read_text()))
    for row, saved, other in zip(
        fresh["results"], published["results"], independent, strict=True
    ):
        assert row["decisions"] == other["decision"] == saved["decisions"]
        np.testing.assert_allclose(
            [c["sum"] for c in row["columns"]],
            np.sum(other["conditional_confusion"], axis=0),
            atol=1e-13,
        )
        for c, s in zip(row["columns"], saved["columns"], strict=True):
            np.testing.assert_allclose(
                [c[k] for k in ("sum", "low", "high")],
                [s[k] for k in ("sum", "low", "high")],
                atol=1e-12,
            )
            assert c["low"] > 1 or c["high"] < 1
        assert row["uniform_output_TV_lower"] >= 0.17


def test_lossless_exchange_with_equal_sources_obeys_the_check():
    banks = {}
    for a in range(3):
        for panel in (-1, 1):
            for role in ("training", "validation"):
                row = [300, 300, 300]
                row[a] += 100
                banks[12, 256, a, panel, role] = row
    sample = evaluate(banks, 12, 256)
    assert all(
        abs(c["sum"] - 1) < 1e-13 and c["low"] <= 1 <= c["high"]
        for c in sample["columns"]
    )
    assert sample["uniform_output_TV_lower"] == 0


def test_ties_fractional_training_and_missing_banks():
    assert decision([[2, 3, 5], [4, 6, 10], [1, 1, 8]]) == [0, 0, 2]
    record = json.loads(SOURCE.read_text())
    record["banks"].append(record["banks"][0])
    with pytest.raises(ValueError, match="uplicate"):
        read_banks(record)
