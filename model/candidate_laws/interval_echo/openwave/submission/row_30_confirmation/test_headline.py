"""Challenge the independent arithmetic path with changed forecast and records."""

import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("ow30_independent_headline", HERE / "headline.py")
headline = importlib.util.module_from_spec(spec)
spec.loader.exec_module(headline)


@pytest.fixture
def data():
    return headline.read(HERE / "validation.json"), headline.meter_from_counts(headline.read(headline.METER))


def test_independent_integer_count_meter_and_22_lines(data):
    receipt, meter = data
    result = headline.evaluate(receipt, meter)
    assert result["frequencies"] == 22
    assert result["low_k_lines"] == 14
    assert meter["m"] == pytest.approx(15.000510893879648, abs=1e-10)
    assert meter["c_squared"] == pytest.approx(16.00155913520649, abs=1e-10)
    assert max(b["lattice_max"] for b in result["boards"]) < 0.000575
    assert max(b["continuum_max"] for b in result["boards"]) < 0.003165


@pytest.mark.parametrize("key,factor", [("m", 0.9), ("c_squared", 0.5)])
def test_changed_prediction_rejects_same_measured_records(data, key, factor):
    receipt, meter = data
    with pytest.raises(headline.Rejected, match="lattice"):
        headline.evaluate(receipt, meter | {key: meter[key] * factor})


@pytest.mark.parametrize("key", ["candidate_q", "launch_label_q"])
def test_spatial_summary_gate_has_a_falsifying_outcome(data, key):
    receipt, meter = data
    receipt["boards"][0]["spatial"][key] = 1e10 if key == "candidate_q" else 0
    with pytest.raises(headline.Rejected):
        headline.evaluate(receipt, meter)


def test_inconsistent_even_training_and_temporal_records_invalid(data):
    receipt, meter = data
    receipt["boards"][0]["dispersion"]["modes"][0]["k"] += 0.01
    with pytest.raises(ValueError, match="inputs disagree"):
        headline.evaluate(receipt, meter)


def test_joint_count_budget_failure():
    bench = headline.read(headline.METER)
    bench["all_counts"][0][0][0][0] += 1
    with pytest.raises(ValueError, match="shot budget"):
        headline.meter_from_counts(bench)


@pytest.mark.parametrize("patch", [
    {"calibration_trials": 0},
    {"hand_calibration": [0, 0, 0, 2**22 + 1]},
    {"hand_calibration": [True, 1, 2, 3]},
])
def test_invalid_calibration_is_rejected_before_division(patch):
    bench = headline.read(headline.METER)
    with pytest.raises(ValueError, match="hand calibration"):
        headline.meter_from_counts(bench | patch)


def test_exit_codes_and_preserved_evidence(data, tmp_path):
    paths = [HERE / "validation.json", HERE / "PACKAGE.json", headline.METER]
    before = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    assert headline.main([]) == 0
    for p in paths[:2]:
        (tmp_path / p.name).write_bytes(p.read_bytes())
    receipt, _ = data
    receipt["boards"][0]["spatial"]["candidate_q"] = 1e10
    (tmp_path / "validation.json").write_text(json.dumps(receipt))
    assert headline.main(["--packet", str(tmp_path)]) == 2
    receipt["draw"]["binding"]["source_commit"] = "0" * 40
    (tmp_path / "validation.json").write_text(json.dumps(receipt))
    assert headline.main(["--packet", str(tmp_path)]) == 4
    assert before == {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
