"""Exact companion replay, strict upstream bindings and bounded quotient checks."""

from __future__ import annotations

import copy
import hashlib
import json
import shlex
import sys
from pathlib import Path

import pytest
import yaml

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))

import gauge_orbit_quotient_instances as quotient  # noqa: E402
import kogut_susskind_fiber_rate_instances as instances  # noqa: E402
import verify_z2_finite_transfer_receipt as verifier  # noqa: E402

WORKFLOW = ROOT / ".github/workflows/finite-gauge-transfer.yml"
REQUIRED_TESTS = {
    "code/yang_mills/tests/test_z2_finite_transfer_receipt.py",
    "code/yang_mills/tests/test_z2_finite_transfer_verifier.py",
    "code/yang_mills/tests/test_z2_transfer_spectrum.py",
    "code/yang_mills/tests/test_z2_full_space_controls.py",
    "code/yang_mills/tests/test_z2_kogut_precision.py",
    "code/yang_mills/tests/test_kogut_susskind_fiber_rate_instances.py",
}


def _rehash(receipt: dict) -> None:
    receipt["sha256_of_runs"] = hashlib.sha256(
        json.dumps(receipt["runs"], sort_keys=True).encode("utf-8")
    ).hexdigest()


@pytest.fixture
def source_specimen(tmp_path, monkeypatch):
    # Rebind retained numbers only to exercise the reader's schema and custody
    # logic during isolated tests. This fixture is not a regenerated live grid.
    receipt = verifier.load_receipt(HERE / "receipts/z2_finite_transfer_receipt.json")
    receipt["producer_sha256"] = hashlib.sha256(
        (HERE / "z2_finite_transfer_receipt.py").read_bytes()
    ).hexdigest()
    _rehash(receipt)
    path = tmp_path / "primary.json"
    path.write_bytes((json.dumps(receipt, sort_keys=True) + "\n").encode("utf-8"))
    monkeypatch.setattr(instances, "COMMITTED_RECEIPT", path)
    return receipt, path


def test_committed_exact_instance_receipt_and_lf_bytes():
    fresh = instances.build()
    serialized = (json.dumps(fresh, indent=2, sort_keys=True, allow_nan=False) + "\n").encode("utf-8")
    assert verifier.load_receipt(instances.OUTPUT) == fresh
    assert instances.OUTPUT.read_bytes() == serialized


def test_companion_binds_the_same_validated_source_bytes(source_specimen):
    receipt, path = source_specimen
    checked = instances.committed_cross_check()
    assert checked["committed_receipt_sha256_at_read"] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert len(checked["kogut_susskind_runs"]) == 6
    # All grid rows are present and validated before conventions are extracted.
    assert len(receipt["runs"]) == 18


def test_companion_writer_is_canonical_lf(source_specimen, tmp_path, monkeypatch):
    output = tmp_path / "instances.json"
    monkeypatch.setattr(instances, "OUTPUT", output)
    instances.main()
    expected = json.dumps(instances.build(), indent=2, sort_keys=True, allow_nan=False) + "\n"
    assert output.read_bytes() == expected.encode("utf-8")
    assert b"\r\n" not in output.read_bytes()


@pytest.mark.parametrize("mutation", [
    "stale_source_hash", "stale_runs_hash", "missing_grid_row", "extra_grid_row",
    "wrong_grid", "unknown_field", "bool_number", "nan_gap", "negative_minimum",
    "false_positive_bound", "duplicate_json_key", "free_zero_forgery", "ks_nonzero_outside",
])
def test_companion_rejects_invalid_primary_before_extracting_conventions(source_specimen, mutation):
    original, path = source_specimen
    changed = copy.deepcopy(original)
    local = next(row for row in changed["runs"] if row["transfer"] == "kogut_susskind")
    if mutation == "stale_source_hash":
        changed["producer_sha256"] = "0" * 64
    elif mutation == "stale_runs_hash":
        changed["runs"][1]["spectral"]["gap_H"] *= 1 + 1e-10
    elif mutation == "missing_grid_row":
        changed["runs"].pop()
    elif mutation == "extra_grid_row":
        changed["runs"].append(copy.deepcopy(local))
    elif mutation == "wrong_grid":
        changed["grid_scope"]["L"] = [2]
    elif mutation == "unknown_field":
        local["certified"] = True
    elif mutation == "bool_number":
        local["spectral"]["gap_H"] = True
    elif mutation == "nan_gap":
        local["spectral"]["gap_H"] = float("nan")
    elif mutation == "negative_minimum":
        local["fiber_dependent_rates"]["rate_min"] = -123.0
    elif mutation == "false_positive_bound":
        fiber = local["fiber_dependent_rates"]
        fiber["rate_min"] = local["variable_rate_floor"]["lower_bound_value"] / 2
        fiber["spread_max_over_min"] = fiber["rate_max"] / fiber["rate_min"]
        assert local["variable_rate_floor"]["numerical_min_respects_bound"] is True
    elif mutation == "free_zero_forgery":
        free = changed["runs"][0]
        free["constant_rate_fit"]["relative_frobenius_residual"] = 0.5
        free["dobrushin"]["eta_star"] = 0.5
        free["dobrushin"]["unit_rate_floor_c_star_times_1_minus_eta"] = 0.5
        free["fiber_dependent_rates"]["offdiagonal_mass_outside_single_flip"] = 100.0
    elif mutation == "ks_nonzero_outside":
        local["fiber_dependent_rates"]["offdiagonal_mass_outside_single_flip"] = 5e-324
    if mutation != "stale_runs_hash":
        _rehash(changed)
    serialized = json.dumps(changed, sort_keys=True)
    if mutation == "duplicate_json_key":
        serialized = '{"schema":"duplicate",' + serialized[1:]
    path.write_bytes((serialized + "\n").encode("utf-8"))
    with pytest.raises(verifier.ReceiptValidationError):
        instances.committed_cross_check()


def test_companion_reads_primary_once(source_specimen, monkeypatch):
    _, path = source_specimen
    original_read = Path.read_bytes
    source_reads = []

    def read_bytes(candidate):
        if candidate == path:
            source_reads.append(candidate)
        return original_read(candidate)

    monkeypatch.setattr(Path, "read_bytes", read_bytes)
    instances.committed_cross_check()
    assert source_reads == [path]


@pytest.mark.parametrize("token,accepted", [("0", True), ("-0.0", True), ("0e-1000", True),
                                            ("1e-1000", False), ("-1e-1000", False),
                                            ("5e-324", False)])
def test_companion_checks_original_number_before_exact_zero(source_specimen, token, accepted):
    receipt, path = source_specimen
    local = next(row for row in receipt["runs"] if row["transfer"] == "kogut_susskind")
    local["fiber_dependent_rates"]["offdiagonal_mass_outside_single_flip"] = json.loads(token)
    _rehash(receipt)
    local["fiber_dependent_rates"]["offdiagonal_mass_outside_single_flip"] = "RAW_NUMBER_TOKEN"
    path.write_text(json.dumps(receipt, sort_keys=True).replace('"RAW_NUMBER_TOKEN"', token),
                    encoding="utf-8")
    if accepted:
        checked = instances.committed_cross_check()
        assert len(checked["kogut_susskind_runs"]) == 6
        assert checked["committed_receipt_sha256_at_read"] == hashlib.sha256(path.read_bytes()).hexdigest()
    else:
        with pytest.raises(verifier.ReceiptValidationError):
            instances.committed_cross_check()


@pytest.mark.parametrize("mutation", ["probability_extrema", "fit_ceiling", "influence_ceiling",
                                      "heat_bath_ceiling", "free_mass"])
def test_companion_rejects_necessary_scalar_contradictions(source_specimen, mutation):
    receipt, path = source_specimen
    run = receipt["runs"][1]
    if mutation == "probability_extrema":
        run["spectral"].update(pi_min=0.03, pi_max=0.5)
    elif mutation == "fit_ceiling":
        run["constant_rate_fit"]["relative_frobenius_residual"] = 2.0
    elif mutation == "influence_ceiling":
        run["dobrushin"].update(eta_star=100.0, dobrushin_condition_holds=False,
                                unit_rate_floor_c_star_times_1_minus_eta=0.0)
    elif mutation == "heat_bath_ceiling":
        run["spectral"]["gap_unit_rate_heat_bath"] = 100.0
    else:
        receipt["runs"][0]["fiber_dependent_rates"]["offdiagonal_mass_single_flip"] *= 0.5
    _rehash(receipt)
    path.write_text(json.dumps(receipt, sort_keys=True), encoding="utf-8")
    with pytest.raises(verifier.ReceiptValidationError):
        instances.committed_cross_check()


@pytest.mark.parametrize("rates", [(2.0, 2.0), (2.5, 10.0 / 3.0)])
def test_abstract_quotient_preserves_scoped_comparison(rates):
    result = quotient.build_instance(*rates)
    assert result["quotient_gap"] == pytest.approx(sum(rates), rel=1e-14, abs=0)
    assert result["product_gap"] == pytest.approx(min(rates), rel=1e-14, abs=0)
    assert result["quotient_gap"] > result["product_gap"]


def test_quotient_script_keeps_nonphysical_abstract_scope(capsys):
    quotient.main()
    result = json.loads(capsys.readouterr().out)
    assert result["physical_lattice_instance"] is False
    assert result["abstract_generator_masks"] == [3]
    assert result["abstract_group_masks"] == [0, 3]
    assert result["committed_rates_lambda1"]["quotient_gap"] == 4.0


def _assert_ci_executes_complete_scope(workflow):
    job = workflow["jobs"]["transfer-and-custody"]
    assert set(job["strategy"]["matrix"]["os"]) == {"ubuntu-latest"}
    assert job["env"]["OPENBLAS_NUM_THREADS"] == "1"
    assert job["env"]["OMP_NUM_THREADS"] == "1"
    assert not job.get("if") and not job.get("continue-on-error")
    assert workflow["permissions"] == {"contents": "read"}
    assert "workflow_dispatch" in workflow["on"]
    assert "code/yang_mills/**" in workflow["on"]["pull_request"]["paths"]
    assert workflow["on"]["push"]["paths"] == workflow["on"]["pull_request"]["paths"]
    dispatched = []
    for step in job["steps"]:
        if "uses" in step:
            assert len(step["uses"].rsplit("@", 1)[1]) == 40
            if step["uses"].startswith("actions/setup-python@"):
                assert step["with"]["python-version"] == "3.12"
        if "pytest" not in step.get("run", ""):
            continue
        assert not step.get("if") and not step.get("continue-on-error")
        command = shlex.split(step["run"])
        assert command[:6] == ["python", "-W", "error", "-m", "pytest", "-q"]
        # Whole-file execution excludes collect-only, selectors and shell
        # success fallbacks; counting also rejects duplicate expensive replays.
        assert all(path in REQUIRED_TESTS for path in command[6:])
        dispatched.extend(command[6:])
    assert set(dispatched) == REQUIRED_TESTS
    assert len(dispatched) == len(REQUIRED_TESTS)


def _workflow():
    return yaml.load(WORKFLOW.read_text(encoding="utf-8"), Loader=yaml.BaseLoader)


def test_ci_executes_controls_once_per_platform():
    _assert_ci_executes_complete_scope(_workflow())


@pytest.mark.parametrize("mutation", ["collect_only", "missing_test", "duplicate", "skip", "allow_failure"])
def test_ci_gate_rejects_false_green_routes(mutation):
    workflow = _workflow()
    steps = workflow["jobs"]["transfer-and-custody"]["steps"]
    step = next(row for row in steps if "pytest" in row.get("run", ""))
    if mutation == "collect_only":
        step["run"] += " --collect-only"
    elif mutation == "missing_test":
        step["run"] = step["run"].replace(sorted(REQUIRED_TESTS)[0], "")
    elif mutation == "duplicate":
        steps.append(copy.deepcopy(step))
    elif mutation == "skip":
        step["if"] = "false"
    else:
        step["continue-on-error"] = "true"
    with pytest.raises(AssertionError):
        _assert_ci_executes_complete_scope(workflow)
