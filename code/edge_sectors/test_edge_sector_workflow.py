"""Keep the full finite edge-sector controls mandatory on both CI platforms."""
from pathlib import Path
import shlex

import pytest
import yaml


WORKFLOW = Path(__file__).resolve().parents[2] / ".github/workflows/finite-gauge-transfer.yml"
INPUTS = {
    "code/edge_sectors/**", "paper/tex_fragments/PAPER.tex",
    "paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.tex",
    "claims/claim_registry.yaml", "claims/falsification_matrix.csv",
    "claims/novelty_matrix.csv", "tools/test_post_r2029_audit_surfaces.py",
    "requirements.txt", "pytest.ini", ".github/workflows/finite-gauge-transfer.yml",
}
COMMAND = ["python", "-W", "error", "-m", "pytest", "-q",
           "code/edge_sectors", "tools/test_post_r2029_audit_surfaces.py"]


def _workflow():
    return yaml.load(WORKFLOW.read_text(encoding="utf-8"), Loader=yaml.BaseLoader)


def _assert_complete_edge_sector_gate(workflow):
    assert "workflow_dispatch" in workflow["on"]
    for event in ("push", "pull_request"):
        assert INPUTS <= set(workflow["on"][event]["paths"])
    job = workflow["jobs"]["edge-sector-ground-states"]
    assert not job.get("if") and not job.get("continue-on-error")
    assert job["runs-on"] == "${{ matrix.os }}"
    assert job["strategy"]["matrix"] == {"os": ["ubuntu-latest"]}
    assert job["env"]["OPENBLAS_NUM_THREADS"] == job["env"]["OMP_NUM_THREADS"] == "1"
    assert not workflow.get("env", {}).get("PYTEST_ADDOPTS")
    assert not job.get("env", {}).get("PYTEST_ADDOPTS")
    runs = [step for step in job["steps"] if "pytest" in step.get("run", "")]
    assert len(runs) == 1
    step = runs[0]
    assert not step.get("if") and not step.get("continue-on-error")
    assert not step.get("env", {}).get("PYTEST_ADDOPTS")
    # Exact arguments reject collection-only, selectors, file-only substitutes,
    # and shell fallbacks that could turn a failed or incomplete run green.
    assert shlex.split(step["run"]) == COMMAND


def test_ci_executes_complete_edge_sector_scope_on_linux():
    _assert_complete_edge_sector_gate(_workflow())


@pytest.mark.parametrize("mutation", [
    "collect_only", "filter", "one_file", "allow_failure", "skip_job",
    "skip_step", "environment_collection", "shell_success",
])
def test_workflow_gate_rejects_false_green_execution(mutation):
    workflow = _workflow()
    job = workflow["jobs"]["edge-sector-ground-states"]
    step = next(step for step in job["steps"] if "pytest" in step.get("run", ""))
    if mutation == "collect_only":
        step["run"] += " --collect-only"
    elif mutation == "filter":
        step["run"] += " -k z2"
    elif mutation == "one_file":
        step["run"] = step["run"].replace("code/edge_sectors", "code/edge_sectors/test_heat_kernel_regressions.py")
    elif mutation == "allow_failure":
        step["continue-on-error"] = "true"
    elif mutation == "skip_job":
        job["if"] = "false"
    elif mutation == "skip_step":
        step["if"] = "false"
    elif mutation == "environment_collection":
        job["env"]["PYTEST_ADDOPTS"] = "--collect-only"
    else:
        step["run"] += " || true"
    with pytest.raises(AssertionError):
        _assert_complete_edge_sector_gate(workflow)


@pytest.mark.parametrize("path", sorted(INPUTS))
def test_workflow_gate_requires_changed_source_and_public_claim_triggers(path):
    workflow = _workflow()
    workflow["on"]["pull_request"]["paths"].remove(path)
    with pytest.raises(AssertionError):
        _assert_complete_edge_sector_gate(workflow)
