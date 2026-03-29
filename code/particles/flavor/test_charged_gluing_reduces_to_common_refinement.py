import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "flavor" / "derive_charged_budget_pushforward.py"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "charged_dirac_scalarization_gluing.json"


def test_charged_gluing_reduces_to_common_refinement() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    status = payload["generatorwise_status"]
    assert payload["generatorwise_reduction_status"] == "common_refinement_only_missing"
    assert payload["artifact"] == "charged_dirac_common_refinement_gluing_certificate"
    assert payload["exact_blocking_clause"] == "ordered_common_refinement_seed_rigidity"
    assert set(status["conjugation"].values()) == {"closed_on_current_metadata"}
    assert set(status["pre_normal_form_equivalence"].values()) == {"closed_on_current_metadata"}
    assert set(status["common_refinement"].values()) == {"candidate_only"}
    assert payload["minimal_missing_witness"] == "common_refinement_transport_equalizer"
