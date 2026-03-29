import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "flavor" / "derive_charged_budget_pushforward.py"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "charged_dirac_scalarization_gluing.json"


def test_shared_charged_budget_stays_open_without_gluing_theorem() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["proof_status"] == "candidate_only"
    assert payload["functional_equalizer_closed"] is False
    assert payload["minimal_missing_bridge_object"] == "charged_dirac_scalarization_gluing_certificate"
