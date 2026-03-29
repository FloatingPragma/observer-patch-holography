import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "neutrino" / "derive_majorana_overlap_defect_scalar_evaluator.py"
OUTPUT = ROOT / "particles" / "runs" / "neutrino" / "majorana_overlap_defect_scalar_evaluator.json"


def test_scalar_evaluator_remains_open_until_cubic_freedom_eliminated() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["proof_status"] == "candidate_only"
    assert payload["invariant_ring_obstruction"]["cubic_freedom_eliminated"] is False
    assert payload["remaining_theorem_object"] == "oph_majorana_scalar_from_centered_edge_norm"
    assert payload["cubic_freedom_status_if_closed"] == "eliminated_exactly"
