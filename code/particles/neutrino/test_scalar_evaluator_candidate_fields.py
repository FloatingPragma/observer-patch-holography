import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "neutrino" / "derive_majorana_overlap_defect_scalar_evaluator.py"
OUTPUT = ROOT / "particles" / "runs" / "neutrino" / "majorana_overlap_defect_scalar_evaluator.json"


def test_scalar_evaluator_exports_constructive_centered_candidate() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["formula_affine_candidate"] == "sum_e mu_e*(1-cos(psi_e-psi_e_star))"
    assert payload["intrinsic_witness_kind"] == "centered_edge_character_norm_defect"
    assert payload["origin_object_candidate"] == "oph_affine_majorana_edge_character_functor"
    assert payload["normalization_witness_status"] == "local_quadratic_match_only"
    assert payload["naive_uncentered_formula_ruled_out"] is True
    assert payload["forced_bundle_norm_formula"] == "Q_nu(xi)=0.5*sum_e mu_e*abs(xi_e)^2"
    assert payload["forced_phase_formula"] == "S_nu(psi)=sum_e mu_e*(1-cos(psi_e-psi_e_star))"
