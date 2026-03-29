import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
ODD_FORM_SCRIPT = ROOT / "particles" / "flavor" / "derive_charged_dirac_odd_deformation_form.py"
RESPONSE_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_odd_response_law.py"
DESCENT_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_sector_descent.py"
FORWARD_SCRIPT = ROOT / "particles" / "flavor" / "build_forward_yukawas.py"
RESPONSE_OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_odd_response_law.json"
FORWARD_OUTPUT = ROOT / "particles" / "runs" / "flavor" / "forward_yukawas.json"
DESCENT_INPUT = ROOT / "particles" / "runs" / "flavor" / "quark_sector_descent.json"


def test_quark_lane_uses_single_kappa_candidate_and_drops_local_shape_blocker() -> None:
    subprocess.run([sys.executable, str(ODD_FORM_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(RESPONSE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(DESCENT_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(FORWARD_SCRIPT), "--input", str(DESCENT_INPUT)], check=True, cwd=ROOT)
    response = json.loads(RESPONSE_OUTPUT.read_text(encoding="utf-8"))
    descent = json.loads(DESCENT_INPUT.read_text(encoding="utf-8"))
    forward = json.loads(FORWARD_OUTPUT.read_text(encoding="utf-8"))
    assert response["coefficient_free"] is True
    assert response["lift_constants"] is None
    assert response["lift_parameterization_kind"] == "single_kappa_signed_projector_ray"
    assert response["delta_logg_q_status"] == "closed_zero_corollary"
    assert descent["shared_norm_origin"] in {
        "charged_dirac_scalarization_restriction",
        "quark_shared_absolute_norm_binding",
        "oph_quark_sector_mean_split",
    }
    assert descent["u_raw_channel_norm_candidate"] == descent["g_u"]
    assert descent["d_raw_channel_norm_candidate"] == descent["g_d"]
    assert "quark_descent_candidate_only" not in forward["promotion_blockers"]
