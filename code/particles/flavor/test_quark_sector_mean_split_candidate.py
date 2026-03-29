import json
import math
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SPREAD_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_spread_map.py"
MEAN_SPLIT_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_sector_mean_split.py"
DESCENT_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_sector_descent.py"
FORWARD_SCRIPT = ROOT / "particles" / "flavor" / "build_forward_yukawas.py"
MEAN_SPLIT_OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_sector_mean_split.json"
DESCENT_OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_sector_descent.json"
FORWARD_OUTPUT = ROOT / "particles" / "runs" / "flavor" / "forward_yukawas.json"


def test_quark_sector_mean_split_candidate_is_compact_and_sector_distinct() -> None:
    subprocess.run([sys.executable, str(SPREAD_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(MEAN_SPLIT_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(DESCENT_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(FORWARD_SCRIPT), "--input", str(DESCENT_OUTPUT)], check=True, cwd=ROOT)

    mean_split = json.loads(MEAN_SPLIT_OUTPUT.read_text(encoding="utf-8"))
    descent = json.loads(DESCENT_OUTPUT.read_text(encoding="utf-8"))
    forward = json.loads(FORWARD_OUTPUT.read_text(encoding="utf-8"))
    spread = json.loads((ROOT / "particles" / "runs" / "flavor" / "quark_spread_map.json").read_text(encoding="utf-8"))

    assert spread["sigma_source_kind"] == "theorem_grade_mean_surface_readback"
    assert spread["spread_emitter_status"] == "closed"
    assert mean_split["candidate_kind"] == "two_scalar_affine_mean_split"
    assert mean_split["theorem_candidate"] == "quark_sector_mean_split_theorem"
    assert mean_split["active_candidate"] == "ordered_affine_mean_readout_candidate"
    assert mean_split["proof_status"] == "closed_current_family_predictive_law"
    assert mean_split["mean_law_kind"] == "theorem_grade_two_scalar_affine_readout_closed"
    assert mean_split["smallest_constructive_missing_object"] is None
    assert mean_split["mean_denominator"] > 0.0
    assert mean_split["skew_denominator"] > 0.0
    assert mean_split["A_ud_candidate"] > 0.0
    assert mean_split["B_ud_candidate"] > 0.0
    assert mean_split["g_u_candidate"] > mean_split["g_d_candidate"] > 0.0
    assert mean_split["g_u_candidate"] < mean_split["shared_norm_value"]
    assert mean_split["g_d_candidate"] < mean_split["shared_norm_value"]
    assert descent["shared_norm_origin"] == "oph_quark_sector_mean_split"
    assert descent["g_u"] == mean_split["g_u_candidate"]
    assert descent["g_d"] == mean_split["g_d_candidate"]

    gm_u = math.prod(forward["singular_values_u"]) ** (1.0 / 3.0)
    gm_d = math.prod(forward["singular_values_d"]) ** (1.0 / 3.0)
    assert gm_u > gm_d
