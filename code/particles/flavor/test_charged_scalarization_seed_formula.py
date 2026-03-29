import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "flavor" / "derive_charged_budget_pushforward.py"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "charged_budget_transport.json"


def test_charged_budget_exports_seed_formula_candidate() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    cert = payload["charged_dirac_scalarization_certificate"]
    assert cert["artifact"] == "charged_dirac_common_refinement_gluing_certificate"
    assert cert["seed_formula"] == "sigma_seed = mean(family_eigenvalues) + min(spectral_gaps)"
    assert cert["law_scope"] == "direct_sum_u_plus_d_plus_e_pre_normal_form_canonical_decomposable_subfamily"
    assert cert["decomposition_independence_status"] == "candidate_only"
    assert cert["presentation_equalizer_kind"] == "sector_projected_sigma_seed_equalizer"
    assert cert["minimal_missing_witness"] == "common_refinement_transport_equalizer"
    assert cert["exact_blocking_clause"] == "ordered_common_refinement_seed_rigidity"
    assert cert["seed_rigidity_reduction"] == ["mean_eigenvalue_invariance", "min_gap_invariance"]
    assert cert["smaller_exact_missing_clause"] == "common_refinement_preserves_mean_eigenvalue_and_min_gap"
    assert cert["mean_gap_invariance_candidate"]["candidate_id"] == "common_refinement_mean_and_min_gap_invariance"
    assert cert["strictly_smaller_next_subclause"] == "common_refinement_preserves_mean_eigenvalue"
    assert cert["mean_side_support_status"] == "tightest_remaining_clause"
    assert cert["gap_side_support_status"] == "proxy_supported_on_current_family"
    assert cert["current_proxy_support"]["same_label_only"] is True
