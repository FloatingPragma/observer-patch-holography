import copy
import json
from collections import Counter

import pytest

import capacity_semantics_menu_certificate as cm


@pytest.fixture(scope="module")
def payload():
    return cm.build()


def test_schema_and_row_count(payload):
    assert payload["schema"] == "oph.capacity_semantics_menu_certificate.v1"
    assert payload["issue"] == 615
    assert payload["menu_row_count"] == len(payload["menu_rows"]) == 17
    assert payload["executed_construction_rows"] == 190
    assert payload["comparison_landed_rows"] == 0


def test_campaign_verdict_fields(payload):
    verdict = payload["campaign_verdict"]
    assert verdict["semantics_enumeration"] == "declared_menu_complete_for_executed_families"
    assert verdict["source_only_fixed_point_selector"] == "conditional_open_interface"
    assert verdict["no_target_import"] == "verified"
    assert verdict["horizon_area_and_ew_bridge"] == "separate_physical_attachments (#589, #547)"


def test_family_coverage(payload):
    families = Counter(row["family"] for row in payload["menu_rows"])
    assert families == {"CAP-K": 4, "CAP-P": 6, "CAP-L": 5, "CAP-B": 1, "coupled": 1}


def test_axes_declared_before_rows(payload):
    axes = {axis["axis"] for axis in payload["semantic_axes"]}
    assert {
        "publicness_family",
        "cell_product_structure",
        "reserve_semantics",
        "reserve_attachment",
        "readback_record_effect",
        "observer_marking",
        "symmetry_quotient",
        "cap_read_family",
        "channel_capacity_semantics",
        "kernel_ontology",
        "continuation_scope",
    } <= axes
    keys = list(payload)
    assert keys.index("semantic_axes") < keys.index("menu_rows")


def test_every_verdict_in_fail_closed_vocabulary(payload):
    for row in payload["menu_rows"]:
        assert row["executed_verdict"] in cm.ALLOWED_VERDICTS


def test_capk_rows_have_no_positive_fixed_point(payload):
    for row in payload["menu_rows"]:
        if row["family"] == "CAP-K":
            assert row["executed_verdict"] == "no_positive_fixed_point"
            assert row["linear_contraction"] == {"map": "F(N) = s*N", "s_below_one": True}


def test_capp_multiplicative_rows_excluded_at_pi(payload):
    verdicts = {
        row["row_id"]: row["executed_verdict"]
        for row in payload["menu_rows"]
        if row["family"] == "CAP-P"
    }
    assert verdicts["capP.s_poisson_port"] == "excluded"
    assert verdicts["capP.s_presence_port"] == "excluded"
    assert verdicts["capP.s_poisson_pair"] == "excluded"
    assert verdicts["capP.s_presence_pair"] == "excluded"
    assert verdicts["capP.add_slot"] == "no_positive_fixed_point"
    assert verdicts["capP.add_port"] == "no_positive_fixed_point"


def test_capl_sublattices_sum_to_the_recorded_lattice(payload):
    capl = [row for row in payload["menu_rows"] if row["family"] == "CAP-L"]
    assert all(row["executed_verdict"] == "excluded" for row in capl)
    total = sum(sum(row["recorded_status_counts"].values()) for row in capl)
    assert total == 180


def test_capb_barred_pre_evaluation_and_coupled_open(payload):
    by_id = {row["row_id"]: row for row in payload["menu_rows"]}
    assert by_id["capB.bridge_constant"]["executed_verdict"] == "excluded_pre_evaluation"
    assert by_id["coupled.cp1_cp2_cp3"]["executed_verdict"] == "conditional_open"


def test_capk_positive_fixed_point_claim_fails_closed(payload):
    row = copy.deepcopy(
        next(r for r in payload["menu_rows"] if r["family"] == "CAP-K")
    )
    row["executed_verdict"] = "excluded"
    with pytest.raises(cm.CertificateError) as excinfo:
        cm.validate_menu_row(row)
    assert excinfo.value.code == "CONTRACTION_CONTRADICTION"


def test_hidden_target_scale_token_fails_closed(payload):
    row = copy.deepcopy(
        next(r for r in payload["menu_rows"] if r["family"] == "CAP-K")
    )
    row["reason"] = "row pins the readback at 9.99e122 nats by hand"
    with pytest.raises(cm.CertificateError) as excinfo:
        cm.validate_menu_row(row)
    assert excinfo.value.code == "HIDDEN_TARGET"


def test_hidden_target_marker_key_fails_closed(payload):
    row = copy.deepcopy(
        next(r for r in payload["menu_rows"] if r["family"] == "CAP-K")
    )
    row["semantics"] = dict(row["semantics"], desired_capacity_nats="pinned by hand")
    with pytest.raises(cm.CertificateError) as excinfo:
        cm.validate_menu_row(row)
    assert excinfo.value.code == "HIDDEN_TARGET"


def test_scan_accepts_digests_and_rejects_extreme_numbers():
    cm.scan_for_targets({"file_sha256": "a" * 64})
    with pytest.raises(cm.CertificateError):
        cm.scan_for_targets({"value": 1e60})
    with pytest.raises(cm.CertificateError):
        cm.scan_for_targets({"value": 1e-60})
    with pytest.raises(cm.CertificateError):
        cm.scan_for_targets({"measured_lambda": "anything"})


def test_capk_contraction_argument_is_symbolic(payload):
    control = next(
        c for c in payload["negative_controls"]
        if c["name"] == "capk_linear_contraction_symbolic"
    )
    detail = control["detail"]
    fractions = {
        row["branch"]: row.get("exact_fraction")
        for row in detail["rows"]
        if "exact_fraction" in row
    }
    assert fractions == {"capK.s_nat_share": "5/6", "capK.s_edge_share": "1/2"}
    assert all(row["recorded_s_hi_below_one"] for row in detail["rows"])


def test_negative_controls_all_passed(payload):
    names = {c["name"] for c in payload["negative_controls"]}
    assert {
        "capk_positive_fixed_point_claim_rejected",
        "hidden_target_scale_token_rejected",
        "hidden_target_marker_key_rejected",
        "capk_linear_contraction_symbolic",
    } <= names
    assert all(c["passed"] for c in payload["negative_controls"])


def test_recorded_mismatch_fails_closed(monkeypatch):
    original = cm.load_json

    def tampered(path):
        data = original(path)
        if "capK" in str(path):
            data["summary"]["by_status"] = {"fixed_point_certified": 4}
        return data

    monkeypatch.setattr(cm, "load_json", tampered)
    with pytest.raises(cm.CertificateError) as excinfo:
        cm.load_recorded()
    assert excinfo.value.code == "RECORDED_MISMATCH"


def test_manifest_on_disk_matches_deterministic_rebuild(payload):
    manifest = json.loads(cm.MANIFEST_PATH.read_text(encoding="utf-8"))
    assert manifest == payload
    cm.scan_for_targets(manifest)


def test_registry_alignment_is_read_only(payload):
    alignment = payload["axiom_registry_alignment"]
    assert alignment["rows"]["capacity_publicness_and_closure"] == "conditional_open_interface"
    assert "read-only" in alignment["access"]
