from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

import pytest

import a5_multipole_fixed_point_certificate as legacy
import a5_multipole_fixed_point_hardening_certificate as hardening


def test_v2_status_and_persistence_boundary_are_explicit() -> None:
    receipt = hardening.build_receipt()
    assert receipt["schema"] == "oph.a5_multipole_fixed_point_receipt.v2"
    assert receipt["status"] == (
        "EXACT_A5_FINGERPRINT_CORE__QUANTITATIVE_PERSISTENCE_OPEN__"
        "PHYSICAL_MAP_OPEN"
    )
    persistence = receipt["critical_points"]["quantitative_persistence"]
    assert persistence["status"] == "OPEN"
    assert len(persistence["missing"]) == 4
    assert "no explicit radius" in persistence["claim_boundary"]


def test_append_only_parent_preserves_fz11_input() -> None:
    receipt = hardening.build_receipt()
    pin = receipt["parent_pins"][0]
    legacy_bytes = hardening.LEGACY_RECEIPT_PATH.read_bytes()
    assert pin["sha256"] == legacy.tagged_sha256(legacy_bytes)
    assert pin["receipt_sha256"] == json.loads(legacy_bytes)["receipt_sha256"]
    assert "FZ-11 pins" in receipt["supersedes"]["reason"]


def test_geometry_is_ordered_oriented_and_equal_trace() -> None:
    geometry = hardening.build_receipt()["cartesian_frame"]["serialized_geometry"]
    ports = geometry["ordered_ports"]
    assert len(ports) == 12
    assert [row["port"] for row in ports] == list(hardening.PORT_IDS)
    assert all(row["norm_squared"] == "5/2+1/2*sqrt5" for row in ports)
    assert all(ports[row["antipode_index"]]["antipode_index"] == row["index"] for row in ports)
    assert len(geometry["oriented_faces"]) == 20
    determinants = []
    for face in geometry["oriented_faces"]:
        rational, radical = face["determinant"].split("+", 1)
        determinants.append(
            (Fraction(rational), Fraction(radical.removesuffix("*sqrt5")))
        )
    assert all(legacy.q5_sign(value) > 0 for value in determinants)
    equal_trace = geometry["equal_trace_declaration"]
    assert equal_trace["declared"] is True
    assert equal_trace["weights"] == ["1"] * 12
    assert geometry["geometry_sha256"].startswith("sha256:")


def test_original_unsquared_equations_accept_exactly_sixty_two_rays() -> None:
    replay = hardening.build_receipt()["critical_points"]["original_equation_replay"]
    assert replay["accepted_directions"] == {
        "vertex": 12,
        "face": 20,
        "edge": 30,
        "total": 62,
    }
    assert replay["antipodal_pairs"] == 31
    assert sorted(replay["accepted_latitude_counts"].values()) == [2, 10, 10, 10, 10, 10, 10]
    signs = replay["squared_meridian_sign_enumeration"]
    assert signs["squared_candidate_directions"] == 100
    assert signs["accepted_by_original_equation"] == 50
    assert signs["rejected_as_squaring_extraneous"] == 50
    assert len(signs["latitudes"]) == 5
    for latitude in signs["latitudes"]:
        assert latitude["squared_candidate_directions"] == 20
        assert latitude["accepted_by_original_equation"] == 10
        assert latitude["rejected_as_squaring_extraneous"] == 10
        assert sum(
            branch["azimuth_count"]
            for branch in latitude["branches"]
            if branch["passes_original_unsquared_equation"]
        ) == 10


def test_stationary_replay_rejects_polynomial_mutation() -> None:
    verts = legacy.cartesian_vertices()
    polynomial = hardening.homogeneous_i6(verts)
    mutated = dict(polynomial)
    mutated[(5, 1, 0)] = legacy.q5(1)
    with pytest.raises(legacy.FingerprintError):
        hardening.exact_stationary_ray_replay(verts, mutated)


def test_response_parent_is_verified_and_correctly_typed() -> None:
    receipt = hardening.build_receipt()
    pin = receipt["parent_pins"][1]
    assert pin["artifact_sha256"].startswith("sha256:")
    provenance = receipt["band_response"]["response_provenance"]
    assert provenance["antipode_J"].startswith("derived exactly")
    assert "issue-599" in provenance["representative_R_minus_J"]


def test_semantically_rehashed_response_tamper_fails(tmp_path: Path) -> None:
    artifact = json.loads(hardening.SEMANTIC_RESPONSE_PATH.read_bytes())
    artifact["source_response"]["operator"] = "positive_graph_antipode_involution"
    artifact["artifact_sha256"] = hardening.artifact_self_hash(artifact)
    path = tmp_path / "semantic_response.json"
    path.write_text(json.dumps(artifact), encoding="utf-8")
    with pytest.raises(legacy.FingerprintError):
        hardening.build_receipt(semantic_response_path=path)


def test_controls_use_on_sphere_geometry_and_matrix_response_mutation() -> None:
    controls = {
        row["control"]: row
        for row in hardening.build_receipt()["fail_closed_controls"]["controls"]
    }
    geometry = controls["one antipodal pair rotated on the exact common sphere"]
    assert geometry["common_radius_preserved"] is True
    assert geometry["antipode_preserved"] is True
    response = controls["response multiplier mutation from -J to +J"]
    assert response["production_signs_levels_0_to_3"] == [-1, 1, -1, 1]
    assert response["mutated_signs_levels_0_to_3"] == [1, -1, 1, -1]
    assert "matrix multiplication" in response["implementation"]


def test_m1_m3_are_local_rows_not_registered_predictions() -> None:
    rules = hardening.build_receipt()["decision_rules_and_local_rows"]
    assert "local_certificate_rows" in rules
    assert "frozen_rows" not in rules
    for row in rules["local_certificate_rows"].values():
        assert row["type"] == "local_nonphysical_conditional_theorem_row"
        assert "not a separately custodied" in row["registration_scope"]
    assert "issue-599" in rules["local_certificate_rows"]["OPH-A5-M3"]["premise_ancestry"]


def test_committed_v2_receipt_is_byte_exact() -> None:
    assert hardening.RECEIPT_PATH.read_bytes() == legacy.canonical_json_bytes(
        hardening.build_receipt()
    )
