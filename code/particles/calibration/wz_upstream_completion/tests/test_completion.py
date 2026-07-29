from pathlib import Path
import copy
import hashlib
import importlib.util
import json

from jsonschema import Draft202012Validator

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "checker", ROOT / "checkers" / "check_completion_bundle.py"
)
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)

manifest_spec = importlib.util.spec_from_file_location(
    "integration_manifest", ROOT / "generate_integration_manifest.py"
)
integration_manifest = importlib.util.module_from_spec(manifest_spec)
manifest_spec.loader.exec_module(integration_manifest)


def base():
    return copy.deepcopy(checker.validate_all())


def schema_errors(key, obj):
    schema = json.loads((ROOT / "schemas" / checker.SCHEMAS[key]).read_text())
    return list(Draft202012Validator(schema).iter_errors(obj))


def test_templates_schema_valid_and_fail_closed():
    d = base()
    assert len(d) == 10
    assert checker.promotion_reasons(d)
    status = checker.aggregate_status(d)
    assert status["distinct_schema_documents"] == 9
    assert status["receipt_instances"] == 10
    assert status["promotion_allowed"] is False


def _forge_all_declarations_without_artifacts(x):
    """Create the known self-attestation attack fixture.

    This deliberately creates no referenced artifact and resolves no digest.
    It exercises the boundary between predicate wiring and evidence checking.
    """
    if isinstance(x, dict):
        for key, value in x.items():
            if isinstance(value, bool):
                x[key] = True
            elif value == checker.ZERO_HASH:
                x[key] = "1" * 64
            else:
                _forge_all_declarations_without_artifacts(value)
    elif isinstance(x, list):
        for value in x:
            _forge_all_declarations_without_artifacts(value)


def test_forged_self_attestation_can_satisfy_candidate_but_never_promotes():
    d = base()
    _forge_all_declarations_without_artifacts(d)
    d["action"]["claim_lane"] = "OPH_NATIVE_PHYSICAL"
    d["action"]["geometry"]["origin"] = "oph_receipt"
    d["matching"]["claim_lane"] = "OPH_NATIVE_PHYSICAL"
    d["fj"]["engines"]["shared_generated_expressions"] = False
    d["fj"]["engines"]["shared_integral_backend"] = False
    d["renorm"]["counterterm_generation"]["handwritten_vertex_list"] = False
    d["brst"]["independence"]["shared_generated_expressions"] = False
    d["brst"]["independence"]["shared_integral_backend"] = False
    d["pole_w"]["current_amplitude"]["positive_residue_required"] = False
    d["pole_z"]["current_amplitude"]["positive_residue_required"] = False
    d["clock"]["dimensionless_gap_over_Estar"] = {
        "lo": 1.0,
        "hi": 1.0,
        "positive_gap_certified": True,
    }

    # This demonstrates why promotion_reasons is only a declaration linter.
    assert checker.promotion_reasons(d) == []
    status = checker.aggregate_status(d)
    assert status["candidate_conjunction_satisfied"] is True
    assert status["production_artifacts_resolved"] is False
    assert status["production_digests_recomputed"] is False
    assert status["promotion_allowed"] is False


def test_target_ancestry_mutation_fails():
    d = base()
    d["action"]["source_ancestry"]["target_blacklist_passed"] = False
    assert any("target ancestry" in x for x in checker.promotion_reasons(d))


def test_imported_geometry_cannot_be_oph_native_physical():
    d = base()
    d["action"]["claim_lane"] = "OPH_NATIVE_PHYSICAL"
    d["action"]["geometry"]["origin"] = "imported_minkowski"
    assert any("geometry is imported" in x for x in checker.promotion_reasons(d))


def test_chart_lane_fails():
    d = base()
    d["action"]["claim_lane"] = "OPH_CHART_ONLY"
    assert any("action lane" in x for x in checker.promotion_reasons(d))


def test_matching_root_uniqueness_required():
    d = base()
    d["matching"]["source_root_unique"] = False
    assert any("source root not unique" in x for x in checker.promotion_reasons(d))


def test_matching_interval_ordering_required():
    d = base()
    d["matching"]["intervals"][0]["q_low_over_Estar"] = 2.0
    d["matching"]["intervals"][0]["q_high_over_Estar"] = 1.0
    assert any("interval ordering" in x for x in checker.promotion_reasons(d))


def test_pure_sm_interval_rejects_mssm_coefficients():
    d = base()
    d["matching"]["intervals"][0]["beta"]["gauge_coefficients"] = {
        "gprime": "33/5", "g": "1", "gs": "-3"
    }
    assert any("non-SM one-loop" in x for x in checker.promotion_reasons(d))


def test_beta_must_be_derived_and_independently_checked():
    d = base()
    d["matching"]["intervals"][0]["beta"]["derived_from_census"] = False
    assert any("census-derived" in x for x in checker.promotion_reasons(d))


def test_top_only_or_incomplete_yukawa_fails():
    d = base()
    d["yukawa"]["evidence"]["complete_open_channels"] = False
    assert any("Yukawa" in x for x in checker.promotion_reasons(d))


def test_massless_approx_without_remainder_fails():
    d = base()
    d["yukawa"]["approximation"] = {
        "mode": "massless_light_with_certified_remainder",
        "remainder_bound_supplied": False,
    }
    assert any("approximation lacks remainder" in x for x in checker.promotion_reasons(d))


def test_fj_complete_map_required():
    d = base()
    d["fj"]["map"]["complete_parameter_map"] = False
    assert any("FJ equivalence" in x for x in checker.promotion_reasons(d))


def test_fj_independent_engines_required():
    d = base()
    d["fj"]["engines"]["shared_integral_backend"] = True
    assert any("FJ engines not independent" in x for x in checker.promotion_reasons(d))


def test_handwritten_counterterms_fail():
    d = base()
    d["renorm"]["counterterm_generation"]["handwritten_vertex_list"] = True
    assert any("handwritten" in x for x in checker.promotion_reasons(d))


def test_chiral_finite_restoration_must_be_declared():
    d = base()
    d["renorm"]["gamma5"]["finite_restoration_declared"] = False
    assert any("renormalization/ST" in x for x in checker.promotion_reasons(d))


def test_uv_or_st_failure_fails():
    d = base()
    d["renorm"]["uv_poles"]["exact_cancellation"] = False
    assert any("renormalization/ST" in x for x in checker.promotion_reasons(d))


def test_counterterms_must_come_from_bare_action():
    d = base()
    d["renorm"]["counterterm_generation"]["from_bare_substitution"] = False
    assert any("complete bare action" in x for x in checker.promotion_reasons(d))


def test_incomplete_diagram_universe_fails():
    d = base()
    d["brst"]["diagram_universe"]["complete"] = False
    assert any("diagram universe" in x for x in checker.promotion_reasons(d))


def test_full_sm_brst_receipt_includes_qcd_gauge_parameter():
    d = base()
    assert "xiS" in d["brst"]["gauge_parameters"]
    obj = copy.deepcopy(d["brst"])
    obj["gauge_parameters"].remove("xiS")
    assert schema_errors("brst", obj)


def test_missing_neutral_mixing_block_fails_schema():
    obj = base()["brst"]
    obj["blocks"].remove("AZ_T")
    assert schema_errors("brst", obj)


def test_order_2_5_cannot_use_strict_one_loop_receipt():
    obj = base()["brst"]
    obj["evidence"]["order"] = "order_2_5"
    assert schema_errors("brst", obj)


def test_nielsen_failure_fails():
    d = base()
    d["brst"]["identities"]["nielsen_matrix"] = False
    assert any("Nielsen" in x or "identities" in x for x in checker.promotion_reasons(d))


def test_nonindependent_pole_engines_fail():
    d = base()
    d["brst"]["independence"]["shared_integral_backend"] = True
    assert any("pole engines not independent" in x for x in checker.promotion_reasons(d))


def test_precision_nesting_required():
    d = base()
    d["brst"]["precision"]["nested_balls"] = False
    assert any("precision" in x for x in checker.promotion_reasons(d))


def test_positive_residue_requirement_for_unstable_w_fails():
    d = base()
    d["pole_w"]["current_amplitude"]["positive_residue_required"] = True
    assert any("W unstable pole" in x for x in checker.promotion_reasons(d))


def test_missing_z_current_amplitude_fails():
    d = base()
    d["pole_z"]["current_amplitude"]["same_pole"] = False
    assert any("Z physical-current" in x for x in checker.promotion_reasons(d))


def test_wrong_boson_tag_fails():
    d = base()
    d["pole_z"]["boson"] = "W"
    assert any("wrong boson tag" in x for x in checker.promotion_reasons(d))


def test_rank_n_minus_one_laurent_hypothesis_required():
    d = base()
    d["pole_z"]["simple_zero"]["rank_at_pole"] = 0
    assert any(
        "rank-n-minus-one" in x for x in checker.promotion_reasons(d)
    )


def test_deterministic_zero_covariance_requires_uniqueness():
    d = base()
    d["law"]["mode"] = "deterministic_delta"
    d["law"]["deterministic"]["covariance_zero_exact"] = True
    d["law"]["deterministic"]["global_root_unique"] = False
    assert any("delta law lacks" in x for x in checker.promotion_reasons(d))


def test_stochastic_weights_must_be_source_justified():
    d = base()
    d["law"]["mode"] = "weighted_ensemble"
    d["law"]["active_mode_consistent"] = True
    d["law"]["stochastic"]["weights_source_justified"] = False
    assert any("stochastic source law" in x for x in checker.promotion_reasons(d))


def test_clock_positive_gap_required():
    d = base()
    d["clock"]["dimensionless_gap_over_Estar"] = {
        "lo": 0.0, "hi": 0.0, "positive_gap_certified": False
    }
    assert any("clock gap" in x for x in checker.promotion_reasons(d))


def test_clock_naturality_required():
    d = base()
    d["clock"]["naturality"]["schedule_invariant"] = False
    assert any("clock naturality" in x for x in checker.promotion_reasons(d))


def test_action_hash_mismatch_fails():
    d = base()
    d["action"]["hashes"]["action_ast"] = "1" * 64
    d["renorm"]["hashes"]["action_ast"] = "2" * 64
    assert any("action hash mismatch" in x for x in checker.promotion_reasons(d))


def test_term_mask_hash_mismatch_fails():
    d = base()
    d["fj"]["hashes"]["term_mask"] = "1" * 64
    d["brst"]["hashes"]["term_mask"] = "2" * 64
    assert any("term-mask hash mismatch" in x for x in checker.promotion_reasons(d))


def test_placeholder_hashes_fail():
    d = base()
    assert any("placeholder hash" in x for x in checker.promotion_reasons(d))


def test_frozen_gauge_grid_has_45_unique_points_and_hash():
    obj = json.loads((ROOT / "data" / "nonlinear_gauge_grid_v1.json").read_text())
    assert obj["scope"] == "strict_one_loop_electroweak_wz_stress_grid"
    assert obj["fixed_parameters"] == {"xiS": 1.0}
    assert obj["mixed_qcd_orders_require_additional_xiS_receipt"] is True
    assert obj["count"] == 45 == len(obj["points"])
    assert len({p["id"] for p in obj["points"]}) == 45
    given = obj.pop("canonical_sha256")
    canonical = json.dumps(obj, separators=(",", ":"), sort_keys=True).encode()
    assert hashlib.sha256(canonical).hexdigest() == given


def test_dependency_dag_is_acyclic_and_promotion_false():
    obj = json.loads((ROOT / "data" / "receipt_dependency_dag_v4.json").read_text())
    edges = obj["edges"]
    nodes = set(obj["nodes"])
    indeg = {n: 0 for n in nodes}
    adj = {n: [] for n in nodes}
    for a, b in edges:
        adj[a].append(b)
        indeg[b] += 1
    q = [n for n, degree in indeg.items() if degree == 0]
    seen = []
    while q:
        n = q.pop()
        seen.append(n)
        for m in adj[n]:
            indeg[m] -= 1
            if indeg[m] == 0:
                q.append(m)
    assert len(seen) == len(nodes)
    given = obj.pop("canonical_sha256")
    canonical = json.dumps(obj, separators=(",", ":"), sort_keys=True).encode()
    assert hashlib.sha256(canonical).hexdigest() == given
    assert obj["nodes"]["OPH_NATIVE_PHYSICAL"]["status"] == "false"


def test_post_audit_integration_manifest_is_complete_and_current():
    manifest = json.loads((ROOT / "INTEGRATION_MANIFEST.json").read_text())
    expected_paths = {
        path.relative_to(ROOT).as_posix()
        for path in integration_manifest.included_files()
    }
    records = {record["path"]: record for record in manifest["files"]}
    assert set(records) == expected_paths
    assert manifest["baseline_specification_schema_documents"] == 9
    assert manifest["diagnostic_schema_documents"] == 1
    assert manifest["total_schema_documents"] == 10
    assert manifest["total_schema_documents"] == len(
        list((ROOT / "schemas").glob("*.schema.json"))
    )
    assert manifest["baseline_specification_receipt_instances"] == 10
    assert manifest["wz_boundary_diagnostic_receipts"] == 6
    assert manifest["baseline_exact_symbolic_checks"] == 8
    assert manifest["wz_exact_finite_corrections"] == 3
    assert manifest["baseline_specification_tests"] == 38
    assert "pytest collection" in manifest["package_test_count_policy"]
    assert manifest["promotion_allowed"] is False
    for rel, record in records.items():
        path = ROOT / rel
        assert path.stat().st_size == record["bytes"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == record["sha256"]
