"""Adversarial input and interpretation checks for conditional RSCC replay."""
import copy
import importlib.util
import json
from pathlib import Path

import pytest

PATH = Path(__file__).with_name("conditional_quark_mass_replay.py")
IMPORT_SPEC = importlib.util.spec_from_file_location("conditional_quark_replay_test_subject", PATH)
replay = importlib.util.module_from_spec(IMPORT_SPEC)
IMPORT_SPEC.loader.exec_module(replay)


@pytest.fixture(scope="module")
def spec():
    return replay.load_spec()


@pytest.fixture(scope="module")
def payload():
    return replay.build_payload()


def test_six_conditional_coordinates_and_fixed_control(payload):
    for name in ("full_rscc", "lower_order_control"):
        assert set(payload[name]["native_mass_chart_GeV"]) == {"u", "d", "s", "c", "b", "t"}
        assert all(float(x) > 0 for x in payload[name]["native_mass_chart_GeV"].values())
    assert payload["model"]["primary"] == "full_rscc"
    assert payload["model"]["control"] == "lower_order_ablation_not_selected"
    assert payload["guards"]["independent_prediction"] is False
    assert payload["guards"]["controls_used_for_model_selection"] is False
    assert payload["guards"]["formula_discovery_target_informed"] is True
    assert payload["guards"]["common_scale_yukawa_matrix"] is False


@pytest.mark.parametrize("path,value", [
    (("calibrations", "alpha_inverse", "value"), "137"),
    (("calibrations", "alpha_inverse", "role"), "source_prediction"),
    (("calibrations", "E_star_GeV", "value"), "2.435e18"),
    (("calibrations", "E_star_GeV", "role"), "derived_absolute_unit"),
    (("coefficient_menu", "mean_d_quadratic"), "1/432"),
    (("coefficient_menu", "a_u_quadratic"), "1/23"),
    (("structural_counts", "oriented_repair_slots"), 12),
    (("charts", "t", "scheme"), "MSbar"),
    (("charts", "c", "scale"), "2 GeV"),
    (("guards", "physical_quark_mass_derivation"), True),
    (("guards", "formula_discovery_target_informed"), False),
    (("hypotheses", "H11_representation_attachment"), "Already proved"),
    (("gauge_solver", "su3_cutoff"), 30),
])
def test_reject_mutated_law_unit_chart_or_claim(spec, tmp_path, path, value):
    changed = copy.deepcopy(spec)
    cursor = changed
    for component in path[:-1]:
        cursor = cursor[component]
    cursor[path[-1]] = value
    destination = tmp_path / "spec.json"
    destination.write_text(json.dumps(changed))
    with pytest.raises(ValueError):
        replay.load_spec(destination)


def test_reject_hidden_target_input(spec, tmp_path):
    changed = copy.deepcopy(spec)
    changed["quark_targets"] = {"top": "172.1"}
    destination = tmp_path / "spec.json"
    destination.write_text(json.dumps(changed))
    with pytest.raises(ValueError, match="hidden input"):
        replay.load_spec(destination)


@pytest.mark.parametrize("field", ["P", "alpha_U", "v_over_E_star", "mz_over_E_star", "alpha3_mz"])
def test_reject_mixed_source_branch(spec, field):
    packet = copy.deepcopy(spec["source_cache"])
    packet[field] = str(float(packet[field]) * 1.0001)
    with pytest.raises(ValueError):
        replay.verify_source_packet(packet, spec)


def test_reject_source_packet_mass_injection(spec):
    packet = copy.deepcopy(spec["source_cache"])
    packet["masses"] = {"u": ".00216"}
    with pytest.raises(ValueError, match="target injection"):
        replay.verify_source_packet(packet, spec)


def test_replay_has_no_numeric_read_from_historical_or_comparison_data(monkeypatch):
    original = Path.read_text
    seen = set()
    allowed = {replay.SPEC_PATH.resolve(), replay.ALPHA_PATH.resolve(), replay.RESPONSE_PATH.resolve(),
               (replay.ROOT/"code/particles/flavor/source_w5_response_constraints_spec.json").resolve(),
               (replay.ROOT/"code/particles/flavor/source_w5_response_constraints.py").resolve(),
               (replay.ROOT/"code/a5_closure/manifests/a3_scheduler_kernel_reference.json").resolve(),
               (replay.ROOT/"code/a5_closure/manifests/record_counting_mechanism_reference.json").resolve()}
    def guarded(path, *args, **kwargs):
        resolved = path.resolve()
        assert resolved in allowed, f"undeclared numeric read: {resolved}"
        seen.add(resolved)
        return original(path, *args, **kwargs)
    monkeypatch.setattr(Path, "read_text", guarded)
    result = replay.build_payload()
    assert seen == allowed
    assert result["guards"]["quark_reference_values_used_at_runtime"] is False


def test_receipt_matches_computation(payload):
    assert json.loads(replay.OUT_PATH.read_text()) == payload


@pytest.mark.parametrize("key,value", [("family_dimension",6),("line_complement_dimension",3),("scalar_fraction","1/6"),("centered_fraction","3/5"),("common_exposure","4/5")])
def test_certified_source_projection_is_fixed(spec,tmp_path,key,value):
    changed=copy.deepcopy(spec);changed["source_response"]["projection"][key]=value
    path=tmp_path/"spec.json";path.write_text(json.dumps(changed))
    with pytest.raises(ValueError):replay.load_spec(path)


def test_tampered_source_certificate_is_not_used(tmp_path):
    cert=json.loads(replay.RESPONSE_PATH.read_text());cert["color"]["normalized_trace"]="4/5"
    path=tmp_path/"source.json";path.write_text(json.dumps(cert))
    with pytest.raises(ValueError):replay.load_source_response(path)
