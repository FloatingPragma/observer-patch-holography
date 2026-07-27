"""Regression and adversarial tests for the issue #328 settling certificate."""

from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path
import sys

import pytest


MODULE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(MODULE_DIR))

import compiled_lattice_settling_certificate as cert  # noqa: E402


MANIFEST_PATH = MODULE_DIR / "manifests" / "compiled_lattice_settling_reference.json"


@pytest.fixture(scope="module")
def manifest() -> dict:
    return json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def recomputed() -> dict:
    return cert.build_manifest()


def test_stored_manifest_is_exactly_recomputable(manifest: dict, recomputed: dict) -> None:
    assert manifest == recomputed


def test_schema_and_issue(manifest: dict) -> None:
    assert manifest["schema"] == "oph.compiled_lattice_settling_certificate.v2"
    assert manifest["issue"] == 328


def test_abstract_compiler_result_stays_separately_labelled(manifest: dict) -> None:
    block = manifest["abstract_compiler_result"]
    assert block["paper"] == "extra/observable_normal_forms.tex"
    assert any("cor:boolean-circuit-compiler" in row for row in block["paper_results"])
    assert block["lean_module"].endswith("ObservableNormalForms/Functional.lean")
    assert any(
        row.endswith("synchronous_depth_settling") for row in block["lean_declarations"]
    )
    assert "separately labelled" in block["label"]


def test_scope_names_open_hardware_attachment(manifest: dict) -> None:
    assert "open" in manifest["scope"]["open"]
    assert "hardware" in manifest["scope"]["open"]
    assert "scatter-project" in manifest["scope"]["open"]
    linkage = manifest["scope"]["scatter_project_linkage"]
    assert linkage["status"] == "blocked_by_same_interface_counterexample"
    assert linkage["uniform_settling_verdict"] is False
    assert "period-two orbit" in linkage["statement"]


def test_scatter_project_source_is_revision_locked(manifest: dict) -> None:
    source = manifest["scatter_project_correspondence_audit"]["source"]
    assert source["repository"] == "muellerberndt/ophminer"
    assert source["commit"] == "afdd1af1491d5bc293595343128c3730637c2a8d"
    assert source["path"] == "crates/fls-settle/src/lib.rs"
    assert source["file_sha256"] == (
        "154abda8798f0ab9c8a95f7f6b74005285d662669768e0fd08d5c74635aa8d21"
    )
    assert source["implementation"] == "LazySettler with UpdateLaw::Synchronous"


def test_exact_scatter_and_projection_match_source_rule() -> None:
    scattered = cert.scatter3_exact((0, 1, 1))
    assert scattered == (Fraction(4, 3), Fraction(1, 3), Fraction(1, 3))
    projected, distances = cert.project_scattered_and_exact((0, 1, 1), {1: 1})
    assert projected == (1, 1, 1)
    assert dict(distances) == {
        (0, 1, 0): Fraction(7, 3),
        (1, 1, 1): Fraction(1),
    }


def test_satisfiable_and_chain_has_nonzero_period_two_orbit(manifest: dict) -> None:
    audit = manifest["scatter_project_correspondence_audit"]
    assert audit["verdict"] == "false_for_source_identified_synchronous_rule"
    assert audit["lattice"]["compiled_depth"] == 2
    assert audit["lattice"]["zero_mismatch_solution_count"] == 1
    assert audit["lattice"]["unique_zero_mismatch_state"] == {"w": 0, "out": 0}
    assert audit["cycle"] == {
        "entry_step": 0,
        "period": 2,
        "states": [{"w": 0, "out": 1}, {"w": 1, "out": 1}],
        "phi": [1, 1],
        "reaches_zero_mismatch": False,
    }
    states = [
        (row["state"]["w"], row["state"]["out"]) for row in audit["trajectory"]
    ]
    assert states == [(0, 1), (1, 1), (0, 1)]
    assert [row["phi"] for row in audit["trajectory"]] == [1, 1, 1]


def test_same_interface_intertwiner_fails_on_counterexample(manifest: dict) -> None:
    witness = manifest["scatter_project_correspondence_audit"][
        "same_interface_intertwiner"
    ]
    assert witness["commutes_on_witness"] is False
    assert witness["scatter_project_next"] == {"w": 1, "out": 1}
    assert witness["ranked_functional_next"] == {"w": 0, "out": 0}
    assert cert.ranked_and_chain_step((0, 1)) == (0, 0)
    assert cert.scatter_project_and_chain_step((0, 1))[0] == (1, 1)


def test_primitive_verification_is_exhaustive(manifest: dict) -> None:
    prims = manifest["primitives"]
    assert set(prims) == {"NAND", "WIRE", "FANOUT2"}
    # domain sizes: 2^(ports + registers)
    assert prims["NAND"]["evidence_bundle"]["checked_port_state_pairs"] == 8
    assert prims["WIRE"]["evidence_bundle"]["checked_port_state_pairs"] == 4
    assert prims["FANOUT2"]["evidence_bundle"]["checked_port_state_pairs"] == 8
    for report in prims.values():
        bundle = report["evidence_bundle"]
        assert bundle["state_independent_update"] is True
        assert bundle["settles_in_rounds"] == 1
        assert bundle["settled_state_is_fixed_point"] is True


def test_settling_constant_and_per_circuit_bounds(manifest: dict) -> None:
    assert manifest["settling_theorem"]["constant_c"] == 1
    expected = {
        "xor_from_nand": (6, 6),
        "nand_tree_depth3": (4, 4),
        "two_bit_adder_from_nand": (13, 13),
    }
    for name, (depth, worst) in expected.items():
        block = manifest["circuits"][name]
        assert block["compiled_depth"] == depth
        settling = block["settling"]
        assert settling["worst_settling_time"] == worst
        assert settling["worst_settling_time"] <= 1 * block["compiled_depth"]
        assert settling["bound_verified"] is True
        assert settling["min_potential_margin"] >= 0


def test_exhaustive_coverage_flags(manifest: dict) -> None:
    assert manifest["circuits"]["xor_from_nand"]["settling"]["coverage"].startswith("exhaustive")
    assert manifest["circuits"]["nand_tree_depth3"]["settling"]["coverage"].startswith(
        "exhaustive"
    )
    assert manifest["circuits"]["two_bit_adder_from_nand"]["settling"]["coverage"].startswith(
        "declared finite family"
    )
    # exhaustive trajectory counts: 2^inputs * 2^registers
    assert manifest["circuits"]["xor_from_nand"]["settling"]["trajectories"] == 4 * 2**11
    assert manifest["circuits"]["nand_tree_depth3"]["settling"]["trajectories"] == 256 * 2**8


def test_depth_weight_countermodel_is_recorded(manifest: dict) -> None:
    record = manifest["potential"]["depth_weight_countermodel"]
    assert record["depth_weight_potential_after"] > record["depth_weight_potential_before"]
    drop = record["path_count_potential_before"] - record["path_count_potential_after"]
    assert drop >= len(record["witness_unsettled_before"])


def test_all_controls_fail_closed(manifest: dict) -> None:
    expected_codes = {
        "noisy_update_table": "INTERTWINER_BROKEN",
        "cross_talk_fanout": "CROSS_TALK",
        "cyclic_layout": "CYCLIC_LAYOUT",
        "fan_in_overload": "PORT_ARITY",
        "multi_consumer_port": "MULTI_CONSUMER",
        "malformed_source_reference": "MALFORMED_SOURCE_REF",
        "unknown_source_primitive": "UNKNOWN_PRIMITIVE",
    }
    controls = manifest["controls"]
    assert set(controls) == set(expected_codes)
    for name, code in expected_codes.items():
        assert controls[name]["verdict"] == "fails_closed"
        assert controls[name]["error_code"] == code


def test_noisy_table_breaks_the_intertwiner() -> None:
    with pytest.raises(cert.CertificateError) as excinfo:
        cert.verify_primitive(cert.noisy_nand_control())
    assert excinfo.value.code == "INTERTWINER_BROKEN"


def test_cross_talking_fanout_is_caught() -> None:
    with pytest.raises(cert.CertificateError) as excinfo:
        cert.verify_primitive(cert.crosstalk_fanout_control())
    assert excinfo.value.code == "CROSS_TALK"


def test_cyclic_layout_is_rejected() -> None:
    prims = cert.reference_primitives()
    cyclic = cert.Netlist(
        "cycle",
        ("x",),
        (
            cert.Instance("W1", "WIRE", (("out", "W2", "z"),)),
            cert.Instance("W2", "WIRE", (("out", "W1", "z"),)),
        ),
        (),
    )
    with pytest.raises(cert.CertificateError) as excinfo:
        cert.check_netlist(cyclic, prims)
    assert excinfo.value.code == "CYCLIC_LAYOUT"


def test_fan_in_overload_is_rejected() -> None:
    prims = cert.reference_primitives()
    overload = cert.Netlist(
        "overload",
        ("x", "y", "z"),
        (cert.Instance("G", "NAND", (("in", "x"), ("in", "y"), ("in", "z"))),),
        (),
    )
    with pytest.raises(cert.CertificateError) as excinfo:
        cert.check_netlist(overload, prims)
    assert excinfo.value.code == "PORT_ARITY"


def test_unmediated_fan_out_is_rejected() -> None:
    prims = cert.reference_primitives()
    multi = cert.Netlist(
        "multi",
        ("x",),
        (cert.Instance("G", "NAND", (("in", "x"), ("in", "x"))),),
        (),
    )
    with pytest.raises(cert.CertificateError) as excinfo:
        cert.check_netlist(multi, prims)
    assert excinfo.value.code == "MULTI_CONSUMER"


@pytest.mark.parametrize(
    "bad_ref",
    [
        (),
        ("out", "W"),
        ("unknown-tag", "x"),
        ("in", 7),
        ("out", "W", 0),
        ["in", "x"],
    ],
)
def test_malformed_source_references_have_stable_error(bad_ref: object) -> None:
    prims = cert.reference_primitives()
    malformed = cert.Netlist(
        "malformed",
        ("x",),
        (cert.Instance("W", "WIRE", (bad_ref,)),),
        (),
    )
    with pytest.raises(cert.CertificateError) as excinfo:
        cert.check_netlist(malformed, prims)
    assert excinfo.value.code == "MALFORMED_SOURCE_REF"


def test_forward_reference_to_unknown_source_kind_has_stable_error() -> None:
    prims = cert.reference_primitives()
    unknown_source_kind = cert.Netlist(
        "unknown_source_kind",
        ("x",),
        (
            cert.Instance("W", "WIRE", (("out", "BAD", "z"),)),
            cert.Instance("BAD", "UNKNOWN_KIND", (("in", "x"),)),
        ),
        (),
    )
    with pytest.raises(cert.CertificateError) as excinfo:
        cert.check_netlist(unknown_source_kind, prims)
    assert excinfo.value.code == "UNKNOWN_PRIMITIVE"


def test_tampered_gate_table_is_detected_end_to_end() -> None:
    prims = cert.reference_primitives()
    broken = cert.build_primitive(
        "NAND",
        in_ports=("x", "y"),
        registers=("q",),
        out_ports=("z",),
        readback=(0,),
        truth_fn=lambda i: (1 - (i[0] & i[1]),),
        next_state_fn=lambda i, s: (1,) if i == (1, 1) else (1 - (i[0] & i[1]),),
    )
    with pytest.raises(cert.CertificateError) as excinfo:
        cert.verify_primitive(broken)
    assert excinfo.value.code == "INTERTWINER_BROKEN"
    assert prims["NAND"].update != broken.update


def test_xor_extension_matches_gate_semantics_and_bound() -> None:
    prims = cert.reference_primitives()
    circuit = cert.xor_circuit()
    net = cert.compile_gates_to_patches(circuit)
    cn = cert.compile_net(net, prims)
    assert cert.verify_extension_semantics(circuit, cn) == 4
    report = cert.exhaustive_settling_report(circuit, cn)
    assert report["worst_settling_time"] == cn.analysis["depth"] == 6


def test_adder_extension_matches_two_bit_addition() -> None:
    prims = cert.reference_primitives()
    circuit = cert.adder_circuit()
    net = cert.compile_gates_to_patches(circuit)
    cn = cert.compile_net(net, prims)
    for a in range(4):
        for b in range(4):
            inputs_vec = (a & 1, (a >> 1) & 1, b & 1, (b >> 1) & 1)
            got = cert.output_values(cn, cert.extension_state(cn, inputs_vec))
            total = a + b
            assert got == {"s0": total & 1, "s1": (total >> 1) & 1, "c2": (total >> 2) & 1}


def test_manifest_emission_is_deterministic(tmp_path: Path, recomputed: dict) -> None:
    out_a = tmp_path / "a.json"
    out_b = tmp_path / "b.json"
    cert.write_json(out_a, recomputed)
    cert.write_json(out_b, cert.build_manifest())
    assert out_a.read_bytes() == out_b.read_bytes()
    assert MANIFEST_PATH.read_bytes() == out_a.read_bytes()


def test_manifest_paths_are_platform_independent(manifest: dict) -> None:
    assert "\\" not in manifest["scope"]["realized"]
    assert "\\" not in manifest["abstract_compiler_result"]["lean_module"]
