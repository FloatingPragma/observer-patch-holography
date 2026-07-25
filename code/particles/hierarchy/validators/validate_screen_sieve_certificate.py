#!/usr/bin/env python3
"""Validate the OPH icosahedral screen-sieve certificate."""

from __future__ import annotations

import json
import importlib.util
import pathlib
import sys
from typing import Any


def load_builder() -> Any:
    builder_path = pathlib.Path(__file__).resolve().parents[1] / "verify_screen_sieve_theorem.py"
    spec = importlib.util.spec_from_file_location(
        "screen_sieve_exact_builder",
        builder_path,
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {builder_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def validate(certificate: dict[str, Any]) -> dict[str, bool]:
    # Rebuild every enumerable/arithmetic field from the antecedent code.
    # The stored ``checks`` object is diagnostic output, never an authority.
    expected = load_builder().build_certificate()
    poly = {item["name"]: item for item in certificate.get("polyhedral_comparison", [])}
    seidel = certificate.get("d_optimal_tomography", {}).get("seidel_uniqueness", {})
    graph_class = seidel.get("normalized_negative_graph_class", {})
    arithmetic = certificate.get("screen_load_arithmetic", {})
    readout_gate = certificate.get("hierarchy_screen_readout_gate", {})
    checks = {
            "certificate_exactly_matches_recomputation": certificate == expected,
            "embedded_check_booleans_match_recomputation": (
                certificate.get("checks") == expected.get("checks")
            ),
            "status_is_conditional": certificate.get("status")
            == "conditional_finite_selector_theorem",
            "all_examples_have_total_charge_12": all(
                item.get("total_charge") == 12 for item in poly.values()
            ),
            "strict_unit_minimum_is_recorded": certificate.get(
                "strict_unit_defect_minimum", {}
            ).get("unit_defect_count")
            == 12,
            "d_optimal_selector_is_exact": certificate.get(
                "d_optimal_tomography", {}
            ).get("axis_squared_inner_product")
            == "1/5",
            "normalized_negative_graphs_are_c5": (
                graph_class.get("distinct_degree_profiles") == [[2, 2, 2, 2, 2]]
                and graph_class.get("every_solution_connected") is True
                and graph_class.get("derived_graph_type") == "C5"
            ),
            "normalized_solutions_exhaust_labelled_c5_graphs": (
                graph_class.get("normalized_solution_count") == 12
                and graph_class.get("derived_labelled_c5_graph_count") == 12
                and graph_class.get("solutions_equal_all_labelled_c5_graphs") is True
            ),
            "seidel_class_count_is_derived": (
                seidel.get("switching_isomorphism_classes")
                == graph_class.get("isomorphism_class_count")
                == len(graph_class.get("canonical_edge_codes", []))
                == 1
            ),
            "orbit_size_is_12": certificate.get("orbit_stabilizer", {}).get("orbit_size") == 12,
            "screen_theorem_output_is_x_over_12": arithmetic.get("local_port_read")
            == "X/12",
            "gamma_screen_is_imported_algebra": (
                arithmetic.get("imported_screen_load") == "X=log(N/pi)"
                and arithmetic.get("imported_cell_entropy") == "P/4"
                and arithmetic.get("imported_beta_EW") == 4
                and arithmetic.get("gamma_screen_definition")
                == "Gamma_screen=beta_EW*(P/4)*(X/12)"
                and arithmetic.get("gamma_screen_simplified")
                == "(P/12)*log(N/pi)"
            ),
            "named_hierarchy_readout_gate_is_open": (
                readout_gate.get("premise_id") == "HIERARCHY-SCREEN-READOUT"
                and readout_gate.get("supplied_by_screen_sieve") is False
                and readout_gate.get("required_identification")
                == "log(E_cell/v)=Gamma_screen"
                and readout_gate.get("status") == "work_in_progress"
            ),
            "screen_certificate_has_no_physical_hierarchy_readout": (
                "hierarchy_readout" not in arithmetic
                and "alpha_U" not in arithmetic
                and "B_EW" not in arithmetic
            ),
        }
    return checks


def main(path: str) -> int:
    cert_path = pathlib.Path(path)
    certificate = json.loads(cert_path.read_text(encoding="utf-8"))
    checks = validate(certificate)
    payload = {"checks": checks, "pass": all(checks.values())}
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0 if payload["pass"] else 1


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: validate_screen_sieve_certificate.py CERTIFICATE", file=sys.stderr)
        raise SystemExit(2)
    raise SystemExit(main(sys.argv[1]))
