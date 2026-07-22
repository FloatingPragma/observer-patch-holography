#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
import json
import math
import pathlib
import sys
import unittest

# Layout note: this repo keeps the certificates flat alongside the suite,
# rather than the upstream bundle's code/ + tests/ split.
ROOT = pathlib.Path(__file__).resolve().parent


def load(name: str):
    path = ROOT / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    sys.modules[name] = module
    try:
        spec.loader.exec_module(module)
    except Exception:
        sys.modules.pop(name, None)
        raise
    return module


class AuditTests(unittest.TestCase):
    def test_gravity_units(self):
        m = load("gravity_units_corrected")
        p = m.evaluate(m.ARCHIVED_G_N)
        self.assertAlmostEqual(p.acceleration_boost_nu, 1.725381227058617, places=12)
        self.assertAlmostEqual(p.velocity_boost_sqrt_nu, 1.3135376762996245, places=12)
        self.assertAlmostEqual(p.gamma_log10_velocity_boost, 0.11844253420192026, places=12)

    def test_a5_decompositions(self):
        m = load("a5_harmonic_decomposition")
        self.assertEqual(m.decompose(3), {"3prime": 1, "4": 1})
        self.assertEqual(m.decompose(5), {"3": 1, "3prime": 1, "5": 1})
        self.assertEqual(m.decompose(6), {"1": 1, "3": 1, "4": 1, "5": 1})
        self.assertTrue(m.payload(15)["port_equals_restriction_H0_plus_H5"])

    def test_a5_selection_certificate(self):
        m = load("a5_selection_certificate")
        p = m.payload()
        self.assertEqual(p["number_of_distances"], 3)
        self.assertEqual(p["spherical_design_strength"], 5)
        self.assertTrue(p["fourth_moment_matches_uniform_S2"])
        self.assertEqual(p["d_optimal_determinants"], ["8", "1024/3125"])
        self.assertEqual(p["switched_seidel_solutions"], 12)

    def test_compact_dimensions(self):
        m = load("a5_compact_lie_classifier")
        self.assertEqual(m.partitions(11), [(3, 8)])
        self.assertEqual(m.partitions(12), [(3, 3, 3, 3)])
        p = m.payload()
        self.assertIn(
            "counterbranches_if_innerness_is_dropped_but_group_integrality_is_retained",
            p,
        )
        self.assertIn(
            "additional_real_module_counterbranch_only_if_group_level_integral_central_lattice_constraints_are_also_dropped",
            p,
        )

    def test_a5_claim_boundary(self):
        m = load("a5_screen_sm_closure")
        state = m.payload["theorem_scope"]
        self.assertIn("coefficient", state["exact_here"])
        self.assertTrue(any("PORT-CURRENT-INNER" in x for x in state["physical_gates"]))
        self.assertTrue(any("MAR" in x for x in state["physical_gates"]))
        port_gate = next(x for x in state["physical_gates"] if "PORT-CURRENT-INNER" in x)
        self.assertIn("conditional", port_gate)
        self.assertIn("remain open", port_gate)
        self.assertNotIn("is closed", port_gate)

        receipt = json.loads(
            (ROOT / "receipts" / "port_current_inner_reference.receipt.json").read_text()
        )
        self.assertTrue(receipt["conditional_algebraic_gate"]["passed"])
        self.assertFalse(receipt["physical_source_gate"]["passed"])
        self.assertFalse(receipt["issue_closure_condition"]["met_locally"])

        registry = json.loads((ROOT.parent.parent / "claims" / "claim_registry.yaml").read_text())
        claim = next(
            row for row in registry["claims"]
            if row["claim_id"] == "OPH-SCREEN-PORT-CURRENT-INNER"
        )
        self.assertEqual(
            claim["status"],
            receipt["claim_boundary"]["status"],
        )

        matter_gate = next(x for x in state["physical_gates"] if "matter lift" in x)
        self.assertIn("conditional", matter_gate)
        self.assertIn("remains open", matter_gate)
        self.assertNotIn("is closed", matter_gate)

        matter_receipt = json.loads(
            (ROOT / "receipts" / "super_tannakian_matter_reference.receipt.json").read_text()
        )
        self.assertTrue(matter_receipt["conditional_algebraic_gate"]["passed"])
        self.assertFalse(matter_receipt["physical_source_gate"]["passed"])
        self.assertFalse(matter_receipt["issue_closure_condition"]["met_locally"])

        matter_claim = next(
            row for row in registry["claims"]
            if row["claim_id"] == "OPH-SCREEN-SUPER-TANNAKIAN-MATTER-LIFT"
        )
        self.assertEqual(
            matter_claim["status"],
            matter_receipt["claim_boundary"]["status"],
        )

    def test_exterior_sm_completion(self):
        m = load("exterior_sm_completion")
        p = m.payload
        self.assertEqual(p["matter_package"]["complex_dimension"], 15)
        self.assertEqual(set(p["anomalies"]["coefficients"].values()), {"0"})
        self.assertIn("Witten", p["anomalies"]["su2_witten_parity"])
        self.assertTrue(all(
            row["su3_singlet_multiplicity"] == 1
            and row["su2_singlet_multiplicity"] == 1
            for row in p["yukawa_channels"].values()
        ))
        self.assertEqual(p["weak_load"]["doublet_multiplicity_per_generation"], 4)
        self.assertEqual(p["face_phase"]["minimal_irreps"], ["3", "3prime"])
        self.assertIn("General family matrices", p["face_phase"]["yukawa_boundary"])
        self.assertEqual(p["abstract_deck_control"]["orbit_count"], 4)
        self.assertIn(
            "not the full even Clifford module",
            p["conditional_closure"]["no_extra_sector_boundary"],
        )
        self.assertIn(
            "not Clifford-stable",
            p["conditional_closure"]["no_extra_sector_boundary"],
        )
        runtime = json.loads((ROOT / "exterior_sm_completion.json").read_text())
        self.assertEqual(runtime, p)

    def test_survival_boundary_certificates(self):
        m = load("survival_boundary_certificates")
        no_go = m.source_completion_nonuniqueness_certificate()
        self.assertEqual(no_go["status"], "EXACT_FINITE_NONIDENTIFIABILITY_THEOREM")
        self.assertEqual(len(no_go["inequivalent_current_completions"]), 2)
        self.assertEqual(len(no_go["inequivalent_matter_completions"]), 2)
        self.assertFalse(no_go["source_only_reconstruction_of_every_completion"])

        p15 = m.rank15_clifford_certificate()
        self.assertEqual(p15["majorana_count"], 10)
        self.assertEqual(p15["rank"], 15)
        self.assertFalse(p15["physical_promotion"])

        gaussian = m.gaussian_composite_1pi_certificate()
        self.assertEqual(gaussian["connected_third_cumulant"], "14")
        self.assertEqual(
            gaussian["composite_effective_action_third_derivative"], "-14/27"
        )
        self.assertEqual(gaussian["fundamental_effective_action_third_derivative"], "0")

        refinement = m.refinement_complement_certificate()
        self.assertTrue(refinement["positive"]["passes_complement_interior_exclusion"])
        self.assertFalse(
            refinement["hidden_zero_mode"]["passes_complement_interior_exclusion"]
        )
        self.assertEqual(refinement["hidden_zero_mode"]["fine_selected_rank"], 4)

        settlement = m.finite_settlement_certificate()
        self.assertEqual(settlement["risk_path"], [16, 16, 12])
        self.assertFalse(settlement["physical_promotion"])

    def test_log_coefficients(self):
        m = load("bh_log_correction")
        p = m.payload()
        self.assertAlmostEqual(p["candidate_coefficients"]["12_port_exact_balance"], 5.5)
        self.assertLess(abs(p["numerical_stirling_checks"]["q12_effective_c_12000_to_120000"] - 5.5), 0.01)

    def test_jacobi_scaling(self):
        m = load("hysteresis_scaling_guard")
        p = m.payload()["worked_ratio_50_to_100_kpc"]
        self.assertAlmostEqual(p["Gamma_50_over_Gamma_100"], 2 ** 1.5, places=12)
        self.assertAlmostEqual(p["Gamma2_50_over_Gamma2_100"], 8.0, places=12)


if __name__ == "__main__":
    unittest.main(verbosity=2)
