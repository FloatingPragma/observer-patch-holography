from __future__ import annotations

import argparse
import sys
import unittest
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import settled_domain_applicability as sda  # noqa: E402


def payload() -> dict:
    return sda.compute(argparse.Namespace(n_scr=sda.DEFAULT_N_SCR))


class SettledDomainApplicabilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = payload()

    def test_declared_scales_match_paper_values(self) -> None:
        scales = self.payload["declared_scales"]
        self.assertAlmostEqual(
            scales["a0_oph_m_s2"], 1.029186271e-10, delta=1e-18
        )
        self.assertAlmostEqual(
            scales["lambda_collar"], 0.932042991, places=9
        )
        self.assertAlmostEqual(
            scales["a_eff_m_s2"], 1.184737388e-10, delta=1e-18
        )

    def test_route_one_is_non_discriminating(self) -> None:
        record = self.payload["route_one_record"]
        self.assertGreater(record["c_over_a_eff_gyr"], 13.8)
        self.assertLess(record["max_repair_cycles_since_big_bang"], 1.0)
        self.assertFalse(record["discriminating"])

    def test_cassini_gate_exempt_with_margin(self) -> None:
        cas = self.payload["gates"]["cassini"]
        self.assertTrue(cas["verdict_exempt"])
        self.assertGreater(cas["lambda_sqrt_x"], 700.0)
        self.assertLessEqual(cas["nu_minus_one_bound"], 1e-320)
        self.assertLess(cas["pull_sigma"], 1.0)
        self.assertGreaterEqual(
            cas["suppression_vs_universal"], cas["suppression_floor"]
        )
        # The transition zone sits inside the Sun's Galactic Jacobi radius,
        # so the exemption rests on the enclosed-flux clause alone.
        self.assertTrue(cas["transition_zone_inside_jacobi_radius"])
        self.assertLess(cas["q2_predicted_bound_s2"], cas["q2_cassini_sigma_s2"])

    def test_sparc_gate_retained(self) -> None:
        spc = self.payload["gates"]["sparc"]
        self.assertTrue(spc["verdict_retained"])
        self.assertLess(spc["deep_limit_rel_err_at_x_1e-8"], 1e-3)
        self.assertLess(spc["newtonian_limit_deviation_at_x_1e6"], 1e-6)
        self.assertGreater(spc["clock_ratio_ext_over_int"], 1.0)
        self.assertAlmostEqual(
            spc["tau_j_internal_outer_disk_gyr"], 0.316276, places=4
        )

    def test_wide_binary_stance_emitted_with_sign(self) -> None:
        wb = self.payload["gates"]["wide_binaries"]
        self.assertTrue(wb["verdict_emitted"])
        self.assertEqual(
            wb["stance"], "anomalous_full_boost_no_external_field_effect"
        )
        self.assertGreater(wb["boost_at_sample"], 1.25)
        self.assertLess(wb["boost_at_sample"], 1.40)
        self.assertGreater(wb["boost_at_g_1e-11"], wb["aqual_efe_boost_ceiling"])
        self.assertGreater(wb["clock_ratio_ext_over_int"], 1.0)

    def test_satellite_rows_split_by_tidal_clock(self) -> None:
        rows = {
            s["name"]: s
            for s in self.payload["gates"]["wide_binaries"]["satellites"]
        }
        self.assertFalse(rows["Crater II"]["settled"])
        self.assertTrue(rows["Fornax"]["settled"])

    def test_all_gates_pass(self) -> None:
        self.assertTrue(self.payload["all_gates_pass"])


if __name__ == "__main__":
    unittest.main()
