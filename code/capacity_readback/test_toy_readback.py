#!/usr/bin/env python3
"""Tests for the toy capacity-readback schema demonstration.

The tests exercise the certificate machinery on the declared toy model only.
Nothing here evaluates or constrains the physical readback map F; CL-7 is open.
"""

from __future__ import annotations

import json
from pathlib import Path
import unittest

from mpmath import iv, mp, mpf

from toy_readback import (
    ARTIFACT_NAME,
    build_certificate,
    enumerate_sector_count,
    toy_nf,
    toy_readback,
    toy_sector_count_closed_form,
)

RUNTIME_CERTIFICATE = Path(__file__).resolve().parent / "runtime" / "toy_readback_certificate.json"

REQUIRED_TOP_LEVEL_KEYS = {
    "artifact",
    "disclaimer",
    "physical_content",
    "moves_cl7",
    "cl7_status",
    "specification",
    "variant",
    "status",
    "interval_backend",
    "toy_model",
    "enumeration_check",
    "contraction_certificate",
    "fixed_point",
    "count_density",
    "blindness",
    "promotion_allowed",
}

REQUIRED_CONTRACTION_KEYS = {
    "interval",
    "image",
    "self_map_pass",
    "derivative_enclosure",
    "lipschitz_bound_L",
    "lipschitz_pass",
    "monotone_nonnegative_pass",
    "banach_pass",
}


def build_default_certificate() -> dict:
    return build_certificate(
        variant="toy",
        interval_lo="40",
        interval_hi="50",
        precision=40,
        enumeration_min=3,
        enumeration_max=8,
        enclosure_half_width="1e-25",
    )


class FinitePipelineTests(unittest.TestCase):
    def test_nf_is_idempotent_and_order_independent(self) -> None:
        self.assertEqual(toy_nf("210"), "012")
        self.assertEqual(toy_nf(toy_nf("210")), toy_nf("210"))
        self.assertEqual(toy_nf("120"), toy_nf("021"))

    def test_enumerated_sector_matches_closed_form(self) -> None:
        for n in range(3, 9):
            self.assertEqual(enumerate_sector_count(n), toy_sector_count_closed_form(n))


class ContractionCertificateTests(unittest.TestCase):
    def test_toy_fixed_point_exists_and_certificate_passes(self) -> None:
        certificate = build_default_certificate()

        self.assertEqual(certificate["status"], "pass")
        contraction = certificate["contraction_certificate"]
        self.assertTrue(contraction["self_map_pass"])
        self.assertTrue(contraction["lipschitz_pass"])
        self.assertTrue(contraction["banach_pass"])
        self.assertLess(mpf(contraction["lipschitz_bound_L"]), 1)

        fixed_point = certificate["fixed_point"]
        self.assertTrue(fixed_point["located"])
        self.assertTrue(fixed_point["box_self_map_pass"])
        lo = mpf(fixed_point["enclosure"]["lo"])
        hi = mpf(fixed_point["enclosure"]["hi"])
        self.assertLess(lo, hi)
        self.assertLess(hi - lo, mpf("1e-20"))

    def test_fixed_point_enclosure_solves_the_readback_equation(self) -> None:
        certificate = build_default_certificate()
        iv.dps = 40
        mp.dps = 40
        enclosure = certificate["fixed_point"]["enclosure"]
        box = iv.mpf([enclosure["lo"], enclosure["hi"]])
        image = toy_readback(box)
        # F(box) and box overlap: |F(mid) - mid| below the box width scale.
        residual = abs(mpf(image.mid) - mpf(box.mid))
        self.assertLess(residual, mpf("1e-20"))

    def test_certificate_schema_keys_present(self) -> None:
        certificate = build_default_certificate()

        self.assertEqual(certificate["artifact"], ARTIFACT_NAME)
        self.assertEqual(REQUIRED_TOP_LEVEL_KEYS - set(certificate), set())
        self.assertEqual(REQUIRED_CONTRACTION_KEYS - set(certificate["contraction_certificate"]), set())
        for key in ("inputs", "reads_measured_lambda", "reads_sl4_estimate", "dependency_cone"):
            self.assertIn(key, certificate["blindness"])

    def test_toy_labeling_is_explicit(self) -> None:
        certificate = build_default_certificate()

        self.assertFalse(certificate["physical_content"])
        self.assertFalse(certificate["moves_cl7"])
        self.assertFalse(certificate["promotion_allowed"])
        self.assertEqual(certificate["cl7_status"], "open")
        self.assertIn("toy model", certificate["disclaimer"])
        self.assertFalse(certificate["blindness"]["reads_measured_lambda"])
        self.assertFalse(certificate["blindness"]["reads_sl4_estimate"])

    def test_non_contracting_variant_is_rejected(self) -> None:
        certificate = build_certificate(
            variant="non_contracting",
            interval_lo="40",
            interval_hi="50",
            precision=40,
            enumeration_min=3,
            enumeration_max=6,
            enclosure_half_width="1e-25",
        )

        self.assertEqual(certificate["status"], "rejected")
        contraction = certificate["contraction_certificate"]
        self.assertFalse(contraction["banach_pass"])
        self.assertFalse(contraction["lipschitz_pass"])
        self.assertGreaterEqual(mpf(contraction["lipschitz_bound_L"]), 1)
        self.assertFalse(certificate["fixed_point"]["located"])

    def test_determinism(self) -> None:
        first = json.dumps(build_default_certificate(), sort_keys=True)
        second = json.dumps(build_default_certificate(), sort_keys=True)
        self.assertEqual(first, second)


class RuntimeArtifactTests(unittest.TestCase):
    def test_runtime_certificate_is_current_and_labeled(self) -> None:
        if not RUNTIME_CERTIFICATE.exists():
            self.skipTest("runtime certificate not generated")
        certificate = json.loads(RUNTIME_CERTIFICATE.read_text(encoding="utf-8"))
        self.assertEqual(certificate["artifact"], ARTIFACT_NAME)
        self.assertFalse(certificate["physical_content"])
        self.assertFalse(certificate["moves_cl7"])
        self.assertEqual(certificate["cl7_status"], "open")
        self.assertEqual(REQUIRED_TOP_LEVEL_KEYS - set(certificate), set())


if __name__ == "__main__":
    unittest.main()
