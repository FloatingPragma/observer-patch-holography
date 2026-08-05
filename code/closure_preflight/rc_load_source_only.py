#!/usr/bin/env python3
"""Recompute the source-only RC-LOAD fixed-point packet.

This is the auditable RER projection of the workspace-root
``n_closure_campaign/closure_certificate.py`` calculation.  It reads the
certified P and alpha_U intervals from the RER P certificate and proves the
declared RC-LOAD map is a contraction on its declared macroscopic domain.
The replay deliberately carries no observational comparison values.

The interval theorem is conditional on the RC-LOAD pricing law.  Whether that
law is source-selected and whether its two readings denote one quantity are
separate H0 gates recorded by the upstream provenance projection.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

from mpmath import iv, mp, mpf


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
P_INTERVAL = (
    ROOT
    / "code/P_derivation/runtime/p_interval_contraction_certificate_2026-07-14.json"
)
PROVENANCE = HERE / "data/rc_load_workspace_provenance_v1.json"
OUTPUT = HERE / "outputs/rc_load_source_only_certificate.json"

SCHEMA = "oph.rc_load_source_only_replay.v1"
PROVENANCE_SCHEMA = "oph.rc_load_workspace_provenance.v1"
MODE = "thomson_structured_running"
PRECISION = 80


def canonical_json_bytes(value: Any) -> bytes:
    return (
        json.dumps(value, ensure_ascii=True, separators=(",", ":"), sort_keys=True)
        + "\n"
    ).encode("ascii")


def tagged(raw: bytes) -> str:
    return "sha256:" + hashlib.sha256(raw).hexdigest()


def load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def pin(path: Path) -> dict[str, Any]:
    raw = path.read_bytes()
    return {
        "path": path.relative_to(ROOT).as_posix(),
        "bytes": len(raw),
        "sha256": tagged(raw),
    }


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _endpoint(value: Any) -> str:
    """Return the scalar endpoint text used by mpmath interval values."""

    text = str(value)
    if text.startswith("[") and text.endswith("]"):
        left, right = text[1:-1].split(",", 1)
        require(left.strip() == right.strip(), "interval endpoint is not scalar")
        return left.strip()
    return text


def _enclosure(value: Any) -> dict[str, str]:
    return {"lo": _endpoint(value.a), "hi": _endpoint(value.b)}


def build() -> dict[str, Any]:
    mp.dps = PRECISION
    iv.dps = PRECISION
    p_certificate = load(P_INTERVAL)
    provenance = load(PROVENANCE)
    require(
        provenance.get("schema") == PROVENANCE_SCHEMA,
        "RC-LOAD provenance schema drift",
    )
    upstream = provenance.get("certificate", {})
    require(
        upstream.get("equation")
        == "ln(N/pi) * (6*pi/(P*alpha_U) - ln(N/pi)) = 6*pi",
        "RC-LOAD equation drift",
    )
    require(upstream.get("verdict") == "PASS", "upstream RC-LOAD certificate failed")

    mode = p_certificate.get("modes", {}).get(MODE, {})
    p_box = mode.get("certified_enclosure", {}).get("P", {})
    au_box = mode.get("interval_diagnostics", {}).get("alpha_u", {})
    for label, box in (("P", p_box), ("alpha_U", au_box)):
        require(set(box) >= {"lo", "hi"}, f"missing {label} enclosure")
        require(mpf(box["lo"]) < mpf(box["hi"]), f"invalid {label} enclosure")

    p_value = iv.mpf([p_box["lo"], p_box["hi"]])
    alpha_u = iv.mpf([au_box["lo"], au_box["hi"]])
    pi_interval = iv.pi

    budget_x = 6 * pi_interval / (p_value * alpha_u)
    discriminant = budget_x**2 - 24 * pi_interval
    require(discriminant.a > 0, "RC-LOAD discriminant is not positive")
    root = (budget_x + iv.sqrt(discriminant)) / 2
    micro_root = (budget_x - iv.sqrt(discriminant)) / 2

    domain_lo = mpf(10)
    domain_point = iv.mpf([domain_lo, domain_lo])
    image_at_lo = budget_x - 6 * pi_interval / domain_point
    maps_into = image_at_lo.a > domain_lo and budget_x.b > image_at_lo.b
    contraction_lipschitz = (6 * pi_interval / domain_point**2).b
    local_lipschitz = (6 * pi_interval / root**2).b
    contraction = contraction_lipschitz < 1
    residual = root * (budget_x - root) - 6 * pi_interval
    residual_contains_zero = residual.a <= 0 <= residual.b
    micro_outside = micro_root.b < domain_lo

    checks = {
        "positive_discriminant": bool(discriminant.a > 0),
        "maps_into_domain": bool(maps_into),
        "contraction": bool(contraction),
        "unique_fixed_point_in_domain": bool(maps_into and contraction),
        "residual_encloses_zero": bool(residual_contains_zero),
        "micro_root_outside_domain": bool(micro_outside),
    }
    verdict = "PASS" if all(checks.values()) else "FAIL"

    value: dict[str, Any] = {
        "schema": SCHEMA,
        "artifact": "oph_rc_load_source_only_replay",
        "status": (
            "CONDITIONAL_FIXED_POINT_CERTIFIED"
            if verdict == "PASS"
            else "INTERVAL_CERTIFICATE_FAILED"
        ),
        "equation": upstream["equation"],
        "scope": {
            "pricing_law_proved_here": upstream["pricing_law_proved_here"],
            "same_quantity_identity_proved_here": upstream[
                "typed_identity_proved_here"
            ],
            "observational_comparison_imported": False,
            "physical_or_laboratory_attachment_evaluated": False,
            "workspace_root_checkout_required": provenance[
                "external_files_required_for_replay"
            ],
        },
        "exposure": {
            "class": provenance["campaign"]["exposure_class"],
            "promotion_from_landing": provenance["campaign"][
                "promotion_from_landing"
            ],
        },
        "inputs": {
            "P": {"lo": p_box["lo"], "hi": p_box["hi"]},
            "alpha_U": {"lo": au_box["lo"], "hi": au_box["hi"]},
            "P_mode": MODE,
        },
        "budget_x": _enclosure(budget_x),
        "banach": {
            "domain_lo": "10",
            "domain_hi": _endpoint(budget_x.b),
            "maps_into_domain": bool(maps_into),
            "sup_lipschitz": _endpoint(contraction_lipschitz),
            "local_lipschitz_at_fixed_point": _endpoint(local_lipschitz),
            "unique_fixed_point_in_domain": bool(maps_into and contraction),
        },
        "fixed_point": {
            "X": _enclosure(root),
            "residual": _enclosure(residual),
            "residual_encloses_zero": bool(residual_contains_zero),
        },
        "excluded_micro_root": {
            "X": _enclosure(micro_root),
            "outside_declared_domain": bool(micro_outside),
        },
        "validation": checks,
        "verdict": verdict,
        "parent_pins": [pin(P_INTERVAL), pin(PROVENANCE)],
        "upstream_workspace_provenance": {
            "source_sha256": provenance["source"]["sha256"],
            "certificate_sha256": provenance["certificate"]["sha256"],
            "campaign_sha256": provenance["campaign"]["sha256"],
        },
    }
    value["certificate_sha256"] = tagged(canonical_json_bytes(value))
    return value


def verify(path: Path = OUTPUT) -> None:
    expected = build()
    actual = load(path)
    if actual != expected:
        raise SystemExit("RC-LOAD source-only replay is stale")
    if path.read_bytes() != canonical_json_bytes(expected):
        raise SystemExit("RC-LOAD source-only replay is not canonical")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.write:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_bytes(canonical_json_bytes(build()))
    if args.verify:
        verify()
        print("RC_LOAD_SOURCE_ONLY_VALID")
    if not args.write and not args.verify:
        print(canonical_json_bytes(build()).decode("ascii"), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
