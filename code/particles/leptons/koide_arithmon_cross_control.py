#!/usr/bin/env python3
"""Pinned external Arithmon negative control for the exact OPH Koide balance.

This is a compare-only same-quantity test, not a prediction freeze. OPH proves
that its nonnegative balanced circulant has Q = 2/3 (with the positive chamber
needed for the physical square-root-mass reading). Independently,
Arithmon/K7-Lean proves Q(27^phi, 3477) < 2/3 for its public mass-pair formulas.
The two receipts therefore cannot be realized simultaneously as one exact
balanced positive-chamber charged-family spectrum.

The exact theorems remain kernel-checked in separate repositories. Nothing in
this file imports Arithmon as an OPH axiom, derives the charged-family
attachment, or adjudicates either framework physically.
"""

from __future__ import annotations

import argparse
import json
from decimal import Decimal, localcontext
from pathlib import Path
from typing import Any

CODE_ROOT = Path(__file__).resolve().parents[2]
REPO_ROOT = CODE_ROOT.parent
FIXTURE = CODE_ROOT / "particles" / "data" / "arithmon_koide_cross_0.json"
OUT = CODE_ROOT / "particles" / "runs" / "leptons" / "koide_arithmon_cross_0.json"
OPH_KOIDE = REPO_ROOT / "Lean" / "ObserverPatchHolography" / "KoideCirculant.lean"

PRECISION = 120
EXPECTED = {
    "fixture_schema": "oph.external_arithmon_koide_cross_fixture.v1",
    "k7_commit": "af9224cf3bd0c7bb8a895b98d4d1146eb138b555",
    "k7_blob": "646500c22e40fcd89d3584df7d64d5a9d452f754",
    "lean_commit": "836296db191c49307f7df59aa0df989399a0b1bb",
    "lean_blob": "47ff66d563e50f64f73ad085490459f28ff0df85",
    "lean_theorem": "GIFT.Relations.KoideAssembly.koideQ_gift_lt_two_thirds",
}

ONE, TWO, THREE = Decimal(1), Decimal(2), Decimal(3)


class CertificateError(ValueError):
    pass


def require(ok: bool, code: str) -> None:
    if not ok:
        raise CertificateError(code)


def load_fixture() -> dict[str, Any]:
    f = json.loads(FIXTURE.read_text(encoding="utf-8"))
    require(f.get("schema") == EXPECTED["fixture_schema"], "FIXTURE_SCHEMA")
    k7, lean = f["freeze"]["arithmon_k7"], f["freeze"]["arithmon_k7_lean"]
    require(
        (k7["commit"], k7["git_blob_sha"])
        == (EXPECTED["k7_commit"], EXPECTED["k7_blob"]),
        "K7_SOURCE_PIN",
    )
    require(
        (lean["commit"], lean["git_blob_sha"], lean["theorem"])
        == (EXPECTED["lean_commit"], EXPECTED["lean_blob"], EXPECTED["lean_theorem"]),
        "K7_LEAN_SOURCE_PIN",
    )
    pair = f["version_pinned_arithmon_mass_pair"]
    require(pair["m_mu_over_m_e"]["formula"] == "27^phi", "MUON_FORMULA")
    require(pair["m_tau_over_m_e"]["formula"] == "3477", "TAU_FORMULA")
    require(f.get("experimental_data_consumed") is False, "EXPERIMENTAL_INPUT")
    require(f.get("prediction_freeze_claimed") is False, "PREDICTION_FREEZE_SCOPE")
    oph = OPH_KOIDE.read_text(encoding="utf-8")
    require(
        "theorem koide_eq_two_thirds_iff" in oph
        and "ratio = 1 / Real.sqrt 2" in oph,
        "OPH_KOIDE_THEOREM",
    )
    return f


def koide_q(m_e: Decimal, m_mu: Decimal, m_tau: Decimal) -> Decimal:
    roots = m_e.sqrt() + m_mu.sqrt() + m_tau.sqrt()
    return (m_e + m_mu + m_tau) / roots**2


def balance_roots(m1: Decimal, m2: Decimal) -> tuple[Decimal, Decimal]:
    p, s = m1.sqrt() + m2.sqrt(), m1 + m2
    d = Decimal(6) * p**2 - THREE * s
    require(d > 0, "BALANCE_DISCRIMINANT")
    r = d.sqrt()
    return (TWO * p + r) ** 2, (TWO * p - r) ** 2


def build_payload() -> dict[str, Any]:
    fixture = load_fixture()
    with localcontext() as ctx:
        ctx.prec = PRECISION
        phi = (ONE + Decimal(5).sqrt()) / TWO
        mu = (phi * Decimal(27).ln()).exp()
        tau = Decimal(3477)
        target = TWO / THREE
        q = koide_q(ONE, mu, tau)

        tau_plus, tau_minus = balance_roots(ONE, mu)
        mu_plus, mu_minus = balance_roots(ONE, tau)
        require(q < target, "LOCAL_SIGN_MISMATCH")
        require(tau_plus > mu > ONE and tau_minus < mu, "TAU_ROOT_ORDER")
        require(ONE < mu_minus < tau and mu_plus > tau, "MU_ROOT_ORDER")

        q_tau = koide_q(ONE, mu, tau_plus)
        q_mu = koide_q(ONE, mu_minus, tau)
        require(abs(q_tau - target) < Decimal("1e-110"), "TAU_REPAIR_CONTROL")
        require(abs(q_mu - target) < Decimal("1e-110"), "MU_REPAIR_CONTROL")

        txt = lambda x: format(+x, "f")
        return {
            "artifact": "oph_arithmon_koide_cross_0",
            "schema": "oph.koide_arithmon_cross_0.v1",
            "status": "EXTERNAL_NEGATIVE_CONTROL__ARITHMON_MASS_PAIR_INCOMPATIBLE_WITH_EXACT_OPH_BALANCE",
            "claim_class": "external_negative_control",
            "issue_context": [546, 736],
            "source_only_physical_prediction": False,
            "prospective_prediction_freeze": False,
            "public_physical_promotion_allowed": False,
            "formal_composition_in_one_kernel": False,
            "experimental_data_consumed": False,
            "claim_boundary": (
                "Compare-only same-quantity compatibility test. OPH and Arithmon exact receipts "
                "are kernel-checked in separate repositories and version-pinned here, but not "
                "composed in one Lean environment. No Arithmon statement is imported as an OPH "
                "premise; no charged-family attachment or physical adjudication is claimed."
            ),
            "same_quantity_bridge": {
                "quantity": "Q = (m_e + m_mu + m_tau) / (sqrt(m_e) + sqrt(m_mu) + sqrt(m_tau))^2",
                "normalization": "m_e = 1 (scale invariance)",
                "oph_scope": "positive-chamber exact balanced circulant implies Q = 2/3",
            },
            "pinned_receipts": {
                "oph": {
                    "path": "Lean/ObserverPatchHolography/KoideCirculant.lean",
                    "theorem": "OPH.KoideCirculant.koide_eq_two_thirds_iff",
                },
                "arithmon": fixture["freeze"],
            },
            "version_pinned_arithmon_inputs": fixture["version_pinned_arithmon_mass_pair"],
            "local_recompute": {
                "decimal_precision": PRECISION,
                "m_mu_over_m_e": txt(mu),
                "m_tau_over_m_e": txt(tau),
                "Q": txt(q),
                "Q_minus_two_thirds": txt(q - target),
                "strictly_below_two_thirds": True,
            },
            "logical_verdict": {
                "oph_exact_balance_requires": "Q = 2/3",
                "arithmon_exact_external_receipt": "Q(27^phi, 3477) < 2/3",
                "joint_exact_realization_possible": False,
                "scope": "version-pinned Arithmon mass pair versus exact balanced positive-chamber OPH spectrum",
            },
            "controls": {
                "replace_tau_by_balance_root": {
                    "replacement_m_tau_over_m_e": txt(tau_plus),
                    "delta": txt(tau_plus - tau),
                    "abs_Q_minus_two_thirds_after_replacement": txt(abs(q_tau - target)),
                    "restores_exact_balance_numerically": True,
                },
                "replace_muon_by_balance_root": {
                    "replacement_m_mu_over_m_e": txt(mu_minus),
                    "delta": txt(mu_minus - mu),
                    "abs_Q_minus_two_thirds_after_replacement": txt(abs(q_mu - target)),
                    "restores_exact_balance_numerically": True,
                },
            },
            "checks": {
                "external_source_pins_match_fixture": True,
                "current_oph_koide_theorem_surface_present": True,
                "no_experimental_mass_input": True,
                "local_sign_matches_external_exact_receipt": True,
                "balance_root_controls_close": True,
            },
            "checks_pass": True,
        }


def canonical(payload: dict[str, Any]) -> bytes:
    return (json.dumps(payload, indent=2, ensure_ascii=False) + "\n").encode()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, default=OUT)
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()
    expected = canonical(build_payload())
    if args.check:
        require(args.output.is_file(), "OUTPUT_MISSING")
        require(args.output.read_bytes() == expected, "OUTPUT_STALE")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_bytes(expected)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
