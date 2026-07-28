#!/usr/bin/env python3
"""Koide balance comparison and the conditional tau interval.

The positive-chamber circulant identity
(``extra/koide_identity_from_positive_c3_face_circulants.tex``, with the
finite tracial GNS balance as a conditional theorem) states, for a
positive-chamber face circulant with square-root masses
sqrt(m_k) = sqrt(s) (a + 2 rho cos(delta + 2 pi k/3)):

    Q := (sum m) / (sum sqrt m)^2 = 1/3 + (2/3)(rho/a)^2,
    Q = 2/3  exactly when  rho/a = 1/sqrt2.

This certificate carries two lanes, both at 100-decimal-digit arithmetic
with outward rounding and explicit corner enumeration.

Lane one (balance comparison, compare-only): the measured charged triple
gives an enclosure for Q; the certificate records the exact distance of
Q = 2/3 from that enclosure in units of the propagated half-width.

Lane two (conditional tau interval): under the balanced-circulant premise
and the mass-ordering premise, the identity determines sqrt(m_tau) from
(m_e, m_mu) as the admissible root of one quadratic.  The certificate
propagates the measured (m_e, m_mu) enclosures through the closed-form
root at every corner, rounds outward, and compares the resulting tau
interval with the measured tau mass.

ANCESTRY (declared, read before citing): the balance condition was first
abstracted from the measured charged triple (Koide, 1981-1983); the
circulant identity gives it a finite structural home and the finite-GNS
construction supplies the balanced modulus as a conditional theorem, but
no source-only derivation of the premise exists in the corpus.  Both
lanes are therefore conditional/compare-only rows with measured-target
ancestry on the premise class: lane two is a conditional postdiction of
one mass from two masses under a premise whose selection history includes
the third.  No lane is a prospective prediction, and the source-only
charged no-go is unchanged.
"""

from __future__ import annotations

import argparse
import json
from decimal import Decimal, localcontext
from functools import wraps
from pathlib import Path
from typing import Any, Callable, TypeVar

ROOT = Path(__file__).resolve().parents[2]
OUT_PATH = ROOT / "particles" / "runs" / "leptons" / "koide_balance_comparison.json"

SCHEMA = "oph.koide_balance_comparison_certificate.v1"
ISSUE_CONTEXT = [546]

DECIMAL_PRECISION = 120
DIGITS = 100

# PDG values in MeV with one-standard-deviation half-widths.
PDG_MEV = {
    "electron": (Decimal("0.51099895069"), Decimal("0.00000000016")),
    "muon": (Decimal("105.6583755"), Decimal("0.0000023")),
    "tau": (Decimal("1776.93"), Decimal("0.09")),
}
PDG_SOURCE = "PDG charged-lepton masses, compare-only measured imports"

TWO = Decimal(2)
THREE = Decimal(3)

ReturnT = TypeVar("ReturnT")


def high_precision(function: Callable[..., ReturnT]) -> Callable[..., ReturnT]:
    """Run certificate arithmetic without reading or mutating global context."""

    @wraps(function)
    def wrapped(*args: Any, **kwargs: Any) -> ReturnT:
        with localcontext() as context:
            context.prec = DECIMAL_PRECISION
            return function(*args, **kwargs)

    return wrapped


class CertificateError(ValueError):
    def __init__(self, code: str, message: str) -> None:
        super().__init__(f"{code}: {message}")
        self.code = code


def require(condition: bool, code: str, message: str) -> None:
    if not condition:
        raise CertificateError(code, message)


def dsqrt(value: Decimal) -> Decimal:
    require(value >= 0, "SQRT_DOMAIN", "square root of a negative enclosure corner")
    return value.sqrt()


@high_precision
def outward(lo: Decimal, hi: Decimal, places: int = 60) -> tuple[str, str]:
    """Round an enclosure outward at the stated decimal places."""

    require(lo <= hi, "ENCLOSURE_ORDER", "enclosure endpoints out of order")
    quantum = Decimal(1).scaleb(-places)
    lo_out = (lo / quantum).to_integral_value(rounding="ROUND_FLOOR") * quantum
    hi_out = (hi / quantum).to_integral_value(rounding="ROUND_CEILING") * quantum
    return str(lo_out), str(hi_out)


@high_precision
def koide_q(m_e: Decimal, m_mu: Decimal, m_tau: Decimal) -> Decimal:
    roots = dsqrt(m_e) + dsqrt(m_mu) + dsqrt(m_tau)
    return (m_e + m_mu + m_tau) / (roots * roots)


@high_precision
def balance_comparison() -> dict[str, Any]:
    """Lane one: the measured enclosure for Q against the exact 2/3."""

    corners: list[Decimal] = []
    for se in (-1, 1):
        for sm in (-1, 1):
            for st in (-1, 1):
                m_e = PDG_MEV["electron"][0] + se * PDG_MEV["electron"][1]
                m_mu = PDG_MEV["muon"][0] + sm * PDG_MEV["muon"][1]
                m_tau = PDG_MEV["tau"][0] + st * PDG_MEV["tau"][1]
                corners.append(koide_q(m_e, m_mu, m_tau))
    q_lo, q_hi = min(corners), max(corners)
    q_central = koide_q(*(PDG_MEV[k][0] for k in ("electron", "muon", "tau")))
    target = TWO / THREE
    half_width = (q_hi - q_lo) / TWO
    require(half_width > 0, "Q_WIDTH", "the Q enclosure must have positive width")
    distance_sigma = abs(target - q_central) / half_width
    inside = q_lo <= target <= q_hi
    lo_str, hi_str = outward(q_lo, q_hi, 15)
    return {
        "measured_Q_central": str(+q_central.quantize(Decimal(1).scaleb(-12))),
        "measured_Q_enclosure_outward": [lo_str, hi_str],
        "target": "2/3",
        "target_inside_enclosure": bool(inside),
        "distance_in_enclosure_half_widths": str(
            +distance_sigma.quantize(Decimal(1).scaleb(-4))
        ),
        "reading": (
            "the exact balance value 2/3 against the measured one-sigma "
            "corner enclosure of Q; a corner enclosure overcovers the "
            "one-sigma ellipsoid, so the distance figure is conservative"
        ),
    }


@high_precision
def conditional_tau_roots(m_e: Decimal, m_mu: Decimal) -> tuple[Decimal, Decimal]:
    """Solve Q(m_e, m_mu, x^2) = 2/3 for x = sqrt(m_tau), in closed form.

    With p = sqrt(m_e) + sqrt(m_mu) and s = m_e + m_mu the equation is
    x^2 - 4 p x + (3 s - 2 p^2) = 0, giving x = 2p ± sqrt(4 p^2 - 3 s + 2 p^2)
    = 2p ± sqrt(6 p^2 - 3 s).
    """

    p = dsqrt(m_e) + dsqrt(m_mu)
    s = m_e + m_mu
    disc = Decimal(6) * p * p - THREE * s
    require(disc > 0, "TAU_DISCRIMINANT", "the balanced quadratic must have real roots")
    root = dsqrt(disc)
    x_plus = TWO * p + root
    x_minus = TWO * p - root
    return x_plus * x_plus, x_minus * x_minus


@high_precision
def conditional_tau_interval() -> dict[str, Any]:
    """Lane two: the tau enclosure from (m_e, m_mu) under the balance."""

    plus_values: list[Decimal] = []
    minus_values: list[Decimal] = []
    for se in (-1, 1):
        for sm in (-1, 1):
            m_e = PDG_MEV["electron"][0] + se * PDG_MEV["electron"][1]
            m_mu = PDG_MEV["muon"][0] + sm * PDG_MEV["muon"][1]
            plus, minus = conditional_tau_roots(m_e, m_mu)
            plus_values.append(plus)
            minus_values.append(minus)
    tau_lo, tau_hi = min(plus_values), max(plus_values)
    lo_str, hi_str = outward(tau_lo, tau_hi, 6)
    central_plus, central_minus = conditional_tau_roots(
        PDG_MEV["electron"][0], PDG_MEV["muon"][0]
    )
    measured, sigma = PDG_MEV["tau"]
    distance_sigma = abs(central_plus - measured) / sigma
    inside_one_sigma = (measured - sigma) <= central_plus <= (measured + sigma)
    spurious = central_minus
    require(
        spurious < PDG_MEV["muon"][0],
        "ROOT_ORDER",
        "the minus root must violate the mass-ordering premise",
    )
    return {
        "premises": [
            "balanced positive-chamber circulant (rho/a = 1/sqrt2; the finite-GNS conditional theorem)",
            "mass ordering m_tau > m_mu (selects the plus root)",
        ],
        "inputs": "measured (m_e, m_mu) with one-sigma corners, declared imports",
        "tau_enclosure_mev_outward": [lo_str, hi_str],
        "tau_central_mev": str(+central_plus.quantize(Decimal(1).scaleb(-6))),
        "measured_tau_mev": [str(measured), str(sigma)],
        "distance_sigma": str(+distance_sigma.quantize(Decimal(1).scaleb(-4))),
        "measured_central_inside_one_sigma_of_output": bool(inside_one_sigma),
        "spurious_root_mev": str(+spurious.quantize(Decimal(1).scaleb(-6))),
        "spurious_root_excluded_by": "mass-ordering premise (below the muon mass)",
        "row_class": "conditional_postdiction_with_measured_premise_ancestry",
    }


# ---------------------------------------------------------------------------
# Controls
# ---------------------------------------------------------------------------


@high_precision
def control_premise_mutation() -> dict[str, Any]:
    """An unbalanced modulus must move the tau output far outside the
    measured band, so the balance premise is load-bearing and falsifiable."""

    # rho/a = 0.72 instead of 1/sqrt2: Q_target = 1/3 + 2/3 * 0.72^2.
    q_target = Decimal(1) / THREE + TWO / THREE * Decimal("0.5184")
    m_e = PDG_MEV["electron"][0]
    m_mu = PDG_MEV["muon"][0]
    p = dsqrt(m_e) + dsqrt(m_mu)
    s = m_e + m_mu
    # General solve: (1/Q_target)(s + x^2) = (p + x)^2.
    a_coeff = Decimal(1) / q_target - Decimal(1)
    b_coeff = -TWO * p
    c_coeff = s / q_target - p * p
    disc = b_coeff * b_coeff - Decimal(4) * a_coeff * c_coeff
    require(disc > 0, "CONTROL_DISC", "mutated quadratic must have real roots")
    x = (-b_coeff + dsqrt(disc)) / (TWO * a_coeff)
    tau_mut = x * x
    measured, sigma = PDG_MEV["tau"]
    shift_sigma = abs(tau_mut - measured) / sigma
    try:
        require(
            shift_sigma < Decimal(5),
            "PREMISE_MUTATION_DETECTED",
            "an unbalanced modulus moves tau far outside the measured band",
        )
    except CertificateError:
        return {
            "expected_failure": True,
            "failed": True,
            "code": "PREMISE_MUTATION_DETECTED",
            "mutated_tau_mev": str(+tau_mut.quantize(Decimal(1).scaleb(-3))),
            "shift_sigma": str(+shift_sigma.quantize(Decimal(1).scaleb(-1))),
        }
    return {"expected_failure": True, "failed": False}


@high_precision
def control_tau_not_consumed_in_lane_two() -> dict[str, Any]:
    """Lane two must not read the tau mass: recomputing it with a doctored
    tau value must produce the identical enclosure."""

    original = conditional_tau_interval()
    saved = PDG_MEV["tau"]
    try:
        PDG_MEV["tau"] = (Decimal("999.0"), Decimal("1.0"))
        doctored = conditional_tau_interval()
    except CertificateError:
        PDG_MEV["tau"] = saved
        return {
            "expected_failure": True,
            "failed": True,
            "code": "ROOT_ORDER",
            "meaning": "the doctored tau reference changes only the comparison fields",
        }
    finally:
        PDG_MEV["tau"] = saved
    same_enclosure = (
        doctored["tau_enclosure_mev_outward"] == original["tau_enclosure_mev_outward"]
        and doctored["tau_central_mev"] == original["tau_central_mev"]
    )
    try:
        require(
            not same_enclosure,
            "TAU_NOT_CONSUMED",
            "the tau enclosure is independent of the tau reference",
        )
    except CertificateError:
        return {
            "expected_failure": True,
            "failed": True,
            "code": "TAU_NOT_CONSUMED",
            "meaning": "the derived enclosure does not read the tau mass; only the comparison fields do",
        }
    return {"expected_failure": True, "failed": False}


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------


@high_precision
def build_payload() -> dict[str, Any]:
    comparison = balance_comparison()
    conditional = conditional_tau_interval()
    controls = {
        "premise_mutation": control_premise_mutation(),
        "tau_not_consumed": control_tau_not_consumed_in_lane_two(),
    }
    for name, verdict in controls.items():
        require(
            verdict["expected_failure"] is True and verdict["failed"] is True,
            "CONTROL_NOT_FAILED",
            f"control {name} did not record its required failure",
        )
    return {
        "schema": SCHEMA,
        "issue_context": ISSUE_CONTEXT,
        "claim_boundary": (
            "Two conditional/compare-only lanes on the positive-chamber "
            "circulant identity: the measured Q enclosure against the exact "
            "balance value 2/3, and the tau interval determined by "
            "(m_e, m_mu) under the balanced-circulant and mass-ordering "
            "premises. The premise class carries measured-target ancestry, "
            "declared above; neither lane is a prospective prediction, no "
            "source-only theorem is claimed, and the source-only charged "
            "no-go is unchanged."
        ),
        "identity": (
            "Q = 1/3 + (2/3)(rho/a)^2 on the positive chamber; Q = 2/3 "
            "exactly at rho/a = 1/sqrt2 "
            "(extra/koide_identity_from_positive_c3_face_circulants.tex)"
        ),
        "arithmetic": {
            "decimal_digits": DIGITS,
            "rounding": "outward at the displayed places",
            "enclosures": "one-sigma corner enumeration, conservative",
        },
        "measured_imports": {
            "source": PDG_SOURCE,
            "values_mev": {
                name: [str(value), str(sigma)]
                for name, (value, sigma) in PDG_MEV.items()
            },
        },
        "balance_comparison": comparison,
        "conditional_tau": conditional,
        "controls": controls,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=OUT_PATH)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    payload = build_payload()
    if args.verify:
        stored = json.loads(args.out.read_text(encoding="utf-8"))
        require(stored == payload, "DRIFT", "stored certificate does not match a rebuild")
        print(json.dumps({"status": "PASS"}))
        return 0
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "status": "WROTE",
                "tau_enclosure_mev": payload["conditional_tau"]["tau_enclosure_mev_outward"],
                "distance_sigma": payload["conditional_tau"]["distance_sigma"],
                "Q_distance": payload["balance_comparison"]["distance_in_enclosure_half_widths"],
            }
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
