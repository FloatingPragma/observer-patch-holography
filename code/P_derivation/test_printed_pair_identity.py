#!/usr/bin/env python3
"""CI check for closure-ledger row CL-6: the printed-pair identity.

Every runtime artifact that prints a closure pair must satisfy

    alpha_root = (P - phi) / sqrt(pi),    alpha_root = 1 / alpha_inv,

because ``P`` is defined by the outer equation ``P = phi + alpha*sqrt(pi)``.
An unconverged solver run prints a probe-side ``P`` and a map-side
``alpha_inv`` that disagree at the truncation depth of the run; row CL-6
recorded exactly that defect. The artifacts are regenerated at solver
precision 100 with enough outer bisection iterations that the identity holds
far below the stated tolerance.

Stated tolerance: relative defect <= 1e-30 (at least 30 significant digits).
"""

from __future__ import annotations

from decimal import Decimal, localcontext
import json
from pathlib import Path

from paper_math import decimal_pi


RUNTIME = Path(__file__).resolve().parent / "runtime"
TRUNK = RUNTIME / "p_closure_trunk_current.json"
FULL_REPORT = RUNTIME / "full_p_alpha_report_current.json"

RELATIVE_TOLERANCE = Decimal("1e-30")
WORK_PRECISION = 160


def _constants() -> tuple[Decimal, Decimal]:
    with localcontext() as ctx:
        ctx.prec = WORK_PRECISION
        pi = +decimal_pi(WORK_PRECISION)
        sqrt_pi = pi.sqrt()
        phi = (Decimal(1) + Decimal(5).sqrt()) / Decimal(2)
        return +phi, +sqrt_pi


def _relative_identity_defect(p: str, alpha_inv: str) -> Decimal:
    phi, sqrt_pi = _constants()
    with localcontext() as ctx:
        ctx.prec = WORK_PRECISION
        alpha_root = Decimal(1) / Decimal(alpha_inv)
        rhs = (Decimal(p) - phi) / sqrt_pi
        return +(abs(alpha_root - rhs) / alpha_root)


def test_decimal_pi_matches_independent_mpmath_pi() -> None:
    """The chain's Chudnovsky pi agrees with an independent mpmath pi."""
    from mpmath.ctx_mp import MPContext

    mp = MPContext()
    mp.dps = 140
    independent = mp.nstr(mp.pi, 120, strip_zeros=False)
    ours = format(decimal_pi(150), "f")[:121]
    assert ours.startswith(independent[:118])


def test_trunk_printed_pair_identity() -> None:
    payload = json.loads(TRUNK.read_text(encoding="utf-8"))
    fixed_point = payload["fixed_point_candidate"]
    defect = _relative_identity_defect(fixed_point["P"], fixed_point["alpha_inv"])
    assert defect <= RELATIVE_TOLERANCE

    # The probe-side alpha printed next to the pair must agree with the
    # map-side root at the same tolerance.
    with localcontext() as ctx:
        ctx.prec = WORK_PRECISION
        alpha = Decimal(fixed_point["alpha"])
        alpha_root = Decimal(1) / Decimal(fixed_point["alpha_inv"])
        assert abs(alpha - alpha_root) / alpha_root <= RELATIVE_TOLERANCE

    # The artifact's own phi and sqrt(pi) strings match independent values.
    phi, sqrt_pi = _constants()
    closed_form = payload["closed_form_candidate"]
    with localcontext() as ctx:
        ctx.prec = WORK_PRECISION
        assert abs(Decimal(closed_form["phi"]) - phi) / phi <= Decimal("1e-30")
        assert abs(Decimal(closed_form["sqrt_pi"]) - sqrt_pi) / sqrt_pi <= Decimal("1e-30")


def test_full_report_printed_pair_identity() -> None:
    payload = json.loads(FULL_REPORT.read_text(encoding="utf-8"))
    defect = _relative_identity_defect(payload["p"], payload["alpha_inv"])
    assert defect <= RELATIVE_TOLERANCE

    with localcontext() as ctx:
        ctx.prec = WORK_PRECISION
        residual = abs(Decimal(payload["alpha_fixed_point_residual"]))
        alpha = Decimal(payload["alpha"])
        assert residual / alpha <= RELATIVE_TOLERANCE
