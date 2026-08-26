#!/usr/bin/env python3
"""Build the local, conditional capacity--alpha tangent receipt.

This module does one deliberately narrow calculation.  It differentiates the
committed ``PaperMathContext.solve_alpha_u_from_p`` root near the CODATA-fed
comparison coordinate and propagates three already-public alpha-constancy
measurements through that *local* tangent.

The calculation is not a derivation of cosmic evolution.  Its physical use
requires all three explicit premises in :data:`PREMISES`:

* B1 identifies the electroweak bridge coordinate with the capacity in the
  conditional w-law;
* B2 makes that bridge a differentiable epoch-by-epoch law near the endpoint;
* B3 says physical evolution follows the same committed solver root and
  co-variation convention whose tangent is evaluated here.

All observations used below predate this calculation.  The resulting receipt
is therefore a retrospective conditional diagnostic, never OPH evidence or a
prediction verdict.  In particular, this module refuses to extrapolate the
local tangent across the large integrated capacity changes suggested by CPL
best fits.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from decimal import Decimal, localcontext
from pathlib import Path


HERE = Path(__file__).resolve().parent
CODE_ROOT = HERE.parents[1]
P_DERIVATION = CODE_ROOT / "P_derivation"
PAPER_MATH = P_DERIVATION / "paper_math.py"
DEFAULT_OUTPUT = HERE / "runtime" / "capacity_alpha_tangent_retrospective.json"

SCHEMA = "oph.capacity_alpha_local_conditional_retrospective.v1"
DECIMAL_PRECISION = 24
WORK_PRECISION = 60
ALPHA_INVERSE_CODATA_2022 = Decimal("137.035999177")
H0_KM_S_MPC = Decimal("67.4")
MPC_IN_KM = Decimal("3.0856775814913673e19")
JULIAN_YEAR_SECONDS = Decimal("31557600")
GAUSSIAN_95_TWO_SIDED = Decimal("1.96")

# This is a reporting guard, not a theorem-level truncation-error bound.  It
# keeps every emitted finite-difference diagnostic below one percent in
# |Delta ln alpha|.  The largest measurement propagated here is about 0.52%.
# Any larger use needs an exact branch solve and an error analysis.
LOCAL_REPORT_MAX_ABS_DELTA_LN_ALPHA = Decimal("0.01")

H_SWEEP = (
    Decimal("1e-3"),
    Decimal("5e-4"),
    Decimal("1e-4"),
    Decimal("5e-5"),
)
REFERENCE_H = Decimal("1e-4")
DERIVATIVE_RELATIVE_SPREAD_LIMIT = Decimal("2e-6")


PREMISES = {
    "B1_same_capacity": {
        "statement": (
            "The conditional electroweak bridge coordinate N_EW is the same "
            "physical capacity N that appears in the fixed-capacity w-law."
        ),
        "status": "undischarged_physical_identification",
    },
    "B2_epochwise_bridge": {
        "statement": (
            "The bridge N = pi exp(6 pi/(P alpha_U(P))) holds differentiably "
            "at each epoch in a neighborhood of the comparison endpoint."
        ),
        "status": "undischarged_epochwise_law",
    },
    "B3_physical_solver_tangent": {
        "statement": (
            "The dimensionless alpha read by local optical clocks is the same "
            "homogeneous cosmological alpha, its present proper-time drift "
            "equals the background cosmic-time drift, and physical epoch "
            "evolution follows the same committed PaperMathContext root branch "
            "and co-variation convention used when differentiating alpha_U(P), "
            "with the remaining solver inputs and conventions held as declared."
        ),
        "status": "undischarged_dynamical_tangent_attachment",
    },
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _decimal_text(value: Decimal) -> str:
    """Serialize a Decimal without converting through binary floating point."""

    return str(value)


def conservative_gaussian_abs_envelope(
    central: Decimal,
    sigma: Decimal,
    *,
    z_value: Decimal = GAUSSIAN_95_TWO_SIDED,
) -> Decimal:
    """Return ``|central| + z sigma``, a conservative symmetric abs envelope."""

    if sigma < 0 or z_value <= 0:
        raise ValueError("sigma must be nonnegative and z_value must be positive")
    return abs(central) + z_value * sigma


def fractional_to_abs_log_envelope(fractional_envelope: Decimal) -> Decimal:
    """Conservatively convert ``|Delta alpha/alpha| <= e`` to log space.

    For ``0 <= e < 1``, the larger of ``ln(1+e)`` and ``-ln(1-e)`` is the
    latter.  Using it avoids silently identifying a finite fractional change
    with a logarithmic change.
    """

    if not Decimal(0) <= fractional_envelope < Decimal(1):
        raise ValueError("fractional envelope must lie in [0, 1)")
    with localcontext() as ctx:
        ctx.prec = WORK_PRECISION
        return +(-(Decimal(1) - fractional_envelope).ln())


def jacobian_from_probe(
    *,
    p: Decimal,
    alpha: Decimal,
    alpha_u: Decimal,
    d_alpha_u_d_p: Decimal,
    pi: Decimal,
    sqrt_pi: Decimal,
) -> dict[str, Decimal]:
    """Evaluate the exact chain-rule expressions at one local probe."""

    if min(p, alpha, alpha_u, pi, sqrt_pi) <= 0:
        raise ValueError("positive probe inputs required")
    with localcontext() as ctx:
        ctx.prec = WORK_PRECISION
        dln_alpha_u_dln_p = d_alpha_u_d_p * p / alpha_u
        dln_p_alpha_u_dln_p = Decimal(1) + dln_alpha_u_dln_p
        log_n_over_pi = Decimal(6) * pi / (p * alpha_u)
        dln_n_dln_p = -log_n_over_pi * dln_p_alpha_u_dln_p
        dln_p_dln_alpha = sqrt_pi * alpha / p
        dln_n_dln_alpha = dln_n_dln_p * dln_p_dln_alpha
        if dln_n_dln_alpha == 0:
            raise ValueError("capacity-alpha tangent is singular")
        return {
            "dln_alpha_u_dln_p": +dln_alpha_u_dln_p,
            "dln_p_alpha_u_dln_p": +dln_p_alpha_u_dln_p,
            "log_n_over_pi": +log_n_over_pi,
            "dln_n_dln_p": +dln_n_dln_p,
            "dln_p_dln_alpha": +dln_p_dln_alpha,
            "dln_n_dln_alpha": +dln_n_dln_alpha,
            "dln_alpha_dln_n": +(Decimal(1) / dln_n_dln_alpha),
        }


def linearized_delta_ln_alpha_from_capacity(
    delta_ln_n: Decimal,
    dln_alpha_dln_n: Decimal,
    *,
    max_abs_delta_ln_alpha: Decimal = LOCAL_REPORT_MAX_ABS_DELTA_LN_ALPHA,
) -> Decimal:
    """Return a guarded first-order response, rejecting large extrapolation."""

    response = delta_ln_n * dln_alpha_dln_n
    if abs(response) > max_abs_delta_ln_alpha:
        raise ValueError(
            "local capacity-alpha tangent extrapolation forbidden: "
            f"|linearized Delta ln alpha|={abs(response)} exceeds "
            f"the reporting guard {max_abs_delta_ln_alpha}; solve the exact "
            "epoch branch with a truncation-error receipt instead"
        )
    return response


def _solver_h_sweep() -> tuple[object, Decimal, Decimal, list[dict[str, str]], Decimal]:
    """Run a deterministic centered h-sweep on the committed solver."""

    sys.path.insert(0, str(P_DERIVATION))
    try:
        from paper_math import PaperMathContext
    finally:
        sys.path.pop(0)

    context = PaperMathContext(precision=DECIMAL_PRECISION)
    with localcontext() as ctx:
        ctx.prec = WORK_PRECISION
        alpha = Decimal(1) / ALPHA_INVERSE_CODATA_2022
        p_c = context.outer_p_from_alpha(alpha)
    alpha_u_center, _report, _mu_u = context.solve_alpha_u_from_p(p_c)

    rows: list[dict[str, str]] = []
    derivative_by_h: dict[Decimal, Decimal] = {}
    for h in H_SWEEP:
        alpha_u_plus, _plus_report, _plus_mu_u = context.solve_alpha_u_from_p(p_c + h)
        alpha_u_minus, _minus_report, _minus_mu_u = context.solve_alpha_u_from_p(p_c - h)
        derivative = (alpha_u_plus - alpha_u_minus) / (Decimal(2) * h)
        derivative_by_h[h] = derivative
        rows.append(
            {
                "h": _decimal_text(h),
                "alpha_u_minus": _decimal_text(alpha_u_minus),
                "alpha_u_plus": _decimal_text(alpha_u_plus),
                "centered_d_alpha_u_d_p": _decimal_text(derivative),
            }
        )

    reference = derivative_by_h[REFERENCE_H]
    relative_spread = max(
        abs((derivative - reference) / reference)
        for derivative in derivative_by_h.values()
    )
    if relative_spread > DERIVATIVE_RELATIVE_SPREAD_LIMIT:
        raise RuntimeError(
            "committed solver h-sweep failed its deterministic stability gate: "
            f"{relative_spread} > {DERIVATIVE_RELATIVE_SPREAD_LIMIT}"
        )
    return context, p_c, alpha_u_center, rows, relative_spread


def _integrated_measurement(
    *,
    name: str,
    redshift: Decimal,
    central_fraction: Decimal,
    sigma_fraction: Decimal,
    dln_n_dln_alpha: Decimal,
    source: str,
) -> dict[str, object]:
    envelope = conservative_gaussian_abs_envelope(
        central_fraction, sigma_fraction
    )
    log_envelope = fractional_to_abs_log_envelope(envelope)
    if log_envelope > LOCAL_REPORT_MAX_ABS_DELTA_LN_ALPHA:
        raise ValueError(f"{name} lies outside the local reporting guard")
    with localcontext() as ctx:
        ctx.prec = WORK_PRECISION
        linearized_delta_ln_n = abs(dln_n_dln_alpha) * log_envelope
        signed_log_average_bound = linearized_delta_ln_n / (
            Decimal(3) * (Decimal(1) + redshift).ln()
        )
    return {
        "name": name,
        "source": source,
        "data_status": "retrospective_public_before_calculation",
        "central_fractional_change": _decimal_text(central_fraction),
        "one_sigma_fractional_change": _decimal_text(sigma_fraction),
        "conservative_95_percent_abs_fractional_envelope": _decimal_text(envelope),
        "conservative_abs_delta_ln_alpha_envelope": _decimal_text(log_envelope),
        "linearized_abs_delta_ln_n_envelope": _decimal_text(linearized_delta_ln_n),
        "linearized_abs_signed_log_average_1_plus_w_envelope": _decimal_text(
            signed_log_average_bound
        ),
        "average_scope": (
            "Bound on the absolute signed log-average. It is a bound on the "
            "average absolute value only under a no-cancellation condition, "
            "such as the monotone w>=-1 branch."
        ),
        "linearization_status": (
            "first_order_local_diagnostic_without_a_global_truncation_error_bound"
        ),
    }


def observational_diagnostics(
    dln_n_dln_alpha: Decimal,
) -> dict[str, object]:
    """Propagate the corrected public-data envelopes through the tangent."""

    with localcontext() as ctx:
        ctx.prec = WORK_PRECISION
        h0_per_year = H0_KM_S_MPC / MPC_IN_KM * JULIAN_YEAR_SECONDS

        filzinger_central = Decimal("1.8e-19")
        filzinger_sigma = Decimal("2.5e-19")
        filzinger_envelope = conservative_gaussian_abs_envelope(
            filzinger_central, filzinger_sigma
        )
        one_plus_w0 = (
            abs(dln_n_dln_alpha)
            * filzinger_envelope
            / (Decimal(3) * h0_per_year)
        )

        espresso_sigma = (
            Decimal("1.3e-6") ** 2 + Decimal("0.4e-6") ** 2
        ).sqrt()

    return {
        "H0_conversion": {
            "H0_km_s_Mpc": _decimal_text(H0_KM_S_MPC),
            "Mpc_in_km": _decimal_text(MPC_IN_KM),
            "Julian_year_seconds": _decimal_text(JULIAN_YEAR_SECONDS),
            "H0_per_Julian_year": _decimal_text(h0_per_year),
        },
        "Filzinger_2023_clock": {
            "source": "Filzinger et al., Phys. Rev. Lett. 130, 253001 (2023)",
            "data_status": "retrospective_public_before_calculation",
            "central_dln_alpha_dt_per_year": _decimal_text(filzinger_central),
            "one_sigma_per_year": _decimal_text(filzinger_sigma),
            "conservative_95_percent_abs_rate_envelope_per_year": _decimal_text(
                filzinger_envelope
            ),
            "conditional_abs_1_plus_w0_envelope": _decimal_text(one_plus_w0),
            "mapping_status": "local_differential_and_conditional_on_B1_B2_B3",
        },
        "ESPRESSO_HE0515_4414": _integrated_measurement(
            name="ESPRESSO HE 0515-4414",
            redshift=Decimal("1.15"),
            central_fraction=Decimal("1.3e-6"),
            sigma_fraction=espresso_sigma,
            dln_n_dln_alpha=dln_n_dln_alpha,
            source=(
                "Murphy et al., Astron. Astrophys. 658, A123 (2022): "
                "1.3 +/- 1.3(stat) +/- 0.4(sys) ppm"
            ),
        ),
        "Hart_Chluba_Planck_2018": _integrated_measurement(
            name="Planck 2018 recombination alpha, Hart-Chluba",
            redshift=Decimal("1100"),
            central_fraction=Decimal("0.0005"),
            sigma_fraction=Decimal("0.0024"),
            dln_n_dln_alpha=dln_n_dln_alpha,
            source=(
                "Hart and Chluba, Mon. Not. R. Astron. Soc. 493, 3255 "
                "(2020): alpha/alpha0 = 1.0005 +/- 0.0024 (CMB only)"
            ),
        ),
    }


def build_receipt() -> dict[str, object]:
    context, p_c, alpha_u, h_rows, relative_spread = _solver_h_sweep()
    reference_row = next(row for row in h_rows if row["h"] == str(REFERENCE_H))
    d_alpha_u_d_p = Decimal(reference_row["centered_d_alpha_u_d_p"])
    with localcontext() as ctx:
        ctx.prec = WORK_PRECISION
        alpha = Decimal(1) / ALPHA_INVERSE_CODATA_2022
        jacobian = jacobian_from_probe(
            p=p_c,
            alpha=alpha,
            alpha_u=alpha_u,
            d_alpha_u_d_p=d_alpha_u_d_p,
            pi=context.pi,
            sqrt_pi=context.sqrt_pi,
        )
        max_abs_delta_ln_n = (
            LOCAL_REPORT_MAX_ABS_DELTA_LN_ALPHA
            / abs(jacobian["dln_alpha_dln_n"])
        )

    return {
        "schema": SCHEMA,
        "classification": {
            "epistemic_status": (
                "retrospective conditional diagnostic; no OPH evidence, "
                "validation, confirmation, falsification, or prediction score"
            ),
            "local_tangent_only": True,
            "measured_target_ancestry": (
                "P_C is built from the CODATA 2022 comparison value of alpha; "
                "it is not the independently source-forward P coordinate"
            ),
            "physical_attachment": "requires every one of B1, B2, and B3",
        },
        "premises": PREMISES,
        "solver": {
            "implementation": "PaperMathContext.solve_alpha_u_from_p",
            "paper_math_path": str(PAPER_MATH.relative_to(HERE.parents[2])),
            "paper_math_sha256": _sha256(PAPER_MATH),
            "receipt_script_sha256": _sha256(Path(__file__).resolve()),
            "decimal_precision": DECIMAL_PRECISION,
            "root_selection_scope": (
                "the first sign-changing bracket on the committed [0.02,0.08] "
                "scan; no global root uniqueness is claimed by this receipt"
            ),
            "reference_h": _decimal_text(REFERENCE_H),
            "relative_derivative_spread": _decimal_text(relative_spread),
            "relative_derivative_spread_limit": _decimal_text(
                DERIVATIVE_RELATIVE_SPREAD_LIMIT
            ),
            "h_sweep_passed": True,
            "h_sweep": h_rows,
        },
        "comparison_coordinate": {
            "alpha_inverse_CODATA_2022": _decimal_text(
                ALPHA_INVERSE_CODATA_2022
            ),
            "alpha": _decimal_text(alpha),
            "P_C": _decimal_text(p_c),
            "alpha_U_at_P_C": _decimal_text(alpha_u),
            "d_alpha_U_d_P": _decimal_text(d_alpha_u_d_p),
        },
        "local_jacobian": {
            key: _decimal_text(value) for key, value in jacobian.items()
        },
        "observational_diagnostics": observational_diagnostics(
            jacobian["dln_n_dln_alpha"]
        ),
        "extrapolation_guard": {
            "max_abs_delta_ln_alpha_for_emitted_local_report": _decimal_text(
                LOCAL_REPORT_MAX_ABS_DELTA_LN_ALPHA
            ),
            "corresponding_linearized_max_abs_delta_ln_n": _decimal_text(
                max_abs_delta_ln_n
            ),
            "integrated_DESI_CPL_extrapolation": "forbidden",
            "reason": (
                "The committed calculation certifies a local tangent only. "
                "Large CPL-integrated capacity changes require an exact "
                "epoch-dependent branch solve, branch/attachment premises, "
                "and a truncation-error receipt."
            ),
        },
        "nonclaims": [
            "The receipt does not derive B1, B2, or B3.",
            "The receipt does not identify P_C with the source-forward P value.",
            "The receipt does not close or retire the monotone-capacity branch.",
            "The receipt does not turn retrospective alpha data into OPH evidence.",
            "The receipt does not propagate DESI CPL fits through a large local-linear extrapolation.",
        ],
    }


def write_receipt(output: Path, receipt: dict[str, object]) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    receipt = build_receipt()
    write_receipt(args.output, receipt)
    print(json.dumps(receipt, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
