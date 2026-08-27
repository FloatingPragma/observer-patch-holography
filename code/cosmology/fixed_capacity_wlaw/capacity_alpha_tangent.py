#!/usr/bin/env python3
"""Build the local, conditional capacity--alpha comparison receipt.

The tangent is supplied by ``capacity_alpha_interval_certificate.py``.  That
certificate proves C1 regularity of the selected declared root on a whole
local box, records sign-definite implicit denominators, encloses
``d log N/d log alpha``, and computes its reciprocal only after excluding
zero.  The older centered finite-difference sweep remains available as a
non-proof regression diagnostic, but it no longer supports the receipt.

The calculation is not a derivation of cosmic evolution.  Its physical use
requires all three explicit premises in :data:`PREMISES`:

* B1 identifies the electroweak bridge coordinate with the capacity in the
  conditional w-law;
* B2 makes that bridge a differentiable epoch-by-epoch law near the endpoint;
* B3 says physical evolution follows the same committed solver root and
  co-variation convention whose tangent is evaluated here.

All observations used below predate this calculation.  The resulting receipt
is therefore a retrospective conditional diagnostic, never OPH evidence or a
prediction verdict.  This module emits a finite-change bound only when the
complete alpha interval lies inside the certified mean-value domain.
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
INTERVAL_PRODUCER = HERE / "capacity_alpha_interval_certificate.py"
INTERVAL_RECEIPT = HERE / "runtime" / "capacity_alpha_interval_certificate.json"

SCHEMA = "oph.capacity_alpha_local_conditional_retrospective.v2"
DECIMAL_PRECISION = 24
WORK_PRECISION = 60
ALPHA_INVERSE_CODATA_2022 = Decimal("137.035999177")
H0_KM_S_MPC = Decimal("67.4")
MPC_IN_KM = Decimal("3.0856775814913673e19")
JULIAN_YEAR_SECONDS = Decimal("31557600")
GAUSSIAN_95_TWO_SIDED = Decimal("1.96")

# Compatibility guard for the legacy point-linear helper below.  The receipt
# uses the narrower rigorous interval domain recorded by the interval
# certificate and never promotes this guard to a theorem.
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
            "The finite relation N = pi exp(6 pi/(P alpha_U(P))) is a physical "
            "epoch-by-epoch law on the domain used by the cosmological model."
        ),
        "status": "undischarged_physical_epochwise_law",
    },
    "B3_physical_solver_tangent": {
        "statement": (
            "Physical evolution selects this committed root branch and "
            "co-variation convention, while a local optical-clock alpha is the "
            "same homogeneous cosmological alpha with the required time map."
        ),
        "status": "undischarged_physical_branch_selection_and_readout",
        "mathematical_subresult": (
            "C1 regularity and a sign-definite tangent of the selected declared "
            "root branch are interval-certified; this does not select the "
            "branch physically."
        ),
    },
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _build_interval_certificate() -> dict[str, object]:
    """Build the rigorous local certificate through its canonical producer."""

    sys.path.insert(0, str(HERE))
    try:
        import capacity_alpha_interval_certificate as interval_certificate
    finally:
        sys.path.pop(0)
    return interval_certificate.build_certificate()


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
    abs_tangent_upper: Decimal,
    certified_log_radius: Decimal,
    source: str,
) -> dict[str, object]:
    envelope = conservative_gaussian_abs_envelope(
        central_fraction, sigma_fraction
    )
    log_envelope = fractional_to_abs_log_envelope(envelope)
    common: dict[str, object] = {
        "name": name,
        "source": source,
        "data_status": "retrospective_public_before_calculation",
        "central_fractional_change": _decimal_text(central_fraction),
        "one_sigma_fractional_change": _decimal_text(sigma_fraction),
        "conservative_95_percent_abs_fractional_envelope": _decimal_text(envelope),
        "conservative_abs_delta_ln_alpha_envelope": _decimal_text(log_envelope),
        "guaranteed_symmetric_log_alpha_inner_radius": _decimal_text(
            certified_log_radius
        ),
    }
    if log_envelope > certified_log_radius:
        return {
            **common,
            "interval_mapping_status": "outside_certified_domain_no_bound_emitted",
            "reason": (
                "The public-data envelope is wider than the rigorous local "
                "mean-value domain.  A partitioned or global branch certificate "
                "is required before propagating it."
            ),
        }
    with localcontext() as ctx:
        ctx.prec = WORK_PRECISION
        delta_ln_n_bound = abs_tangent_upper * log_envelope
        signed_log_average_bound = delta_ln_n_bound / (
            Decimal(3) * (Decimal(1) + redshift).ln()
        )
    return {
        **common,
        "interval_mean_value_abs_delta_ln_n_envelope": _decimal_text(
            delta_ln_n_bound
        ),
        "interval_mean_value_abs_signed_log_average_1_plus_w_envelope": _decimal_text(
            signed_log_average_bound
        ),
        "average_scope": (
            "Bound on the absolute signed log-average. It is a bound on the "
            "average absolute value only under a no-cancellation condition, "
            "such as the monotone w>=-1 branch."
        ),
        "interval_mapping_status": "inside_certified_domain_mean_value_bound",
    }


def observational_diagnostics(
    abs_tangent_upper: Decimal,
    certified_log_radius: Decimal,
) -> dict[str, object]:
    """Propagate only envelopes covered by the rigorous local domain."""

    with localcontext() as ctx:
        ctx.prec = WORK_PRECISION
        h0_per_year = H0_KM_S_MPC / MPC_IN_KM * JULIAN_YEAR_SECONDS

        filzinger_central = Decimal("1.8e-19")
        filzinger_sigma = Decimal("2.5e-19")
        filzinger_envelope = conservative_gaussian_abs_envelope(
            filzinger_central, filzinger_sigma
        )
        one_plus_w0 = (
            abs_tangent_upper
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
            abs_tangent_upper=abs_tangent_upper,
            certified_log_radius=certified_log_radius,
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
            abs_tangent_upper=abs_tangent_upper,
            certified_log_radius=certified_log_radius,
            source=(
                "Hart and Chluba, Mon. Not. R. Astron. Soc. 493, 3255 "
                "(2020): alpha/alpha0 = 1.0005 +/- 0.0024 (CMB only)"
            ),
        ),
    }


def build_receipt() -> dict[str, object]:
    interval_certificate = _build_interval_certificate()
    try:
        committed_interval = json.loads(INTERVAL_RECEIPT.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError("canonical interval certificate is missing or malformed") from exc
    if committed_interval != interval_certificate:
        raise RuntimeError(
            "canonical interval certificate is stale; regenerate and verify it first"
        )

    branch = interval_certificate["branch_certificate"]
    derivatives = branch["local_derivatives"]
    mean_value = branch["mean_value_certificate"]
    certified_log_radius = Decimal(
        branch["certified_domain"][
            "guaranteed_symmetric_log_alpha_inner_radius"
        ]
    )
    abs_tangent_upper = Decimal(mean_value["abs_slope_upper"])

    return {
        "schema": SCHEMA,
        "classification": {
            "epistemic_status": (
                "retrospective conditional diagnostic; no OPH evidence, "
                "validation, confirmation, falsification, or prediction score"
            ),
            "rigorous_local_interval_tangent": True,
            "mathematical_branch_differentiability": "attained_on_certified_domain",
            "physical_epoch_evolution": "undischarged",
            "measured_target_ancestry": (
                "P_C is built from the CODATA 2022 comparison value of alpha; "
                "it is not the independently source-forward P coordinate"
            ),
            "physical_attachment": "requires every one of B1, B2, and B3",
        },
        "premises": PREMISES,
        "solver": {
            "implementation": (
                "outward-rounded interval arithmetic with sign-change root "
                "brackets and implicit-function derivatives"
            ),
            "receipt_script_sha256": _sha256(Path(__file__).resolve()),
            "paper_math_path": str(PAPER_MATH.relative_to(HERE.parents[2])),
            "paper_math_sha256": _sha256(PAPER_MATH),
            "interval_certificate_path": str(
                INTERVAL_RECEIPT.relative_to(HERE.parents[2])
            ),
            "interval_certificate_sha256": _sha256(INTERVAL_RECEIPT),
            "interval_producer_path": str(
                INTERVAL_PRODUCER.relative_to(HERE.parents[2])
            ),
            "interval_producer_sha256": _sha256(INTERVAL_PRODUCER),
            "root_selection_scope": (
                branch["implicit_function_certificate"]["selection_scope"]
            ),
            "finite_difference_role": (
                "optional regression diagnostic only; no finite difference "
                "supports the emitted derivative enclosure"
            ),
        },
        "comparison_coordinate": branch["comparison_coordinate"],
        "certified_domain": branch["certified_domain"],
        "implicit_function_certificate": branch[
            "implicit_function_certificate"
        ],
        "local_jacobian_interval": derivatives,
        "mean_value_certificate": mean_value,
        "observational_diagnostics": observational_diagnostics(
            abs_tangent_upper, certified_log_radius
        ),
        "extrapolation_guard": {
            "guaranteed_symmetric_log_alpha_inner_radius": _decimal_text(
                certified_log_radius
            ),
            "integrated_DESI_CPL_extrapolation": "forbidden",
            "reason": (
                "Finite-change propagation is emitted only inside the rigorous "
                "mean-value domain.  Larger changes require a partitioned or "
                "global interval branch certificate plus the physical premises."
            ),
        },
        "nonclaims": [
            "The receipt does not derive B1, B2, or B3.",
            "The receipt does not identify P_C with the source-forward P value.",
            "The receipt does not close or retire the monotone-capacity branch.",
            "The receipt does not turn retrospective alpha data into OPH evidence.",
            "The receipt does not propagate DESI CPL fits outside the certified interval domain.",
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
