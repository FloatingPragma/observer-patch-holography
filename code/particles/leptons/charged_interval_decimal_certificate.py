#!/usr/bin/env python3
"""Independent decimal enclosures for the charged-lepton closure intervals.

The production builders use binary floating point for their diagnostic
centrals.  This module recomputes the interval endpoints from the serialized
decimal inputs at 100-digit precision, brackets pi, adds an explicit numerical
safety pad, and rounds the reported bounds outward.  Its certificate concerns
the arithmetic implication from the recorded inputs.  It does not change their
epistemic status: the ratios are target anchored, the electromagnetic endpoint
is empirical, and the higher-order remainder is calibrated at the measured
charged-lepton triple.
"""

from __future__ import annotations

from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
from typing import Any, Mapping

PRECISION = 100
PI_LO = Decimal(
    "3.1415926535897932384626433832795028841971693993751058209749445923078164"
)
PI_HI = Decimal(
    "3.1415926535897932384626433832795028841971693993751058209749445923078165"
)
SAFETY_PAD = Decimal("1e-60")
ROUNDING_ERROR_BOUND = Decimal("1e-70")
OUTWARD_QUANTUM = Decimal("1e-15")
MZ_GEV = Decimal("91.1876")
KERNEL_TRUNCATION_PACKET = Decimal("0.0005")
ALPHA_INV_0_UNCERTAINTY = Decimal("2.1e-8")
MASS_ORDER = ("electron", "muon", "tau")


def _d(value: Any) -> Decimal:
    return Decimal(str(value))


def _shape_inputs(
    readout: Mapping[str, Any],
) -> tuple[tuple[Decimal, Decimal], tuple[Decimal, Decimal, Decimal]]:
    centered = tuple(_d(value) for value in readout["centered_log_shape_exact"])
    witness = tuple(_d(value) for value in readout["predicted_singular_values_abs"])
    if len(centered) != 3 or len(witness) != 3:
        raise ValueError("charged readout must contain three centered logs and masses")
    ratios = (
        (centered[1] - centered[0]).exp(),
        (centered[2] - centered[0]).exp(),
    )
    return ratios, witness


def _lepton_packet(
    pi: Decimal,
    witness_e: Decimal,
    ratios: tuple[Decimal, Decimal],
) -> Decimal:
    masses = (witness_e, ratios[0] * witness_e, ratios[1] * witness_e)
    return sum(
        (Decimal(1) / (Decimal(3) * pi))
        * ((MZ_GEV * MZ_GEV / (mass * mass)).ln() - Decimal(5) / Decimal(3))
        for mass in masses
    )


def _invert_packet(
    pi: Decimal,
    packet: Decimal,
    ratios: tuple[Decimal, Decimal],
) -> Decimal:
    log_ratio_sum = ratios[0].ln() + ratios[1].ln()
    return (
        MZ_GEV.ln()
        - (
            Decimal(3) * pi * packet
            + Decimal(2) * log_ratio_sum
            + Decimal(5)
        )
        / Decimal(6)
    ).exp()


def _outward(value: Decimal, *, lower: bool) -> str:
    rounding = ROUND_FLOOR if lower else ROUND_CEILING
    return format(value.quantize(OUTWARD_QUANTUM, rounding=rounding), "f")


def _certificate_payload(
    kappa_lo: Decimal,
    kappa_hi: Decimal,
    mass_bounds: Mapping[str, tuple[Decimal, Decimal]],
) -> dict[str, Any]:
    if SAFETY_PAD <= ROUNDING_ERROR_BOUND:
        raise ArithmeticError("numerical safety pad must exceed the rounding bound")
    return {
        "method": (
            "100-digit Decimal recomputation at both endpoints of an explicit "
            "pi bracket, padded by 1e-60 and rounded outward"
        ),
        "input_semantics": (
            "serialized decimal inputs are treated as exact values; the "
            "certificate encloses the arithmetic consequence of those inputs"
        ),
        "epistemic_scope": (
            "target-anchored empirical-closure diagnostic, not a source-only "
            "or prospective charged-lepton prediction"
        ),
        "precision_decimal_digits": PRECISION,
        "pi_interval": [str(PI_LO), str(PI_HI)],
        "safety_pad": str(SAFETY_PAD),
        "rounding_error_bound": str(ROUNDING_ERROR_BOUND),
        "rounding_error_bound_basis": (
            "the fixed expression has fewer than 200 operations on magnitudes "
            "below 1e4; correctly rounded 100-digit Decimal elementary "
            "operations contribute less than 1e-90 in total, so 1e-70 is a "
            "conservative absolute bound"
        ),
        "safety_margin_over_rounding_error": str(
            SAFETY_PAD / ROUNDING_ERROR_BOUND
        ),
        "outward_quantum": str(OUTWARD_QUANTUM),
        "monotonicity": (
            "the one-loop packet is strictly decreasing in the common log-mass "
            "shift; the declared endpoint extremes therefore give the interval "
            "endpoints"
        ),
        "kappa_interval": [
            _outward(kappa_lo, lower=True),
            _outward(kappa_hi, lower=False),
        ],
        "mass_intervals_gev": {
            particle: [
                _outward(bounds[0], lower=True),
                _outward(bounds[1], lower=False),
            ]
            for particle, bounds in mass_bounds.items()
        },
    }


def coherent_certificate(
    readout: Mapping[str, Any],
    endpoint: Mapping[str, Any],
    bridge: Mapping[str, Any],
) -> dict[str, Any]:
    """Enclose the payload-coherent charged interval."""

    with localcontext() as context:
        context.prec = PRECISION
        ratios, witness = _shape_inputs(readout)
        reference = bridge["reference_decomposition_compare_only"]
        alpha_inv = _d(reference["alpha_inv_0"])
        delta_top = _d(reference["Delta_top"])
        delta_lep_3loop = _d(reference["Delta_lep"])
        a_lep = _d(endpoint["transport_split"]["lepton_transport_delta_inv_alpha"])
        alpha_width_packet = (
            a_lep * ALPHA_INV_0_UNCERTAINTY / alpha_inv
        )

        def evaluate(
            pi: Decimal,
        ) -> tuple[Decimal, Decimal, dict[str, Decimal]]:
            packet_witness = _lepton_packet(pi, witness[0], ratios)
            remainder = delta_lep_3loop - packet_witness / alpha_inv
            if remainder <= 0:
                raise ValueError("higher-order charged-lepton remainder must be positive")

            def solve(
                higher_order: Decimal,
                kernel_slack: Decimal,
                alpha_slack: Decimal,
            ) -> Decimal:
                required = a_lep / alpha_inv - delta_top - higher_order
                packet = required * alpha_inv + kernel_slack + alpha_slack
                mass_e = _invert_packet(pi, packet, ratios)
                return (mass_e / witness[0]).ln()

            kappa_lo = solve(
                Decimal(0),
                KERNEL_TRUNCATION_PACKET,
                alpha_width_packet,
            )
            kappa_hi = solve(
                Decimal(2) * remainder,
                -KERNEL_TRUNCATION_PACKET,
                -alpha_width_packet,
            )
            factors = (Decimal(1), ratios[0], ratios[1])
            mass_bounds = {
                particle: witness[0] * factor * kappa.exp()
                for particle, factor, kappa in zip(
                    MASS_ORDER,
                    factors,
                    (kappa_lo, kappa_lo, kappa_lo),
                    strict=True,
                )
            }
            mass_upper = {
                particle: witness[0] * factor * kappa_hi.exp()
                for particle, factor in zip(MASS_ORDER, factors, strict=True)
            }
            return kappa_lo, kappa_hi, {
                particle: mass_bounds[particle]
                for particle in MASS_ORDER
            } | {
                f"{particle}__upper": mass_upper[particle]
                for particle in MASS_ORDER
            }

        # Keep lower and upper mass values distinct while reusing the generic
        # pi-bracket logic.
        evaluations = [evaluate(PI_LO), evaluate(PI_HI)]
        kappa_lo = min(value[0] for value in evaluations) - SAFETY_PAD
        kappa_hi = max(value[1] for value in evaluations) + SAFETY_PAD
        mass_bounds = {
            particle: (
                min(value[2][particle] for value in evaluations) - SAFETY_PAD,
                max(value[2][f"{particle}__upper"] for value in evaluations)
                + SAFETY_PAD,
            )
            for particle in MASS_ORDER
        }
        return _certificate_payload(kappa_lo, kappa_hi, mass_bounds)


def rectangle_certificate(
    readout: Mapping[str, Any],
    endpoint: Mapping[str, Any],
    bridge: Mapping[str, Any],
) -> dict[str, Any]:
    """Enclose the independent gap-by-payload rectangle interval."""

    with localcontext() as context:
        context.prec = PRECISION
        ratios, witness = _shape_inputs(readout)
        reference = bridge["reference_decomposition_compare_only"]
        alpha_inv = _d(reference["alpha_inv_0"])
        delta_top = _d(reference["Delta_top"])
        delta_lep_3loop = _d(reference["Delta_lep"])
        a0 = _d(endpoint["transport_split"]["a0_anchor_inv_alpha"])
        gaps = tuple(
            _d(value)
            for value in endpoint["compare_only"][
                "same_scheme_anchor_gap_interval_inv_alpha"
            ]
        )
        hadronic = _d(endpoint["inputs"]["delta_alpha_had_5_MZ"])
        hadronic_uncertainty = _d(
            endpoint["inputs"]["delta_alpha_had_5_MZ_uncertainty"]
        )

        def evaluate(
            pi: Decimal,
        ) -> tuple[Decimal, Decimal, dict[str, Decimal]]:
            packet_witness = _lepton_packet(pi, witness[0], ratios)
            remainder = delta_lep_3loop - packet_witness / alpha_inv
            if remainder <= 0:
                raise ValueError("higher-order charged-lepton remainder must be positive")

            def solve(
                gap: Decimal,
                hadronic_value: Decimal,
                higher_order: Decimal,
                kernel_slack: Decimal,
            ) -> Decimal:
                target = (
                    Decimal(1)
                    - hadronic_value
                    - delta_top
                    - (a0 + gap) / alpha_inv
                    - higher_order
                )
                packet = target * alpha_inv + kernel_slack
                mass_e = _invert_packet(pi, packet, ratios)
                return (mass_e / witness[0]).ln()

            kappa_lo = solve(
                gaps[0],
                hadronic - hadronic_uncertainty,
                Decimal(0),
                KERNEL_TRUNCATION_PACKET,
            )
            kappa_hi = solve(
                gaps[1],
                hadronic + hadronic_uncertainty,
                Decimal(2) * remainder,
                -KERNEL_TRUNCATION_PACKET,
            )
            factors = (Decimal(1), ratios[0], ratios[1])
            return kappa_lo, kappa_hi, {
                particle: witness[0] * factor * kappa_lo.exp()
                for particle, factor in zip(MASS_ORDER, factors, strict=True)
            } | {
                f"{particle}__upper": witness[0] * factor * kappa_hi.exp()
                for particle, factor in zip(MASS_ORDER, factors, strict=True)
            }

        evaluations = [evaluate(PI_LO), evaluate(PI_HI)]
        kappa_lo = min(value[0] for value in evaluations) - SAFETY_PAD
        kappa_hi = max(value[1] for value in evaluations) + SAFETY_PAD
        mass_bounds = {
            particle: (
                min(value[2][particle] for value in evaluations) - SAFETY_PAD,
                max(value[2][f"{particle}__upper"] for value in evaluations)
                + SAFETY_PAD,
            )
            for particle in MASS_ORDER
        }
        return _certificate_payload(kappa_lo, kappa_hi, mass_bounds)


def interval_as_floats(values: list[str]) -> list[float]:
    """Convert an outward decimal pair to its stable JSON number display."""

    return [float(value) for value in values]


__all__ = [
    "coherent_certificate",
    "interval_as_floats",
    "rectangle_certificate",
]
