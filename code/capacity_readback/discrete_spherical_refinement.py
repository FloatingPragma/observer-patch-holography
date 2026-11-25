"""Conditional level selection and capacity sensitivities, with no physical promotion.

Run from the repository root:
    python3 code/capacity_readback/discrete_spherical_refinement.py

The cell-to-area identification is a hypothesis. The selected level uses the
observational benchmark. Decimal precision controls arithmetic only.
"""

from decimal import Decimal, localcontext
import json


def report() -> dict:
    with localcontext() as ctx:
        ctx.prec = 60
        pi = Decimal("3.14159265358979323846264338327950288419716939937510582097494")
        pixel = Decimal("1.6309682414")
        lp = Decimal("1.616255e-35")
        c = Decimal(299792458)
        mpc = Decimal("3.0856775814913673e22")
        h0 = Decimal("67.36") * 1000 / mpc
        omega_lambda = Decimal("0.6847")
        planck_lambda = 3 * omega_lambda * h0**2 / c**2
        planck_capacity = 3 * pi / (planck_lambda * lp**2)
        capacity = 5 * pixel * Decimal(4)**202

        def fmt(x: Decimal) -> str:
            return format(x, ".17E")

        def change(ratio: Decimal, exponent: Decimal) -> str:
            return fmt(100 * ((exponent * ratio.ln()).exp() - 1))

        powers = {
            "capacity": Decimal(1),
            "Lambda_and_vacuum_density_fixed_constants": Decimal(-1),
            "H_dS_and_a_dS_and_temperature": Decimal("-0.5"),
            "R_dS_and_horizon_time": Decimal("0.5"),
            "conditional_EW_power_law_fixed_pixel": -pixel / 12,
            "conditional_one_tick_multiplier": -Decimal(1) / 48,
        }
        baselines = {
            "rounded_horizon_display": Decimal("3.31e122"),
            "planck2018_central_values": planck_capacity,
        }
        sensitivities = {
            name: {
                "capacity_nats": fmt(baseline),
                "percent_changes": {
                    quantity: change(capacity / baseline, power)
                    for quantity, power in powers.items()
                },
            }
            for name, baseline in baselines.items()
        }
        # The source-forward EW branch has its own pixel and exact bridge
        # definition. This comparison is a substitution diagnostic at fixed P.
        forward_pixel = Decimal("1.630972095858897")
        forward_capacity = 5 * forward_pixel * Decimal(4)**202
        ew_bridge_display = Decimal("3.5321315434e122")
        ew_change = change(forward_capacity / ew_bridge_display, -forward_pixel / 12)
        levels = []
        for level in (201, 202, 203):
            cells = 20 * 4**level
            n = Decimal(cells) * pixel / 4
            levels.append({
                "level": level,
                "cells_exact": str(cells),
                "conditional_capacity_nats": fmt(n),
                "conditional_Lambda_lP2": fmt(3 * pi / n),
            })
        return {
            "status": "CONDITIONAL_OBSERVATION_SELECTED_DISCRETE_REFINEMENT",
            "physical_record_capacity_derived": False,
            "independent_Lambda_prediction": False,
            "frozen_targets_modified": False,
            "inputs": {
                "pixel": str(pixel),
                "pixel_branch": "gauge-width variant, displayed root",
                "pixel_source": "paper/tex_fragments/CLOSURE_PRINCIPLE.tex",
                "capacity_identification": "one carrier cell has area-law allocation P/4",
                "hierarchy": "20*4**L cells; uniform fourfold icosahedral refinement",
                "planck_length_m": str(lp),
                "planck_length_source": "https://physics.nist.gov/cuu/pdf/all.pdf",
                "benchmark_source": "https://arxiv.org/abs/1807.06209",
                "benchmark_H0_km_s_Mpc": "67.36",
                "benchmark_Omega_Lambda": str(omega_lambda),
                "arithmetic_precision_digits": ctx.prec,
                "printed_precision_is_physical_uncertainty": False,
            },
            "levels": levels,
            "level202": {
                "capacity_nats": fmt(capacity),
                "Lambda_lP2": fmt(3 * pi / capacity),
                "Lambda_m_minus_2": fmt(3 * pi / (capacity * lp**2)),
                "R_dS_m": fmt(lp * (capacity / pi).sqrt()),
                "H_dS_s_minus_1": fmt(c / (lp * (capacity / pi).sqrt())),
            },
            "sensitivities": sensitivities,
            "source_forward_EW_substitution_diagnostic": {
                "pixel": str(forward_pixel),
                "bridge_capacity_display_nats": str(ew_bridge_display),
                "level202_capacity_nats": fmt(forward_capacity),
                "v_over_Ecell_percent_change": ew_change,
                "same_capacity_identification_proved": False,
            },
            "unchanged_at_fixed_pixel": [
                "scalar tilt 1-P/48 and the diagnostic clock tilt",
                "fixed-capacity equation of state w=-1",
                "local cell incidence and dimensionless representation data",
                "empirically fitted a0 without a Lambda attachment",
            ],
            "requires_separate_background_calculation": [
                "present-day H0, distances, age, structure growth and CMB transfer",
            ],
        }


if __name__ == "__main__":
    print(json.dumps(report(), indent=2, allow_nan=False))
