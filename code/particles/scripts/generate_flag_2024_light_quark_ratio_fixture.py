#!/usr/bin/env python3
"""Generate the hand-transcribed FLAG 2024 light-quark ratio fixture.

The fixture is compare-only.  FLAG publishes the two marginal uncertainties
used here, but this repository has no covariance for their derived ratio.
Consumers must therefore keep the independent and maximally positively
correlated propagation cases separate.
"""

from __future__ import annotations

import json
from decimal import Decimal, localcontext
from pathlib import Path
from typing import Any


OUTPUT = (
    Path(__file__).resolve().parents[1]
    / "data/flag_2024_light_quark_ratio_fixture.json"
)


def build_average(
    nf: str,
    ms_over_mud: str,
    sigma_ms_over_mud: str,
    ms_over_mud_notation: str,
    mu_over_md: str,
    sigma_mu_over_md: str,
    mu_over_md_notation: str,
) -> dict[str, Any]:
    with localcontext() as context:
        context.prec = 60
        ratio_a = Decimal(ms_over_mud)
        sigma_a = Decimal(sigma_ms_over_mud)
        ratio_u = Decimal(mu_over_md)
        sigma_u = Decimal(sigma_mu_over_md)
        central = ratio_a * (Decimal(1) + ratio_u) / Decimal(2)
        contribution_a = (
            (Decimal(1) + ratio_u) / Decimal(2) * sigma_a
        )
        contribution_u = ratio_a / Decimal(2) * sigma_u
        independent = (
            contribution_a * contribution_a
            + contribution_u * contribution_u
        ).sqrt()
        rho_plus_one = contribution_a + contribution_u
    return {
        "nf": nf,
        "ms_over_mud": {
            "value": ms_over_mud,
            "standard_uncertainty": sigma_ms_over_mud,
            "published_notation": ms_over_mud_notation,
        },
        "mu_over_md": {
            "value": mu_over_md,
            "standard_uncertainty": sigma_mu_over_md,
            "published_notation": mu_over_md_notation,
        },
        "derived_ms_over_md": {
            "value": format(central, "f"),
            "independent_standard_uncertainty": format(independent, "f"),
            "rho_plus_one_standard_uncertainty": format(
                rho_plus_one, "f"
            ),
        },
    }


def build_payload() -> dict[str, Any]:
    return {
        "artifact": "oph_flag_2024_light_quark_ratio_fixture",
        "schema": "oph.flag_2024_light_quark_ratio_fixture.v1",
        "status": "COMPARE_ONLY_HAND_TRANSCRIBED_REFERENCE",
        "source": {
            "publisher": "Flavour Lattice Averaging Group (FLAG)",
            "edition": "FLAG Review 2024, arXiv:2411.04268v3",
            "citation": (
                "Y. Aoki et al., FLAG Review 2024, arXiv:2411.04268; "
                "light-quark averages summarized in Tables 11 and 12, with "
                "the higher-precision Nf=2+1+1 ms/mud recommendation in "
                "Eq. (37)."
            ),
            "url": "https://arxiv.org/abs/2411.04268",
            "transcription_note": (
                "Values are transcribed from the published review. The "
                "repository does not vendor the source PDF or claim a raw "
                "table-payload hash."
            ),
        },
        "derived_quantity": {
            "name": "ms_over_md",
            "identity": "ms/md = (ms/mud) * (1 + mu/md) / 2",
            "input_covariance_available": False,
            "uncertainty_policy": {
                "independent": (
                    "Propagate the two quoted marginal standard uncertainties "
                    "in quadrature."
                ),
                "maximally_positive_correlation": (
                    "Propagate with rho=+1 by adding the two absolute "
                    "first-order contributions."
                ),
                "conservative_rejection_gate": (
                    "Use the larger rho=+1 propagated uncertainty. This is a "
                    "conservative compare-only gate, not a reconstructed FLAG "
                    "covariance."
                ),
                "conservative_rejection_threshold_sigma": 5.0,
            },
        },
        "averages": [
            build_average(
                "2+1+1",
                "27.227",
                "0.081",
                "27.227(81)",
                "0.465",
                "0.024",
                "0.465(24)",
            ),
            build_average(
                "2+1",
                "27.42",
                "0.12",
                "27.42(12)",
                "0.485",
                "0.019",
                "0.485(19)",
            ),
        ],
        "claim_boundary": {
            "comparison_only": True,
            "oph_fit_or_selection_input": False,
            "oph_theory_uncertainty_supplied": False,
            "prediction_preexisted_audit": True,
            "significance_gate_preregistered": False,
            "comparison_is_retrospective": True,
            "note": (
                "These external values may falsify a conditional lane output; "
                "they do not select its Clebsch assignment or supply an OPH "
                "theory uncertainty."
            ),
        },
    }


def canonical_bytes(payload: dict[str, Any]) -> bytes:
    return (json.dumps(payload, indent=2, ensure_ascii=False) + "\n").encode(
        "utf-8"
    )


def main() -> int:
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_bytes(canonical_bytes(build_payload()))
    print(f"wrote {OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
