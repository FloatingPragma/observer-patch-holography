#!/usr/bin/env python3
"""Audit whether the declared OPH MaxEnt branch forces Koide balance.

The result is a conditional theorem plus a countermodel.  It separates the
algebraic Z3 identity from the additional ensemble and deterministic-readout
bridges needed to select the balanced carrier.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = (
    ROOT / "particles" / "runs" / "leptons"
    / "charged_z3_maxent_balance_audit.json"
)


def effective_ratio(m_zero: float, m_charged: float, action_gap: float) -> float:
    """Return w_charged/w_zero for S_charged-S_zero=action_gap."""

    return (m_charged / m_zero) * math.exp(-action_gap)


def koide_from_power_ratio(charged_over_singlet: float) -> float:
    """For Z3 square-root masses, Q=(1+E_charged/E_singlet)/3."""

    return (1.0 + charged_over_singlet) / 3.0


def build_artifact() -> dict[str, Any]:
    natural_ratio = effective_ratio(1.0, 2.0, 0.0)
    compensated_ratio = effective_ratio(1.0, 2.0, math.log(2.0))
    return {
        "artifact": "oph_charged_z3_maxent_balance_audit",
        "status": "CLOSED_CONDITIONAL_THEOREM_CURRENT_OPH_INPUT_OPEN",
        "public_koide_promotion_allowed": False,
        "algebraic_theorem": {
            "carrier": "Phi=a I+b P+conj(b) P^2, P^3=I",
            "singlet_power": "E_0=3 a^2",
            "charged_power": "E_c=6 |b|^2",
            "koide_identity": "Q=(1+E_c/E_0)/3=(1+2|b/a|^2)/3",
            "balance_equivalence": "Q=2/3 iff E_c=E_0 iff |b|/a=1/sqrt(2)",
        },
        "conditional_maxent_theorem": {
            "premises": [
                "the MaxEnt variables are the two normalized block powers p_j=E_j/(E_0+E_c)",
                "the relative-MaxEnt base weights m_j and block actions S_j are source-emitted",
                "the ensemble probabilities descend to deterministic charged-Yukawa block powers",
            ],
            "law": "E_c/E_0=(m_c/m_0)*exp(-(S_c-S_0))",
            "conclusion": "Koide follows exactly iff the two effective block weights are equal",
        },
        "natural_trace_countermodel": {
            "real_irrep_dimensions": {"singlet": 1, "charged": 2},
            "base_weight_ratio_m_c_over_m_0": 2,
            "equal_action_gap": 0,
            "emitted_power_ratio": natural_ratio,
            "emitted_Q": koide_from_power_ratio(natural_ratio),
            "conclusion": (
                "Trace/Haar/Lebesgue MaxEnt counts the charged doublet twice and does not emit Koide balance."
            ),
        },
        "smallest_balance_input": {
            "required_action_gap_S_c_minus_S_0": "ln(2)",
            "effective_power_ratio": compensated_ratio,
            "emitted_Q": koide_from_power_ratio(compensated_ratio),
            "equivalent_option": "declare equal base weights on the two coarse blocks",
        },
        "current_oph_boundary": {
            "one_bit_degeneracy_compensation_derived": False,
            "maxent_probability_to_yukawa_power_bridge_derived": False,
            "twelve_twenty_four_attachment": (
                "The oriented 24-slot register supplies a bookkeeping doubling, but no current theorem "
                "attaches it as an ln(2) action penalty on the charged Z3 block."
            ),
            "claim": (
                "Current OPH MaxEnt symmetry inheritance preserves Z3 invariance but does not select "
                "the balanced singlet/charged power ratio."
            ),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": artifact["status"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
