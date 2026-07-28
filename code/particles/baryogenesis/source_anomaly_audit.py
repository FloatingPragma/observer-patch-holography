#!/usr/bin/env python3
"""Exact finite audit for the OPH baryogenesis source theorem.

This script checks only representation/anomaly arithmetic.  It does not
construct a physical repair generator or use the baryon abundance to choose
an attachment.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import json
from pathlib import Path
from typing import Callable, Dict

N_G = 3
N_C = 3


@dataclass(frozen=True)
class WeakDoublet:
    name: str
    multiplicity: int
    hypercharge: Fraction
    baryon: Fraction
    lepton: Fraction


DOUBLETS = (
    WeakDoublet("Q_L", N_C, Fraction(1, 6), Fraction(1, 3), Fraction(0)),
    WeakDoublet("L_L", 1, Fraction(-1, 2), Fraction(0), Fraction(1)),
)


def anomaly_per_generation(charge: Callable[[WeakDoublet], Fraction]) -> Fraction:
    """Return sum r_psi * 2T_2(R_psi); 2T_2(doublet)=1."""
    return sum((d.multiplicity * charge(d) for d in DOUBLETS), Fraction(0))


attachments: Dict[str, Callable[[WeakDoublet], Fraction]] = {
    "Y": lambda d: d.hypercharge,
    "B-L": lambda d: d.baryon - d.lepton,
    "B": lambda d: d.baryon,
    "L": lambda d: d.lepton,
    "B+L": lambda d: d.baryon + d.lepton,
    "unit_left_doublet": lambda d: Fraction(1),
}

results = {}
for name, charge in attachments.items():
    per_gen = anomaly_per_generation(charge)
    total = N_G * per_gen
    results[name] = {
        "per_generation": str(per_gen),
        "k_R_Ng3": str(total),
        "k_R_Ng3_float": float(total),
    }

assert Fraction(results["Y"]["k_R_Ng3"]) == 0
assert Fraction(results["B-L"]["k_R_Ng3"]) == 0
assert Fraction(results["B"]["k_R_Ng3"]) == 3
assert Fraction(results["L"]["k_R_Ng3"]) == 3
assert Fraction(results["B+L"]["k_R_Ng3"]) == 6
assert Fraction(results["unit_left_doublet"]["k_R_Ng3"]) == 12

receipt = {
    "theorem": "finite-quotient anomaly-and-current theorem",
    "scope": "mixed SU(2)_L^2-U(1)_R anomaly arithmetic only",
    "inputs": {
        "N_g": N_G,
        "N_c": N_C,
        "weak_doublets": [
            {
                "name": d.name,
                "multiplicity": d.multiplicity,
                "Y": str(d.hypercharge),
                "B": str(d.baryon),
                "L": str(d.lepton),
                "two_T2": "1",
            }
            for d in DOUBLETS
        ],
    },
    "outputs": results,
    "conclusions": {
        "direct_hypercharge_deck_attachment": "k_R=0 (no-go)",
        "B_plus_L_attachment": "k_R=6 (conditional attachment)",
        "normal_form_selects_sign": False,
        "repair_generator_constructed": False,
    },
}

out = Path(__file__).with_name("source_anomaly_receipt.json")
out.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
print(json.dumps(receipt, indent=2, sort_keys=True))
print(f"\nWrote {out}")
