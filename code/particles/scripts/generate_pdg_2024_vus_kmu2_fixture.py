#!/usr/bin/env python3
"""Generate the compare-only PDG ``|V_us|`` coordinate used by the axis no-go.

The value is the explicitly identified K_mu2-channel determination in the
2024 PDG CKM review.  It is not represented as the global CKM-fit value and it
is never an input to the construction or selection of icosahedral axes.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


OUTPUT = (
    Path(__file__).resolve().parents[1]
    / "data/pdg_2024_vus_kmu2_fixture.json"
)


def build_payload() -> dict[str, Any]:
    return {
        "artifact": "oph_pdg_2024_vus_kmu2_fixture",
        "schema": "oph.pdg_2024_vus_kmu2_fixture.v1",
        "status": "COMPARE_ONLY_HAND_TRANSCRIBED_REFERENCE",
        "source": {
            "publisher": "Particle Data Group",
            "edition": (
                "Review of Particle Physics 2024, CKM Quark-Mixing Matrix, "
                "revised April 2024; PDF dated 31 May 2024"
            ),
            "citation": (
                "S. Navas et al. (Particle Data Group), Phys. Rev. D 110, "
                "030001 (2024), Sec. 12.2.2: the K -> mu nu(gamma) over "
                "pi -> mu nu(gamma) determination using the lattice-QCD "
                "decay-constant ratio"
            ),
            "url": (
                "https://pdg.lbl.gov/2024/reviews/"
                "rpp2024-rev-ckm-matrix.pdf"
            ),
            "transcription_note": (
                "The published central value and quoted uncertainty are "
                "transcribed into this deterministic fixture. The source PDF "
                "is not vendored or assigned an invented raw-payload hash."
            ),
        },
        "coordinate": {
            "name": "abs_Vus",
            "determination": "Kmu2_decay_constant_ratio",
            "value": "0.2250",
            "standard_uncertainty": "0.0004",
            "published_notation": "0.2250 +/- 0.0004",
            "uncertainty_semantics": (
                "quoted PDG uncertainty for this channel-specific "
                "determination; no covariance or likelihood is reconstructed"
            ),
        },
        "claim_boundary": {
            "comparison_only": True,
            "used_to_construct_or_select_axes": False,
            "global_ckm_fit_value": False,
            "oph_fit_or_selection_input": False,
            "note": (
                "This channel-specific coordinate only checks the scale of "
                "the exact residual-axis gap. It is not the PDG global-fit "
                "value and does not enter the icosahedral enumeration."
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
