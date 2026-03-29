#!/usr/bin/env python3
"""Emit a constructive plan for full isospin-symmetric nucleon contractions."""

from __future__ import annotations

import argparse
import json
import pathlib


def main() -> int:
    ap = argparse.ArgumentParser(description="Write a proton contraction plan artifact.")
    ap.add_argument("--out", default="particles/runs/hadron/proton_contraction_plan.json")
    args = ap.parse_args()

    out_path = pathlib.Path(args.out)
    if not out_path.is_absolute():
        out_path = pathlib.Path(__file__).resolve().parents[3] / out_path
    out_path.parent.mkdir(parents=True, exist_ok=True)

    payload = {
        "artifact": "N_iso.full_baryon_contractions",
        "status": "formula_closed",
        "proof_status": "formula_closed",
        "predictive_promotion_allowed": False,
        "objective": "upgrade the stable-channel N_iso readout from direct-term-only to full baryon contractions",
        "formula_id": "N_local_Cgamma5_direct_minus_exchange",
        "channel_scope": "isospin_symmetric_N_iso",
        "wick_terms": ["direct", "exchange"],
        "color_antisymmetrization": "epsilon_sink_times_epsilon_source",
        "operator_formula": "sum_x eps^{abc}(u^T C gamma5 d)u",
        "parity_projector": "(1+gamma_0)/2",
        "required_tensors": [
            "Q[:,:,a,ap]",
            "Q[:,:,b,bp]",
            "Q[:,:,c,cp]",
            "Q[:,:,c,ap]",
            "Q[:,:,a,cp]",
        ],
        "contraction_rules": {
            "direct": "Gd = Sc * trace(Gamma @ Sb @ GammaB @ Sa.T)",
            "exchange": "Gx = Sca @ GammaB.T @ Sb.T @ Gamma.T @ Sac",
            "full": "Gfull = Gd - Gx",
        },
        "output_sequences": [
            "corr_direct_t",
            "corr_exchange_t",
            "corr_t",
            "am_eff_t",
        ],
        "deliverables": [
            "raw contraction tensors or npz reference",
            "full parity-projected N_iso correlator sequence",
            "direct-minus-exchange bookkeeping",
            "fit-ready metadata",
        ],
        "promotion_gate": "StableChannelForwardWindowConvergence",
        "notes": [
            "This artifact closes the direct-minus-exchange contraction rule for the isospin-symmetric nucleon channel.",
            "What remains open is sequence population on seeded unquenched ensembles, not the algebraic contraction identity itself.",
        ],
    }
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
