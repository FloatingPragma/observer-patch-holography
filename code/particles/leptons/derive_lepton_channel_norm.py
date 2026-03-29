#!/usr/bin/env python3
"""Expose the charged-lepton channel norm as a separate artifact."""

from __future__ import annotations

import argparse
import json
import pathlib
from datetime import datetime, timezone


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "leptons" / "lepton_descent_tensor.json"
DEFAULT_CHARGED_BUDGET = ROOT / "particles" / "runs" / "flavor" / "charged_budget_transport.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "leptons" / "lepton_channel_norm.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _stream_map(stream: list[dict[str, object]]) -> dict[str, float]:
    mapped: dict[str, float] = {}
    for item in stream:
        refinement = str(item.get("refinement", "snapshot"))
        value = item.get("value")
        if value is None:
            continue
        mapped[refinement] = float(value)
    return mapped


def _stream_list(mapped: dict[str, float]) -> list[dict[str, float | str]]:
    return [{"refinement": refinement, "value": float(value)} for refinement, value in sorted(mapped.items())]


def main() -> int:
    parser = argparse.ArgumentParser(description="Build a charged-lepton channel-norm artifact.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT), help="Input charged-lepton channel-tensor artifact.")
    parser.add_argument(
        "--charged-budget",
        default=str(DEFAULT_CHARGED_BUDGET),
        help="Optional shared charged-budget artifact.",
    )
    parser.add_argument("--output", default=str(DEFAULT_OUT), help="Output channel-norm JSON path.")
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    summary = dict(payload.get("sector_response_summary", {}))
    raw_candidate = summary.get("raw_channel_norm_candidate")
    raw_candidate = None if raw_candidate is None else float(raw_candidate)
    refinement_stream = list(summary.get("raw_channel_norm_by_refinement", []))
    scale_scope_candidate = summary.get("scale_scope_candidate")
    shared_budget_key = summary.get("shared_budget_key")
    charged_budget_path = pathlib.Path(args.charged_budget)
    charged_budget_payload = (
        json.loads(charged_budget_path.read_text(encoding="utf-8")) if charged_budget_path.exists() else {}
    )
    charged_budget_proof_status = str(charged_budget_payload.get("proof_status", "open"))
    charged_scalarization_law = dict(charged_budget_payload.get("charged_dirac_scalarization_certificate", {}))
    beta_e_stream = list(dict(charged_budget_payload.get("beta_by_sector_by_refinement", {})).get("e", []))
    charged_budget_total_stream = list(charged_budget_payload.get("B_ch_by_refinement", []))
    beta_e_map = _stream_map(beta_e_stream)
    charged_budget_total_map = _stream_map(charged_budget_total_stream)
    common_refinements = sorted(set(beta_e_map) & set(charged_budget_total_map))
    g_e_candidate_stream = {
        refinement: beta_e_map[refinement] * charged_budget_total_map[refinement]
        for refinement in common_refinements
    }
    g_e_candidate = None
    if common_refinements:
        g_e_candidate = float(g_e_candidate_stream[common_refinements[-1]])
    elif raw_candidate is not None:
        g_e_candidate = raw_candidate

    proof_status = "open"
    scale_scope = "sector_local_unproven"
    shared_budget_share_e = None
    closure_route = None
    channel_norm_closed = False
    g_e = None
    g_e_by_refinement = refinement_stream

    if charged_budget_payload:
        proof_status = charged_budget_proof_status if charged_budget_proof_status in {
            "shared_budget_only",
            "shared_budget_closed",
            "sector_local_closed",
        } else "shared_budget_only"
        scale_scope = "shared_charged_budget"
        closure_route = "shared_charged_budget"
        shared_budget_share_e = beta_e_map[common_refinements[-1]] if common_refinements else None
        g_e_by_refinement = _stream_list(g_e_candidate_stream)
        if proof_status == "shared_budget_closed" and g_e_candidate is not None:
            g_e = g_e_candidate
            channel_norm_closed = True
    elif scale_scope_candidate == "shared_budget_only":
        proof_status = "shared_budget_only"
        scale_scope = "shared_budget_only"
        shared_budget_share_e = raw_candidate
    elif scale_scope_candidate == "sector_local_closed":
        proof_status = "sector_local_closed"
        scale_scope = "sector_local"
        closure_route = "sector_local"
        g_e = raw_candidate
        channel_norm_closed = g_e is not None
    else:
        proof_status = "open"
        scale_scope = "sector_local_unproven"

    artifact = {
        "artifact": "oph_lepton_channel_norm",
        "generated_utc": _timestamp(),
        "labels": payload.get("labels"),
        "channel": payload.get("channel", "L-H-E"),
        "g_e_candidate": g_e_candidate,
        "g_e_by_refinement": g_e_by_refinement,
        "g_e": g_e,
        "channel_norm_closed": channel_norm_closed,
        "scale_scope": scale_scope,
        "proof_status": proof_status,
        "closure_route": closure_route,
        "shared_budget_key": shared_budget_key,
        "shared_budget_share_e": shared_budget_share_e,
        "beta_e_by_refinement": beta_e_stream,
        "charged_budget_total_by_refinement": charged_budget_total_stream,
        "channel_norm_refinement_certificate": {
            "status": "snapshot_only" if len(g_e_by_refinement) <= 1 else "candidate_stream",
            "samples": len(g_e_by_refinement),
        },
        "metadata": {
            "channel_tensor_artifact": payload.get("artifact", "unknown"),
            "charged_budget_transport_artifact": charged_budget_payload.get("artifact"),
            "charged_scalarization_law": charged_scalarization_law,
            "note": "Channel norm is kept explicit. Shared charged-budget and sector-local routes stay distinct until a theorem closes one of them.",
            **dict(payload.get("metadata", {})),
        },
    }

    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
