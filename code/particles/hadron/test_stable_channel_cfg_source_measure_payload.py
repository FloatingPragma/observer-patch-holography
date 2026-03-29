#!/usr/bin/env python3
"""Smoke-test the stable-channel cfg/source measure payload artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[2]
FULL_SCRIPT = ROOT / "particles" / "hadron" / "derive_full_unquenched_correlator.py"
POP_SCRIPT = ROOT / "particles" / "hadron" / "derive_stable_channel_sequence_population.py"
SCRIPT = ROOT / "particles" / "hadron" / "derive_stable_channel_cfg_source_measure_payload.py"
RECEIPT_SCRIPT = ROOT / "particles" / "hadron" / "derive_runtime_schedule_receipt_n_therm_and_n_sep.py"
def main() -> int:
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = pathlib.Path(tmp_dir)
        payload_path = tmp_path / "payload.json"
        receipt_path = tmp_path / "receipt.json"
        subprocess.run([sys.executable, str(FULL_SCRIPT)], check=True, cwd=ROOT)
        subprocess.run([sys.executable, str(POP_SCRIPT)], check=True, cwd=ROOT)
        subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--runtime-receipt",
                str(tmp_path / "missing_receipt.json"),
                "--output",
                str(payload_path),
            ],
            check=True,
            cwd=ROOT,
        )
        subprocess.run(
            [
                sys.executable,
                str(RECEIPT_SCRIPT),
                "--payload",
                str(payload_path),
                "--output",
                str(receipt_path),
            ],
            check=True,
            cwd=ROOT,
        )
        subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--runtime-receipt",
                str(receipt_path),
                "--output",
                str(payload_path),
            ],
            check=True,
            cwd=ROOT,
        )
        payload = json.loads(payload_path.read_text(encoding="utf-8"))
    if payload.get("artifact") != "oph_hadron_stable_channel_cfg_source_measure_payload":
        print("wrong stable-channel cfg/source payload artifact id", file=sys.stderr)
        return 1
    if payload.get("status") != "law_closed_waiting_measure_realization":
        print("cfg/source payload artifact should wait on measure realization", file=sys.stderr)
        return 1
    ensembles = payload.get("ensemble_payloads", [])
    if not ensembles:
        print("cfg/source payload artifact should expose per-ensemble payload shells", file=sys.stderr)
        return 1
    if "cfg_source_corr_t" not in ensembles[0]["pi_iso"]:
        print("pi channel payload should expose cfg/source correlator arrays", file=sys.stderr)
        return 1
    if "cfg_source_corr_direct_t" not in ensembles[0]["N_iso"] or "cfg_source_corr_exchange_t" not in ensembles[0]["N_iso"]:
        print("N channel payload should expose direct and exchange cfg/source correlator arrays", file=sys.stderr)
        return 1
    if payload.get("smallest_constructive_missing_object") != "runtime_schedule_receipt_N_therm_and_N_sep":
        print("cfg/source payload artifact should reduce to the external schedule receipt on the seeded 2+1 family", file=sys.stderr)
        return 1
    if ensembles[0].get("n_cfg") != 2 or ensembles[0].get("n_src_per_cfg") != 2:
        print("cfg/source payload artifact should populate the deterministic 2 cfg x 2 source support contract", file=sys.stderr)
        return 1
    if len(ensembles[0].get("cfg_ids", [])) != 2:
        print("cfg/source payload artifact should emit concrete cfg ids", file=sys.stderr)
        return 1
    if "cfg_realization_contract" not in ensembles[0]:
        print("cfg/source payload artifact should expose the cfg realization contract", file=sys.stderr)
        return 1
    if "cfg_seed_hash_formula" not in ensembles[0]["cfg_realization_contract"]:
        print("cfg/source payload artifact should expose the deterministic cfg seed hash law", file=sys.stderr)
        return 1
    if len(ensembles[0]["cfg_realization_contract"].get("cfg_seed_hashes", {})) != 2:
        print("cfg/source payload artifact should emit concrete cfg seed hashes", file=sys.stderr)
        return 1
    if payload.get("support_realization_schedule", {}).get("status") != "external_receipt_required_before_execution":
        print("cfg/source payload artifact should expose the external runtime receipt gate before schedule execution", file=sys.stderr)
        return 1
    if payload.get("support_realization_schedule", {}).get("runtime_receipt_artifact") != "runtime_schedule_receipt_N_therm_and_N_sep":
        print("cfg/source payload artifact should point at the runtime receipt artifact", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
