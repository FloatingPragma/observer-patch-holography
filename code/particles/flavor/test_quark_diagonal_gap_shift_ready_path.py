#!/usr/bin/env python3
"""Ensure quark descent can consume a populated diagonal gap-shift artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[2]
SPREAD_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_spread_map.py"
MEAN_SPLIT_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_sector_mean_split.py"
DESCENT_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_sector_descent.py"
SPREAD_OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_spread_map.json"
DIAGONAL_OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_gap_shift_map.json"


def main() -> int:
    subprocess.run([sys.executable, str(SPREAD_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(MEAN_SPLIT_SCRIPT)], check=True, cwd=ROOT)

    spread = json.loads(SPREAD_OUTPUT.read_text(encoding="utf-8"))
    diagonal = json.loads(DIAGONAL_OUTPUT.read_text(encoding="utf-8"))
    diagonal["tau_u_log_per_side"] = 0.05
    diagonal["tau_d_log_per_side"] = 0.08
    diagonal["proof_status"] = "closed"

    with tempfile.TemporaryDirectory() as tmpdir:
        diag_path = pathlib.Path(tmpdir) / "diagonal_gap_shift.json"
        out_path = pathlib.Path(tmpdir) / "quark_sector_descent.json"
        absent_scalar_path = pathlib.Path(tmpdir) / "absent_scalar.json"
        absent_tau_path = pathlib.Path(tmpdir) / "absent_tau.json"
        absent_emitter_path = pathlib.Path(tmpdir) / "absent_emitter.json"
        absent_b_odd_path = pathlib.Path(tmpdir) / "absent_b_odd.json"
        diag_path.write_text(json.dumps(diagonal, indent=2) + "\n", encoding="utf-8")
        subprocess.run(
            [
                sys.executable,
                str(DESCENT_SCRIPT),
                "--diagonal-gap-shift",
                str(diag_path),
                "--diagonal-gap-shift-scalar-evaluator",
                str(absent_scalar_path),
                "--diagonal-gap-shift-tau-map",
                str(absent_tau_path),
                "--diagonal-gap-shift-emitter",
                str(absent_emitter_path),
                "--b-odd-source-scalar-evaluator",
                str(absent_b_odd_path),
                "--output",
                str(out_path),
            ],
            check=True,
            cwd=ROOT,
        )
        descent = json.loads(out_path.read_text(encoding="utf-8"))

    b_ord = diagonal["B_ord"]
    expected_u = [float(base) + 0.05 * float(shift) for base, shift in zip(spread["E_u_log"], b_ord)]
    expected_d = [float(base) + 0.08 * float(shift) for base, shift in zip(spread["E_d_log"], b_ord)]
    if descent.get("tau_u_log_per_side") != 0.05 or descent.get("tau_d_log_per_side") != 0.08:
        print("descent did not preserve populated diagonal gap-shift coefficients", file=sys.stderr)
        return 1
    if descent.get("even_excitation_proof_status") != "closed":
        print("closed diagonal gap-shift should promote the even-excitation proof status", file=sys.stderr)
        return 1
    if descent.get("tau_value_source_closed") is not True:
        print("closed source map should mark the tau source as closed", file=sys.stderr)
        return 1
    if any(abs(float(a) - float(b)) > 1.0e-12 for a, b in zip(descent["E_u_log"], expected_u)):
        print("up-sector logs did not absorb the diagonal gap shift", file=sys.stderr)
        return 1
    if any(abs(float(a) - float(b)) > 1.0e-12 for a, b in zip(descent["E_d_log"], expected_d)):
        print("down-sector logs did not absorb the diagonal gap shift", file=sys.stderr)
        return 1

    # Merely populating the same arithmetic in a comparison-only tau map must
    # not reproduce the source-closed result.
    with tempfile.TemporaryDirectory() as tmpdir:
        diag_path = pathlib.Path(tmpdir) / "diagonal_gap_shift.json"
        tau_path = pathlib.Path(tmpdir) / "comparison_tau_map.json"
        scalar_path = pathlib.Path(tmpdir) / "comparison_scalar.json"
        out_path = pathlib.Path(tmpdir) / "quark_sector_descent.json"
        absent_emitter_path = pathlib.Path(tmpdir) / "absent_emitter.json"
        absent_b_odd_path = pathlib.Path(tmpdir) / "absent_b_odd.json"
        diag_path.write_text(json.dumps(diagonal, indent=2) + "\n", encoding="utf-8")
        scalar_path.write_text(
            json.dumps(
                {
                    "artifact": "oph_family_excitation_diagonal_gap_shift_scalar_evaluator",
                    "proof_status": "comparison_only_scalar_projection_value_source_open",
                    "smallest_constructive_missing_object": "source_derived_c_d_over_c_u_or_t1_value",
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        tau_path.write_text(
            json.dumps(
                {
                    "artifact": "oph_family_excitation_diagonal_gap_shift_tau_map",
                    "proof_status": "comparison_only_scalar_projection_value_source_open",
                    "source_value_promotion_ready": False,
                    "value_classification": "comparison_only",
                    "tau_u_log_per_side": 0.05,
                    "tau_d_log_per_side": 0.08,
                    "smallest_constructive_missing_object": "source_derived_c_d_over_c_u_or_t1_value",
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        subprocess.run(
            [
                sys.executable,
                str(DESCENT_SCRIPT),
                "--diagonal-gap-shift",
                str(diag_path),
                "--diagonal-gap-shift-scalar-evaluator",
                str(scalar_path),
                "--diagonal-gap-shift-tau-map",
                str(tau_path),
                "--diagonal-gap-shift-emitter",
                str(absent_emitter_path),
                "--b-odd-source-scalar-evaluator",
                str(absent_b_odd_path),
                "--output",
                str(out_path),
            ],
            check=True,
            cwd=ROOT,
        )
        comparison_descent = json.loads(out_path.read_text(encoding="utf-8"))
    if comparison_descent.get("tau_value_source_closed") is not False:
        print("comparison-only tau values must not be source closed", file=sys.stderr)
        return 1
    if comparison_descent.get("even_excitation_proof_status") != "comparison_only_diagonal_gap_shift_value_source_open":
        print("comparison-only tau values must preserve the open source status", file=sys.stderr)
        return 1
    if comparison_descent.get("exact_missing_object") != "source_derived_c_d_over_c_u_or_t1_value":
        print("comparison-only tau values should preserve the D12 value-source gap", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
