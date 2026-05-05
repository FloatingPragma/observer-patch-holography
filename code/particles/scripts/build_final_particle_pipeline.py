#!/usr/bin/env python3
"""Run the simplified final particle-prediction pipeline.

This is the source-tree companion to compute_current_output_table.py.  It does
not try to regenerate every historical exploratory artifact.  It refreshes the
current constructive contracts, the public status/provenance ledgers, and the
final prediction bundle in dependency order.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import subprocess
import sys


CODE_ROOT = Path(__file__).resolve().parents[2]


@dataclass(frozen=True)
class PipelineStep:
    id: str
    script: str
    description: str

    def command(self) -> list[str]:
        return [sys.executable, self.script]


FINAL_PIPELINE_STEPS: tuple[PipelineStep, ...] = (
    PipelineStep(
        "p_closure_trunk",
        "P_derivation/emit_p_closure_trunk.py",
        "Refresh the compressed candidate P-trunk.",
    ),
    PipelineStep(
        "thomson_endpoint_contract",
        "P_derivation/thomson_endpoint_contract.py",
        "Refresh the Ward-projected Thomson endpoint contract.",
    ),
    PipelineStep(
        "rg_matching_threshold_contract",
        "P_derivation/rg_matching_threshold_contract.py",
        "Refresh the RG/matching/threshold contract.",
    ),
    PipelineStep(
        "hadron_spectral_measure_contract",
        "particles/hadron/derive_ward_projected_spectral_measure_contract.py",
        "Refresh the closed-out-of-scope hadronic spectral-measure contract.",
    ),
    PipelineStep(
        "direct_top_bridge_contract",
        "particles/calibration/derive_direct_top_bridge_contract.py",
        "Refresh the direct-top codomain conversion contract.",
    ),
    PipelineStep(
        "charged_end_to_end_nonclosure",
        "particles/leptons/derive_charged_end_to_end_impossibility_theorem.py",
        "Refresh the charged-lepton end-to-end nonclosure theorem packet.",
    ),
    PipelineStep(
        "derivation_gap_ledger",
        "particles/scripts/build_derivation_gap_ledger.py",
        "Refresh the coupled derivation-gap ledger.",
    ),
    PipelineStep(
        "results_status_table",
        "particles/scripts/build_results_status_table.py",
        "Refresh the current public results table.",
    ),
    PipelineStep(
        "exact_nonhadron_mass_bundle",
        "particles/scripts/build_exact_nonhadron_mass_bundle.py",
        "Refresh the exact non-hadron mass bundle.",
    ),
    PipelineStep(
        "exact_fit_surface",
        "particles/scripts/build_exact_fit_surface.py",
        "Refresh the exact-fit-only surface.",
    ),
    PipelineStep(
        "pipeline_closure_status_bootstrap",
        "particles/scripts/build_particle_pipeline_closure_status.py",
        "Refresh closure status before provenance rebuild.",
    ),
    PipelineStep(
        "blind_prediction_provenance",
        "particles/scripts/build_blind_prediction_provenance.py",
        "Refresh row-level blind-prediction provenance.",
    ),
    PipelineStep(
        "pipeline_closure_status_finalize",
        "particles/scripts/build_particle_pipeline_closure_status.py",
        "Refresh closure status after provenance rebuild.",
    ),
    PipelineStep(
        "final_end_to_end_predictions",
        "particles/scripts/build_final_end_to_end_predictions.py",
        "Refresh the final current prediction bundle.",
    ),
    PipelineStep(
        "derivation_chain_closure_matrix",
        "particles/scripts/build_derivation_chain_closure_matrix.py",
        "Refresh the derivation-chain closure matrix.",
    ),
    PipelineStep(
        "mass_derivation_svg",
        "particles/scripts/generate_mass_derivation_svg.py",
        "Refresh the rendered mass-derivation graph.",
    ),
)


def build_steps() -> tuple[PipelineStep, ...]:
    return FINAL_PIPELINE_STEPS


def _run_step(step: PipelineStep, *, verbose: bool) -> None:
    completed = subprocess.run(
        step.command(),
        cwd=CODE_ROOT,
        check=False,
        text=True,
        capture_output=not verbose,
    )
    if completed.returncode != 0:
        if completed.stdout:
            sys.stdout.write(completed.stdout)
        if completed.stderr:
            sys.stderr.write(completed.stderr)
        raise subprocess.CalledProcessError(
            completed.returncode,
            step.command(),
            output=completed.stdout,
            stderr=completed.stderr,
        )


def _dry_run_payload() -> dict[str, object]:
    return {
        "artifact": "oph_final_particle_pipeline_plan",
        "cwd": str(CODE_ROOT),
        "steps": [
            {
                "id": step.id,
                "script": step.script,
                "description": step.description,
            }
            for step in build_steps()
        ],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the simplified final particle-prediction pipeline.")
    parser.add_argument("--dry-run", action="store_true", help="Print the ordered steps without executing them.")
    parser.add_argument("--verbose", action="store_true", help="Stream each builder's output.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.dry_run:
        print(json.dumps(_dry_run_payload(), indent=2, sort_keys=True))
        return 0

    for step in build_steps():
        print(f"[particle-final] {step.id}: {step.description}")
        _run_step(step, verbose=args.verbose)
    print("[particle-final] complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
