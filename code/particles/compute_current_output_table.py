#!/usr/bin/env python3
"""Build the current particle status surface in a disposable runtime workspace.

This wrapper keeps `reverse-engineering-reality/code/particles` source-only.
It stages a temporary working copy under the top-level `temp/particles_runtime`
directory, runs the active build chain there, and writes the latest rendered
status outputs into `temp/particles_runtime/current`.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path


CODE_ROOT = Path(__file__).resolve().parent
RER_ROOT = CODE_ROOT.parents[1]
WORKSPACE_ROOT = CODE_ROOT.parents[2]
DEFAULT_RUNTIME_ROOT = WORKSPACE_ROOT / "temp" / "particles_runtime"

EXCLUDE_NAMES = {
    "__pycache__",
    ".pytest_cache",
    "runs",
    "RESULTS_STATUS.md",
    "results_status.json",
    "particle_mass_derivation_graph.svg",
    "HADRON_SYSTEMATICS_STATUS.md",
}

CURRENT_QUARK_FORWARD = {
    "artifact": "oph_forward_yukawas",
    "public_surface_candidate_allowed": True,
    "source_mode": "factorized_descent",
    "singular_values_u": [
        0.002287663517422587,
        1.2566920642121726,
        172.38864559527133,
    ],
    "singular_values_d": [
        0.005121438948009003,
        0.0916111424053981,
        3.82011787992005,
    ],
}

CURRENT_QUARK_MEAN_SPLIT = {
    "artifact": "oph_quark_sector_mean_split",
    "active_candidate": "ordered_affine_mean_readout_candidate",
}


def _ignore(_src: str, names: list[str]) -> set[str]:
    return {name for name in names if name in EXCLUDE_NAMES}


def _run(cmd: list[str], *, cwd: Path) -> None:
    subprocess.run(cmd, check=True, cwd=cwd)


def _copy_outputs(work_particles: Path, current_dir: Path) -> None:
    current_dir.mkdir(parents=True, exist_ok=True)
    outputs = [
        "RESULTS_STATUS.md",
        "results_status.json",
        "particle_mass_derivation_graph.svg",
        "runs/status/status_table_forward_current.json",
    ]
    for rel in outputs:
        src = work_particles / rel
        dst = current_dir / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


def _seed_curated_quark_surface(work_particles: Path) -> None:
    flavor_dir = work_particles / "runs" / "flavor"
    flavor_dir.mkdir(parents=True, exist_ok=True)
    (flavor_dir / "forward_yukawas.json").write_text(
        json.dumps(CURRENT_QUARK_FORWARD, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (flavor_dir / "quark_sector_mean_split.json").write_text(
        json.dumps(CURRENT_QUARK_MEAN_SPLIT, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def build_runtime(runtime_root: Path, *, with_hadrons: bool) -> None:
    work_code = runtime_root / "work" / "code"
    work_particles = work_code / "particles"
    current_dir = runtime_root / "current"

    if work_code.exists():
        shutil.rmtree(work_code)
    work_code.mkdir(parents=True, exist_ok=True)
    shutil.copytree(CODE_ROOT, work_particles, ignore=_ignore, dirs_exist_ok=True)

    build_steps = [
        ["python3", "particles/calibration/derive_d10_ew_observable_family.py"],
        ["python3", "particles/calibration/derive_d10_ew_source_transport_pair.py"],
        ["python3", "particles/calibration/derive_d10_ew_source_transport_readout.py"],
        ["python3", "particles/calibration/derive_d10_ew_population_evaluator.py"],
        ["python3", "particles/calibration/derive_d10_ew_exact_closure_beyond_current_carrier.py"],
        ["python3", "particles/calibration/derive_d10_ew_fiberwise_population_tree_law_beneath_single_tree_identity.py"],
        ["python3", "particles/calibration/derive_d10_ew_tau2_current_carrier_obstruction.py"],
        ["python3", "particles/calibration/derive_d10_ew_exact_wz_coordinate_beyond_single_tree_identity.py"],
        ["python3", "particles/calibration/derive_d10_ew_exact_mass_pair_chart_current_carrier.py"],
        ["python3", "particles/calibration/derive_d10_ew_repair_branch_beyond_current_carrier.py"],
        ["python3", "particles/calibration/derive_d10_ew_repair_target_point_diagnostic.py"],
        ["python3", "particles/calibration/derive_d10_ew_w_anchor_neutral_shear_factorization.py"],
        ["python3", "particles/calibration/derive_d10_ew_minimal_conditional_promotion.py"],
        ["python3", "particles/calibration/derive_d10_ew_target_emitter_candidate.py"],
        ["python3", "particles/calibration/derive_d10_ew_source_transport_readout.py"],
        ["python3", "particles/calibration/derive_d10_ew_exactness_audit.py"],
        ["python3", "particles/calibration/derive_d11_critical_surface_readout.py"],
        ["python3", "particles/calibration/derive_d11_forward_seed.py"],
        ["python3", "particles/calibration/derive_d11_forward_seed_promotion_certificate.py"],
    ]

    if with_hadrons:
        build_steps.extend(
            [
                ["python3", "particles/qcd/derive_lambda_msbar_descendant.py"],
                ["python3", "particles/hadron/derive_full_unquenched_correlator.py"],
                ["python3", "particles/hadron/derive_runtime_schedule_receipt_n_therm_and_n_sep.py"],
                ["python3", "particles/hadron/derive_stable_channel_cfg_source_measure_payload.py"],
                ["python3", "particles/hadron/derive_stable_channel_sequence_population.py"],
                ["python3", "particles/hadron/derive_stable_channel_sequence_evaluation.py"],
                ["python3", "particles/hadron/derive_current_hadron_lane_audit.py"],
            ]
        )

    for step in build_steps:
        _run(step, cwd=work_code)

    _seed_curated_quark_surface(work_particles)

    status_cmd = ["python3", "particles/scripts/build_results_status_table.py"]
    svg_cmd = ["python3", "particles/scripts/generate_mass_derivation_svg.py"]
    if with_hadrons:
        status_cmd.append("--with-hadrons")
        svg_cmd.append("--with-hadrons")
    _run(status_cmd, cwd=work_code)
    _run(svg_cmd, cwd=work_code)
    _copy_outputs(work_particles, current_dir)


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the current particle status surface in a disposable runtime workspace.")
    parser.add_argument("--runtime-root", default=str(DEFAULT_RUNTIME_ROOT))
    parser.add_argument("--with-hadrons", action="store_true")
    args = parser.parse_args()

    runtime_root = Path(args.runtime_root).resolve()
    build_runtime(runtime_root, with_hadrons=args.with_hadrons)
    print(f"runtime work tree: {runtime_root / 'work' / 'code' / 'particles'}")
    print(f"current outputs: {runtime_root / 'current'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
