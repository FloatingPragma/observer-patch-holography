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
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any


CODE_ROOT = Path(__file__).resolve().parent
RER_ROOT = CODE_ROOT.parents[1]
WORKSPACE_ROOT = CODE_ROOT.parents[2]
DEFAULT_RUNTIME_ROOT = WORKSPACE_ROOT / "temp" / "particles_runtime"
NEUTRINO_COMPARE_ONLY_FIT = Path("runs/neutrino/neutrino_compare_only_scale_fit.json")
NEUTRINO_TWO_PARAMETER_EXACT_ADAPTER = Path("runs/neutrino/neutrino_two_parameter_exact_adapter.json")
NEUTRINO_EXACT_ADAPTER_BRIDGE_COORDINATE = Path("runs/neutrino/neutrino_exact_adapter_bridge_coordinate.json")
NEUTRINO_BRIDGE_RIGIDITY_THEOREM = Path("runs/neutrino/neutrino_bridge_rigidity_theorem.json")
NEUTRINO_ABSOLUTE_ATTACHMENT_THEOREM = Path("runs/neutrino/neutrino_absolute_attachment_theorem.json")

EXCLUDE_NAMES = {
    "__pycache__",
    ".pytest_cache",
    "runs",
    "RESULTS_STATUS.md",
    "results_status.json",
    "particle_mass_derivation_graph.svg",
    "HADRON_SYSTEMATICS_STATUS.md",
}

RUNTIME_SURFACED_ARTIFACTS = (
    Path("runs/calibration/d10_ew_w_anchor_neutral_shear_factorization_official_pdg_2025_update.json"),
    Path("runs/calibration/d11_reference_exact_adapter.json"),
    Path("runs/flavor/forward_yukawas.json"),
    Path("runs/flavor/quark_current_family_exact_readout.json"),
    Path("runs/flavor/quark_current_family_selected_sheet_closure.json"),
    Path("runs/flavor/quark_current_family_quadratic_readout_theorem.json"),
    Path("runs/flavor/quark_current_premise_no_go_theorem.json"),
    Path("runs/flavor/quark_maximal_theorem_emitted_package.json"),
    Path("runs/flavor/quark_minimal_extension_closure_theorem.json"),
    Path("runs/flavor/light_quark_isospin_overlap_defect_selector_law.json"),
    Path("runs/flavor/quark_d12_t1_value_law.json"),
    Path("runs/flavor/quark_physical_branch_repair_theorem.json"),
    Path("runs/flavor/quark_relative_sheet_selector.json"),
    Path("runs/flavor/quark_sector_mean_split.json"),
    Path("runs/leptons/forward_charged_leptons.json"),
    Path("runs/leptons/lepton_current_family_affine_anchor_theorem.json"),
    Path("runs/leptons/lepton_current_family_exact_readout.json"),
    Path("runs/leptons/lepton_current_family_quadratic_readout_theorem.json"),
    Path("runs/neutrino/forward_neutrino_closure_bundle.json"),
    Path("runs/neutrino/neutrino_bridge_rigidity_theorem.json"),
    Path("runs/neutrino/neutrino_absolute_attachment_theorem.json"),
    Path("runs/neutrino/exact_blocking_items.json"),
    Path("runs/neutrino/neutrino_weighted_cycle_theorem_object.json"),
    Path("runs/neutrino/realized_same_label_gap_defect_readback.json"),
    Path("runs/neutrino/same_label_overlap_defect_weight_normalizer.json"),
    Path("runs/neutrino/same_label_scalar_certificate.json"),
    Path("runs/neutrino/neutrino_weighted_cycle_shared_basis_representation.json"),
    Path("runs/neutrino/neutrino_physical_majorana_phase_theorem.json"),
    Path("runs/neutrino/pmns_from_shared_basis.json"),
    Path("runs/neutrino/shared_charged_lepton_left_basis.json"),
    Path("runs/neutrino/intrinsic_neutrino_eta_demo_payload.json"),
    Path("runs/neutrino/intrinsic_neutrino_exact_eta_map.json"),
    Path("runs/neutrino/intrinsic_neutrino_exact_mixing_law_validation.json"),
    Path("runs/neutrino/intrinsic_neutrino_mass_eigenstate_bundle_from_scalar_certificate.json"),
    Path("runs/neutrino/majorana_overlap_defect_hessian.json"),
    Path("runs/neutrino/majorana_overlap_defect_scalar_evaluator.json"),
    Path("runs/neutrino/neutrino_absolute_attachment_scaffold.json"),
    Path("runs/neutrino/neutrino_compare_only_scale_fit.json"),
    Path("runs/neutrino/neutrino_scale_anchor.json"),
    Path("runs/neutrino/defect_weighted_mu_e_family.json"),
    Path("runs/neutrino/neutrino_attachment_normalizer_candidate_audit.json"),
    Path("runs/neutrino/neutrino_residual_amplitude_candidate_audit.json"),
    Path("runs/neutrino/neutrino_weighted_cycle_repair.json"),
    Path("runs/neutrino/neutrino_weighted_cycle_absolute_amplitude_bridge.json"),
    Path("runs/neutrino/neutrino_lambda_nu_bridge_candidate.json"),
    Path("runs/neutrino/neutrino_exact_adapter_bridge_coordinate.json"),
    Path("runs/uv/bw_fixed_local_collar_modular_transport_common_floor_scaffold.json"),
    Path("runs/flavor/overlap_edge_transport_cocycle.json"),
)

GROUP_ORDER = ["Bosons", "Leptons", "Quarks", "Hadrons"]
STATUS_COLORS = {
    "structural": "\033[96m",
    "calibration": "\033[36m",
    "secondary_quantitative": "\033[95m",
    "continuation": "\033[33m",
    "simulation_dependent": "\033[91m",
}
RESET = "\033[0m"
BOLD = "\033[1m"
DIM = "\033[2m"
GREEN = "\033[92m"
RED = "\033[91m"
YELLOW = "\033[93m"
BLUE = "\033[94m"


def _ignore(_src: str, names: list[str]) -> set[str]:
    return {name for name in names if name in EXCLUDE_NAMES}


def _run(cmd: list[str], *, cwd: Path, verbose: bool) -> None:
    completed = subprocess.run(
        cmd,
        check=False,
        cwd=cwd,
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
            cmd,
            output=completed.stdout,
            stderr=completed.stderr,
        )


def _copy_outputs(work_particles: Path, current_dir: Path) -> None:
    current_dir.mkdir(parents=True, exist_ok=True)
    outputs = [
        "RESULTS_STATUS.md",
        "results_status.json",
        "EXACT_NONHADRON_MASSES.md",
        "exact_nonhadron_masses.json",
        "EXACT_FITS_ONLY.md",
        "exact_fits_only.json",
        "particle_mass_derivation_graph.svg",
        "runs/status/status_table_forward_current.json",
        "runs/status/exact_fits_only_current.json",
        "runs/status/exact_nonhadron_masses_current.json",
        "runs/neutrino/neutrino_compare_only_scale_fit.json",
        "runs/neutrino/neutrino_two_parameter_exact_adapter.json",
        "runs/neutrino/neutrino_exact_adapter_bridge_coordinate.json",
        "runs/neutrino/neutrino_bridge_rigidity_theorem.json",
        "runs/neutrino/neutrino_absolute_attachment_theorem.json",
        "runs/uv/bw_internalization_scaffold.json",
    ]
    for rel in outputs:
        src = work_particles / rel
        dst = current_dir / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


def _mount_ancillary_particle_code(runtime_root: Path) -> None:
    src = WORKSPACE_ROOT / "arXiv" / "RC1" / "ancillary" / "code" / "particles"
    dst = runtime_root / "work" / "arXiv" / "RC1" / "ancillary" / "code" / "particles"
    if not src.exists():
        raise FileNotFoundError(f"missing ancillary particle code tree: {src}")
    dst.parent.mkdir(parents=True, exist_ok=True)
    try:
        dst.symlink_to(src, target_is_directory=True)
    except FileExistsError:
        return
    except OSError:
        shutil.copytree(src, dst, dirs_exist_ok=True)


def _seed_canonical_surface_artifacts(work_particles: Path) -> None:
    for rel in RUNTIME_SURFACED_ARTIFACTS:
        src = CODE_ROOT / rel
        dst = work_particles / rel
        if not src.exists():
            raise FileNotFoundError(f"missing canonical surfaced artifact: {src}")
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


def _read_status_markdown(current_dir: Path) -> str:
    return (current_dir / "RESULTS_STATUS.md").read_text(encoding="utf-8").rstrip()


def _read_status_json(current_dir: Path) -> dict[str, Any]:
    return json.loads((current_dir / "results_status.json").read_text(encoding="utf-8"))


def _read_optional_neutrino_fit(current_dir: Path) -> dict[str, Any] | None:
    path = current_dir / NEUTRINO_COMPARE_ONLY_FIT
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def _use_color(mode: str) -> bool:
    if mode == "always":
        return True
    if mode == "never":
        return False
    return sys.stdout.isatty() and os.environ.get("TERM", "") != "dumb"


def _style(text: str, code: str, *, enabled: bool) -> str:
    if not enabled:
        return text
    return f"{code}{text}{RESET}"


def _status_cell(status: str, *, enabled: bool) -> str:
    color = STATUS_COLORS.get(status, "")
    label = status.replace("_", " ")
    return _style(label, color, enabled=enabled) if color else label


def _flag_cell(value: bool, *, enabled: bool) -> str:
    if value:
        return _style("yes", GREEN, enabled=enabled)
    return _style("no", RED, enabled=enabled)


def _terminal_width() -> int:
    return shutil.get_terminal_size((120, 40)).columns


def _truncate(text: str, width: int) -> str:
    if len(text) <= width:
        return text
    if width <= 1:
        return text[:width]
    return text[: width - 1] + "…"


def _render_box_table(headers: list[str], rows: list[list[str]]) -> str:
    widths = [len(header) for header in headers]
    for row in rows:
        for index, cell in enumerate(row):
            widths[index] = max(widths[index], len(cell))

    def hline(left: str, mid: str, right: str) -> str:
        return left + mid.join("─" * (width + 2) for width in widths) + right

    def render_row(values: list[str]) -> str:
        return "│ " + " │ ".join(value.ljust(widths[index]) for index, value in enumerate(values)) + " │"

    parts = [
        hline("┌", "┬", "┐"),
        render_row(headers),
        hline("├", "┼", "┤"),
    ]
    parts.extend(render_row(row) for row in rows)
    parts.append(hline("└", "┴", "┘"))
    return "\n".join(parts)


def _render_terminal_report(payload: dict[str, Any], *, color: bool) -> str:
    width = _terminal_width()
    lines: list[str] = []

    title = "OPH Particle Output Table"
    lines.append(_style(title, BOLD + BLUE, enabled=color))
    lines.append(_style("Current disposable-runtime surface vs PDG/reference values", DIM, enabled=color))
    lines.append("")

    inputs = payload["inputs"]
    lines.append(
        "Inputs: "
        f"P={inputs['P']} | log_dim_H={inputs['log_dim_H']:.0e} | loops={inputs['loops']} | "
        f"with_hadrons={inputs['with_hadrons']} | hadron_profile={inputs['hadron_profile']}"
    )
    surface = payload["surface_state"]
    lines.append(
        "Surface: "
        f"{surface['public_surface_kind']} | policy={surface['surface_policy']}"
    )
    candidates = surface["active_local_public_candidates"]
    lines.append(
        "Candidates: "
        f"D10={_flag_cell(candidates['d10_mass_pair'], enabled=color)} | "
        f"D11={_flag_cell(candidates['d11_forward_seed'], enabled=color)} | "
        f"charged={_flag_cell(candidates['charged_local_candidate'], enabled=color)} | "
        f"neutrinos={_flag_cell(candidates['neutrino_local_candidate'], enabled=color)} | "
        f"quarks={_flag_cell(candidates['quark_forward_candidate'], enabled=color)} | "
        f"hadrons={_flag_cell(candidates['hadrons_enabled'], enabled=color)}"
    )
    lines.append(
        _style(
            "Hadrons stay hidden here because that lane is execution-bound and deferred from the active exact-spectrum program.",
            DIM,
            enabled=color,
        )
    )
    premise_boundaries = payload.get("premise_boundaries", {})
    uv_boundary = premise_boundaries.get("uv_bw_internalization")
    if uv_boundary:
        lines.append("")
        lines.append(_style("Premise Boundaries", BOLD, enabled=color))
        lines.append(
            "UV/BW: "
            f"{uv_boundary['status']} | next={uv_boundary['remaining_object']} | follow_on={uv_boundary['follow_on_object']}"
        )
        lines.append(_style(uv_boundary["reason_current_corpus_fails"], DIM, enabled=color))
        if uv_boundary.get("remaining_objects"):
            lines.append(
                _style(
                    f"Split objects: {uv_boundary['remaining_objects'][0]} + {uv_boundary['remaining_objects'][1]}",
                    DIM,
                    enabled=color,
                )
            )
        if uv_boundary.get("candidate_extension_route"):
            lines.append(
                "Candidate extension: "
                f"{uv_boundary['candidate_extension_status']} | target={uv_boundary['candidate_extension_target']}"
            )
            lines.append(_style(uv_boundary["candidate_extension_route"], DIM, enabled=color))

    grouped: dict[str, list[dict[str, Any]]] = {group: [] for group in GROUP_ORDER}
    for row in payload["rows"]:
        grouped.setdefault(row["group"], []).append(row)

    headers = ["Particle", "Status", "OPH (GeV)", "Reference", "Delta"]
    max_ref = max(20, min(34, width - 72))
    max_delta = max(14, min(24, width - 96))
    for group in GROUP_ORDER:
        rows = grouped.get(group, [])
        if not rows:
            continue
        lines.append("")
        lines.append(_style(group, BOLD, enabled=color))
        table_rows: list[list[str]] = []
        for row in rows:
            table_rows.append(
                [
                    row["particle"],
                    _status_cell(row["status"], enabled=False),
                    row["prediction_display_gev"],
                    _truncate(row["reference_display"], max_ref),
                    _truncate(row["delta_display"], max_delta),
                ]
            )
        table = _render_box_table(headers, table_rows)
        if color:
            colored_rows = []
            for original, row in zip(rows, table_rows):
                colored_rows.append(
                    [
                        row[0],
                        _status_cell(original["status"], enabled=True),
                        row[2],
                        row[3],
                        row[4],
                    ]
                )
            # Re-render with ANSI status cells while keeping column widths based on raw text.
            raw_widths = [len(header) for header in headers]
            for row in table_rows:
                for index, cell in enumerate(row):
                    raw_widths[index] = max(raw_widths[index], len(cell))

            def hline(left: str, mid: str, right: str) -> str:
                return left + mid.join("─" * (w + 2) for w in raw_widths) + right

            def render_row(values: list[str], raw_values: list[str]) -> str:
                padded = []
                for index, value in enumerate(values):
                    pad = raw_widths[index] - len(raw_values[index])
                    padded.append(value + (" " * pad))
                return "│ " + " │ ".join(padded) + " │"

            parts = [
                hline("┌", "┬", "┐"),
                render_row(headers, headers),
                hline("├", "┼", "┤"),
            ]
            for color_row, raw_row in zip(colored_rows, table_rows):
                parts.append(render_row(color_row, raw_row))
            parts.append(hline("└", "┴", "┘"))
            table = "\n".join(parts)
        lines.append(table)

    generated = payload["generated_utc"]
    lines.append("")
    lines.append(_style(f"Generated: {generated}", DIM, enabled=color))
    return "\n".join(lines).rstrip()


def _render_neutrino_fit_section(fit_payload: dict[str, Any], *, color: bool) -> str:
    lines: list[str] = []
    lines.append(_style("Neutrino Compare-Only Fit", BOLD, enabled=color))
    mismatch = fit_payload["central_ratio_mismatch"]
    lines.append(
        "Single-scale exact central match: "
        + _style("no", RED, enabled=color)
        + f" | fixed ratio={mismatch['predicted_ratio_21_over_32']:.7f} "
        + f"vs PDG central ratio={mismatch['reference_ratio_21_over_32']:.7f}"
    )
    lines.append(
        _style(
            f"Relative ratio mismatch: {100.0 * mismatch['relative_difference']:.3f}%",
            DIM,
            enabled=color,
        )
    )
    weighted = fit_payload["fits"]["weighted_least_squares"]
    headers = ["Fit", "lambda_nu", "m1 eV", "m2 eV", "m3 eV", "Δ21 (eV²)", "Δ32 (eV²)", "σ21", "σ32"]
    rows = []
    for key in ["solar_only", "atmospheric_only", "weighted_least_squares"]:
        fit = fit_payload["fits"][key]
        rows.append(
            [
                key.replace("_", " "),
                f"{fit['lambda_nu']:.6f}",
                f"{fit['masses_eV'][0]:.6f}",
                f"{fit['masses_eV'][1]:.6f}",
                f"{fit['masses_eV'][2]:.6f}",
                f"{fit['delta_m_sq_eV2']['21']:.6e}",
                f"{fit['delta_m_sq_eV2']['32']:.6e}",
                f"{fit['residual_sigma']['21']:+.2f}",
                f"{fit['residual_sigma']['32']:+.2f}",
            ]
        )
    lines.append(_render_box_table(headers, rows))
    lines.append(
        _style(
            "This section is compare-only. It tests whether the repaired scale-free branch can be matched by one fitted positive scale, without promoting lambda_nu to theorem-grade.",
            DIM,
            enabled=color,
        )
    )
    return "\n".join(lines)


def build_runtime(runtime_root: Path, *, with_hadrons: bool, verbose: bool) -> Path:
    work_root = runtime_root / "work"
    work_code = work_root / "code"
    work_particles = work_code / "particles"
    current_dir = runtime_root / "current"

    if work_root.exists():
        shutil.rmtree(work_root)
    work_code.mkdir(parents=True, exist_ok=True)
    shutil.copytree(CODE_ROOT, work_particles, ignore=_ignore, dirs_exist_ok=True)
    _mount_ancillary_particle_code(runtime_root)

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
        ["python3", "particles/calibration/derive_d10_ew_target_free_repair_value_law.py"],
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
        _run(step, cwd=work_code, verbose=verbose)

    _seed_canonical_surface_artifacts(work_particles)
    _run(["python3", "particles/uv/derive_bw_internalization_scaffold.py"], cwd=work_code, verbose=verbose)
    _run(["python3", "particles/neutrino/derive_neutrino_compare_only_scale_fit.py"], cwd=work_code, verbose=verbose)
    _run(["python3", "particles/neutrino/derive_neutrino_two_parameter_exact_adapter.py"], cwd=work_code, verbose=verbose)
    _run(["python3", "particles/neutrino/derive_neutrino_bridge_correction_candidate_audit.py"], cwd=work_code, verbose=verbose)
    _run(["python3", "particles/neutrino/derive_neutrino_attachment_irreducibility.py"], cwd=work_code, verbose=verbose)
    _run(["python3", "particles/neutrino/derive_neutrino_bridge_rigidity_theorem.py"], cwd=work_code, verbose=verbose)
    _run(["python3", "particles/neutrino/derive_neutrino_absolute_attachment_theorem.py"], cwd=work_code, verbose=verbose)
    _run(["python3", "particles/neutrino/derive_neutrino_weighted_cycle_shared_basis_representation.py"], cwd=work_code, verbose=verbose)
    _run(["python3", "particles/neutrino/derive_neutrino_physical_majorana_phase_theorem.py"], cwd=work_code, verbose=verbose)
    _run(["python3", "particles/neutrino/derive_exact_blocking_items.py"], cwd=work_code, verbose=verbose)
    _run(["python3", "particles/neutrino/export_forward_neutrino_closure_bundle.py"], cwd=work_code, verbose=verbose)
    _run(
        ["python3", "particles/leptons/derive_lepton_current_family_affine_anchor_theorem.py"],
        cwd=work_code,
        verbose=verbose,
    )
    _run(
        ["python3", "particles/flavor/derive_quark_current_family_selected_sheet_closure.py"],
        cwd=work_code,
        verbose=verbose,
    )

    status_cmd = ["python3", "particles/scripts/build_results_status_table.py"]
    exact_nonhadron_cmd = ["python3", "particles/scripts/build_exact_nonhadron_mass_bundle.py"]
    exact_fit_cmd = ["python3", "particles/scripts/build_exact_fit_surface.py"]
    svg_cmd = ["python3", "particles/scripts/generate_mass_derivation_svg.py"]
    if with_hadrons:
        status_cmd.append("--with-hadrons")
        svg_cmd.append("--with-hadrons")
    _run(status_cmd, cwd=work_code, verbose=verbose)
    _run(exact_nonhadron_cmd, cwd=work_code, verbose=verbose)
    _run(exact_fit_cmd, cwd=work_code, verbose=verbose)
    _run(svg_cmd, cwd=work_code, verbose=verbose)
    _copy_outputs(work_particles, current_dir)
    return current_dir


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the current particle status surface in a disposable runtime workspace.")
    parser.add_argument("--runtime-root", default=str(DEFAULT_RUNTIME_ROOT))
    parser.add_argument("--with-hadrons", action="store_true")
    parser.add_argument(
        "--format",
        choices=["terminal", "markdown", "json"],
        default="terminal",
        help="Output format for the rendered surface.",
    )
    parser.add_argument(
        "--color",
        choices=["auto", "always", "never"],
        default="auto",
        help="ANSI color mode for terminal output.",
    )
    parser.add_argument(
        "--no-print-table",
        action="store_true",
        help="Build the runtime surface without printing the rendered output.",
    )
    parser.add_argument(
        "--show-paths",
        action="store_true",
        help="Print the runtime work tree and output directory after the build.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Stream individual builder output instead of running quietly.",
    )
    args = parser.parse_args()

    runtime_root = Path(args.runtime_root).resolve()
    current_dir = build_runtime(runtime_root, with_hadrons=args.with_hadrons, verbose=args.verbose)
    if not args.no_print_table:
        if args.format == "markdown":
            print(_read_status_markdown(current_dir))
        elif args.format == "json":
            print(json.dumps(_read_status_json(current_dir), indent=2))
        else:
            color = _use_color(args.color)
            payload = _read_status_json(current_dir)
            print(_render_terminal_report(payload, color=color))
            neutrino_fit = _read_optional_neutrino_fit(current_dir)
            if neutrino_fit is not None:
                print()
                print(_render_neutrino_fit_section(neutrino_fit, color=color))
    if args.show_paths:
        if not args.no_print_table:
            print()
        print(f"runtime work tree: {runtime_root / 'work' / 'code' / 'particles'}")
        print(f"current outputs: {runtime_root / 'current'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
