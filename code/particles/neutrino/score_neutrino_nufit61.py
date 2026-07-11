#!/usr/bin/env python3
"""Score a frozen scale-free neutrino candidate on NuFIT 6.1 profiles.

The NuFIT release contains marginalized profile surfaces, not a public
six-dimensional likelihood.  Overlapping profiles must not be added.  This
scorer reports each selected surface separately and uses their maximum as a
lower bound on the unavailable full fixed-candidate Delta chi-square.

The absolute oscillation scale is not promoted in the OPH candidate.  The
``DMS/DMA`` surface is therefore minimized along the source-predicted
``Delta m21^2 / Delta m32^2`` curve, with the common scale profiled out.
"""

from __future__ import annotations

import argparse
import bisect
import contextlib
import hashlib
import json
import math
import pathlib
import re
import shutil
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Any, Iterator, TextIO


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_CANDIDATE = ROOT / "particles" / "runs" / "neutrino" / "neutrino_weighted_cycle_repair.json"
DEFAULT_KERNEL = ROOT / "particles" / "runs" / "flavor" / "family_transport_kernel.json"
DEFAULT_MANIFEST = ROOT / "particles" / "neutrino" / "nufit61_sources.json"
DEFAULT_OUTPUT = ROOT / "particles" / "runs" / "neutrino" / "nufit61_weighted_cycle_retrospective_score.json"

REQUIRED_SECTIONS = ("T13/T12", "T23/DCP", "DMS/DMA")
SECTION_RE = re.compile(r"^# ([A-Z0-9/]+) projection:")
THREE_SIGMA_TWO_DOF_DELTA_CHI2 = 11.829158081900795


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_json(path: pathlib.Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


@contextlib.contextmanager
def _table_lines(path: pathlib.Path) -> Iterator[TextIO]:
    if path.suffix != ".xz":
        with path.open("r", encoding="utf-8") as handle:
            yield handle
        return

    xz = shutil.which("xz")
    if xz is None:
        raise RuntimeError("reading a .xz NuFIT table requires the xz executable")
    process = subprocess.Popen(
        [xz, "-dc", str(path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
    )
    assert process.stdout is not None
    try:
        yield process.stdout
    finally:
        process.stdout.close()
        stderr = process.stderr.read() if process.stderr is not None else ""
        return_code = process.wait()
        if process.stderr is not None:
            process.stderr.close()
        if return_code != 0:
            raise RuntimeError(f"xz failed for {path.name}: {stderr.strip()}")


@dataclass(frozen=True)
class RectilinearGrid:
    """A complete two-dimensional grid with bilinear interpolation."""

    x_axis: tuple[float, ...]
    y_axis: tuple[float, ...]
    values: dict[tuple[float, float], float]

    @classmethod
    def from_rows(cls, rows: list[tuple[float, float, float]], section: str) -> "RectilinearGrid":
        if not rows:
            raise ValueError(f"NuFIT section {section} is empty")
        x_axis = tuple(sorted({row[0] for row in rows}))
        y_axis = tuple(sorted({row[1] for row in rows}))
        values: dict[tuple[float, float], float] = {}
        for x_value, y_value, delta_chi2 in rows:
            key = (x_value, y_value)
            if key in values and values[key] != delta_chi2:
                raise ValueError(f"NuFIT section {section} has conflicting duplicate point {key}")
            values[key] = delta_chi2
        expected = len(x_axis) * len(y_axis)
        if len(values) != expected:
            raise ValueError(
                f"NuFIT section {section} is not a complete grid: {len(values)} points, expected {expected}"
            )
        return cls(x_axis=x_axis, y_axis=y_axis, values=values)

    @staticmethod
    def _bracket(axis: tuple[float, ...], value: float, name: str) -> tuple[float, float, float]:
        if len(axis) < 2:
            raise ValueError(f"{name} axis needs at least two points")
        tolerance = 1.0e-12 * max(1.0, abs(value), abs(axis[0]), abs(axis[-1]))
        if value < axis[0] - tolerance or value > axis[-1] + tolerance:
            raise ValueError(f"{name}={value} lies outside [{axis[0]}, {axis[-1]}]; extrapolation is forbidden")
        value = min(max(value, axis[0]), axis[-1])
        upper = bisect.bisect_right(axis, value)
        if upper == 0:
            lower = 0
            upper = 1
        elif upper == len(axis):
            lower = len(axis) - 2
            upper = len(axis) - 1
        else:
            lower = upper - 1
        lo_value = axis[lower]
        hi_value = axis[upper]
        weight = (value - lo_value) / (hi_value - lo_value)
        return lo_value, hi_value, weight

    def interpolate(self, x_value: float, y_value: float) -> dict[str, Any]:
        x0, x1, tx = self._bracket(self.x_axis, x_value, "x")
        y0, y1, ty = self._bracket(self.y_axis, y_value, "y")
        z00 = self.values[(x0, y0)]
        z10 = self.values[(x1, y0)]
        z01 = self.values[(x0, y1)]
        z11 = self.values[(x1, y1)]
        interpolated = (
            (1.0 - tx) * (1.0 - ty) * z00
            + tx * (1.0 - ty) * z10
            + (1.0 - tx) * ty * z01
            + tx * ty * z11
        )
        return {
            "delta_chi2": float(interpolated),
            "interpolation": "bilinear_no_extrapolation",
            "point": [float(x_value), float(y_value)],
            "cell": {
                "x": [x0, x1],
                "y": [y0, y1],
                "corner_delta_chi2": [z00, z10, z01, z11],
            },
        }


def _read_grids(path: pathlib.Path) -> dict[str, RectilinearGrid]:
    rows: dict[str, list[tuple[float, float, float]]] = {name: [] for name in REQUIRED_SECTIONS}
    active: str | None = None
    with _table_lines(path) as lines:
        for raw_line in lines:
            if raw_line.startswith("# "):
                match = SECTION_RE.match(raw_line)
                active = match.group(1) if match and match.group(1) in rows else None
                continue
            if active is None or not raw_line.strip() or raw_line.startswith("#"):
                continue
            fields = raw_line.split()
            if len(fields) != 3:
                raise ValueError(f"NuFIT section {active} expected three columns, got {len(fields)}")
            rows[active].append(tuple(float(field) for field in fields))
    return {name: RectilinearGrid.from_rows(section_rows, name) for name, section_rows in rows.items()}


def _wrap_delta_degrees(value: float) -> float:
    wrapped = (value + 180.0) % 360.0 - 180.0
    return 180.0 if math.isclose(wrapped, -180.0) and value > 0.0 else wrapped


def _candidate_coordinates(candidate: dict[str, Any]) -> dict[str, float]:
    pmns = candidate["pmns_observables"]
    return {
        "sin2_theta12": math.sin(math.radians(float(pmns["theta12_deg"]))) ** 2,
        "sin2_theta13": math.sin(math.radians(float(pmns["theta13_deg"]))) ** 2,
        "sin2_theta23": math.sin(math.radians(float(pmns["theta23_deg"]))) ** 2,
        "delta_cp_deg_wrapped": _wrap_delta_degrees(float(pmns["delta_deg"])),
        "ratio_dm21_over_dm32": float(candidate["dimensionless_ratio_dm21_over_dm32"]),
    }


def _ratio_profile(
    grid: RectilinearGrid,
    ratio_dm21_over_dm32: float,
    samples: int,
) -> dict[str, Any]:
    if ratio_dm21_over_dm32 <= 0.0:
        raise ValueError("the normal-ordering mass-squared ratio must be positive")
    if samples < 1001:
        raise ValueError("ratio profiling requires at least 1001 samples")
    fraction_dm21_over_dm31 = ratio_dm21_over_dm32 / (1.0 + ratio_dm21_over_dm32)
    dma_min = max(
        grid.y_axis[0],
        (10.0 ** grid.x_axis[0]) / (fraction_dm21_over_dm31 * 1.0e-3),
    )
    dma_max = min(
        grid.y_axis[-1],
        (10.0 ** grid.x_axis[-1]) / (fraction_dm21_over_dm31 * 1.0e-3),
    )
    if dma_min >= dma_max:
        raise ValueError("the source ratio curve does not cross the NuFIT DMS/DMA grid")

    best: dict[str, Any] | None = None
    for index in range(samples):
        dma = dma_min + (dma_max - dma_min) * index / (samples - 1)
        dms_eV2 = fraction_dm21_over_dm31 * dma * 1.0e-3
        log10_dms = math.log10(dms_eV2)
        point = grid.interpolate(log10_dms, dma)
        if best is None or point["delta_chi2"] < best["delta_chi2"]:
            best = {
                **point,
                "delta_m21_sq_eV2": dms_eV2,
                "delta_m31_sq_eV2": dma * 1.0e-3,
                "delta_m32_sq_eV2": (dma * 1.0e-3) - dms_eV2,
            }
    assert best is not None
    best["curve"] = "Delta_m21^2 = [r/(1+r)] Delta_m31^2, r = Delta_m21^2/Delta_m32^2"
    best["profiled_nuisance"] = "one positive common oscillation scale"
    best["samples"] = samples
    return best


def _verify_table(path: pathlib.Path, source: dict[str, Any], source_id: str) -> dict[str, Any]:
    actual_hash = _sha256(path)
    actual_bytes = path.stat().st_size
    expected_hash = str(source["sha256"])
    expected_bytes = int(source["bytes"])
    if actual_hash != expected_hash or actual_bytes != expected_bytes:
        raise ValueError(
            f"{source_id} does not match the pinned NuFIT source: "
            f"sha256={actual_hash}, bytes={actual_bytes}"
        )
    return {
        "source_id": source_id,
        "filename": path.name,
        "sha256": actual_hash,
        "bytes": actual_bytes,
        "url": source["url"],
        "atmospheric_treatment": source.get("atmospheric_treatment"),
    }


def score_table(
    path: pathlib.Path,
    source: dict[str, Any],
    source_id: str,
    coordinates: dict[str, float],
    profile_samples: int,
) -> dict[str, Any]:
    receipt = _verify_table(path, source, source_id)
    grids = _read_grids(path)
    t13_t12 = grids["T13/T12"].interpolate(
        coordinates["sin2_theta13"], coordinates["sin2_theta12"]
    )
    t23_dcp = grids["T23/DCP"].interpolate(
        coordinates["sin2_theta23"], coordinates["delta_cp_deg_wrapped"]
    )
    ratio_profile = _ratio_profile(
        grids["DMS/DMA"], coordinates["ratio_dm21_over_dm32"], profile_samples
    )
    coarse_samples = max(1001, (profile_samples + 1) // 2)
    coarse_ratio_profile = _ratio_profile(
        grids["DMS/DMA"], coordinates["ratio_dm21_over_dm32"], coarse_samples
    )
    ratio_profile["half_resolution_samples"] = coarse_samples
    ratio_profile["half_resolution_delta_chi2_difference"] = abs(
        ratio_profile["delta_chi2"] - coarse_ratio_profile["delta_chi2"]
    )
    component_values = {
        "T13/T12": t13_t12["delta_chi2"],
        "T23/DCP": t23_dcp["delta_chi2"],
        "DMS/DMA_ratio_profile": ratio_profile["delta_chi2"],
    }
    limiting_projection = max(component_values, key=component_values.__getitem__)
    lower_bound = component_values[limiting_projection]
    return {
        "source": receipt,
        "profiles": {
            "T13/T12": t13_t12,
            "T23/DCP": t23_dcp,
            "DMS/DMA_ratio_profile": ratio_profile,
        },
        "joint_fixed_candidate_delta_chi2_lower_bound": lower_bound,
        "lower_bound_reason": (
            "Each published surface profiles over undisplayed parameters. Their maximum, not their sum, "
            "is therefore a lower bound on the unavailable full fixed-candidate Delta chi-square."
        ),
        "limiting_projection": limiting_projection,
        "passes_published_3sigma_2d_compatibility": lower_bound <= THREE_SIGMA_TWO_DOF_DELTA_CHI2,
        "wilks_two_dof_reference_p_for_limiting_projection": math.exp(-0.5 * lower_bound),
    }


def _source_boundary(kernel: dict[str, Any]) -> dict[str, Any]:
    kernel_status = str(kernel.get("status", "missing"))
    kernel_proof_status = str(kernel.get("proof_status", "missing"))
    eligible = kernel_status == "source_only_frozen" and kernel_proof_status == "closed_source_emitted"
    return {
        "family_transport_kernel_status": kernel_status,
        "family_transport_kernel_proof_status": kernel_proof_status,
        "historical_target_exposure": True,
        "historical_target_exposure_evidence": [
            "The weighted-cycle builder entered git with PDG windows and an atmospheric comparison anchor.",
            "Commit 802d93a ranked exponent-law candidates against the PDG mass-squared ratio before commit decd4a9 promoted the midpoint law.",
        ],
        "source_only_prediction_eligible": eligible,
        "prospective_evidence_eligible": False,
        "allowed_claim": "retrospective target-informed template candidate stress test",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Score the frozen weighted-cycle candidate on NuFIT 6.1 profiles.")
    parser.add_argument("--candidate", default=str(DEFAULT_CANDIDATE))
    parser.add_argument("--family-kernel", default=str(DEFAULT_KERNEL))
    parser.add_argument("--source-manifest", default=str(DEFAULT_MANIFEST))
    parser.add_argument("--tb-off-no", required=True, help="Pinned v61.release-TBoff-NO.txt.xz")
    parser.add_argument("--tb-yes-no", required=True, help="Pinned v61.release-TByes-NO.txt.xz")
    parser.add_argument("--profile-samples", type=int, default=20001)
    parser.add_argument("--generated-utc", default="")
    parser.add_argument("--output", default=str(DEFAULT_OUTPUT))
    args = parser.parse_args()

    candidate_path = pathlib.Path(args.candidate)
    kernel_path = pathlib.Path(args.family_kernel)
    manifest_path = pathlib.Path(args.source_manifest)
    candidate = _load_json(candidate_path)
    kernel = _load_json(kernel_path)
    manifest = _load_json(manifest_path)
    coordinates = _candidate_coordinates(candidate)

    scores = {
        "TBoff-NO": score_table(
            pathlib.Path(args.tb_off_no),
            manifest["files"]["TBoff-NO"],
            "TBoff-NO",
            coordinates,
            args.profile_samples,
        ),
        "TByes-NO": score_table(
            pathlib.Path(args.tb_yes_no),
            manifest["files"]["TByes-NO"],
            "TByes-NO",
            coordinates,
            args.profile_samples,
        ),
    }
    all_fail = all(
        not result["passes_published_3sigma_2d_compatibility"] for result in scores.values()
    )
    boundary = _source_boundary(kernel)
    payload = {
        "artifact": "oph_neutrino_nufit61_retrospective_profile_score",
        "generated_utc": args.generated_utc or _timestamp(),
        "candidate": {
            "artifact": candidate.get("artifact"),
            "filename": candidate_path.name,
            "sha256": _sha256(candidate_path),
            "coordinates": coordinates,
            "ordering": "normal",
        },
        "source_boundary": {
            **boundary,
            "family_transport_kernel_filename": kernel_path.name,
            "family_transport_kernel_sha256": _sha256(kernel_path),
        },
        "nufit_release": {
            "version": manifest["version"],
            "data_cutoff": manifest["data_cutoff"],
            "official_result_page": manifest["official_result_page"],
            "source_manifest_filename": manifest_path.name,
            "source_manifest_sha256": _sha256(manifest_path),
        },
        "scores": scores,
        "decision": {
            "criterion": (
                "The candidate fails the published 3-sigma two-parameter compatibility gate if the "
                "conservative lower bound exceeds Delta chi-square 11.829158 in both declared atmospheric treatments."
            ),
            "threshold_delta_chi2_2d_3sigma": THREE_SIGMA_TWO_DOF_DELTA_CHI2,
            "current_weighted_cycle_candidate_rejected_by_declared_gate": all_fail,
            "oph_core_falsified": False,
            "reason_core_not_falsified": (
                "The candidate is target-informed and descends from a template flavor kernel, so rejection applies "
                "to this weighted-cycle continuation branch rather than to finite OPH."
            ),
        },
        "statistical_scope": {
            "full_six_dimensional_likelihood_publicly_available": False,
            "profiles_summed": False,
            "wilks_p_values_are_descriptive_only": True,
            "retrospective_p_values_are_prospective_evidence": False,
            "absolute_scale_policy": "profiled on the DMS/DMA ratio curve; no OPH absolute-mass attachment is used",
        },
    }

    output_path = pathlib.Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {output_path}")
    print(
        json.dumps(
            {
                source_id: {
                    "T23/DCP": result["profiles"]["T23/DCP"]["delta_chi2"],
                    "lower_bound": result["joint_fixed_candidate_delta_chi2_lower_bound"],
                    "passes_3sigma_2d": result["passes_published_3sigma_2d_compatibility"],
                }
                for source_id, result in scores.items()
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
