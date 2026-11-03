#!/usr/bin/env python3
"""Score a frozen scale-free neutrino candidate on NuFIT 6.1 profiles.

The NuFIT release contains marginalized profile surfaces, not a public
six-dimensional likelihood.  Overlapping profiles must not be added.  This
scorer reports each selected surface separately. Under its declared bilinear
interpolation, certified component lower bounds are combined by maximum;
no interpolation-error bound for the unavailable full likelihood is claimed.

The absolute oscillation scale is not promoted in the OPH candidate.  The
``DMS/DMA`` surface is therefore minimized along the source-predicted
``Delta m21^2 / Delta m32^2`` curve, with the common scale profiled out.
"""

from __future__ import annotations

import argparse
import bisect
import contextlib
from decimal import Context, Decimal, ROUND_CEILING, ROUND_FLOOR
from fractions import Fraction
import hashlib
import heapq
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


# Receipts in this package stamp generated_utc into their own bytes, so a file
# digest of a receipt is a byte identity, not a content identity: regenerating
# an unchanged receipt moves it. Content digests drop that one key before
# canonicalising, so provenance survives a re-stamp, while the published JSON
# keeps the timestamp for human use.
TIMESTAMP_KEY = "generated_utc"


def _content_sha256(path: pathlib.Path) -> str:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(payload, dict):
        payload.pop(TIMESTAMP_KEY, None)
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


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
    """Sample the nuisance curve for diagnostics, never for a lower bound."""
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


def _log10_interval(value: Fraction) -> tuple[Fraction, Fraction]:
    """Enclose log10 of an exact rational using correctly rounded Decimal.log10.

    Directed division encloses the argument; log10 is monotone and correctly
    rounded to nearest (Python's decimal contract). One adjacent decimal on
    each side therefore encloses each endpoint, independently of float libm.
    See https://docs.python.org/3/library/decimal.html#decimal.Decimal.log10.
    """
    if value <= 0:
        raise ValueError("logarithm argument must be positive")
    down, up, near = Context(prec=70, rounding=ROUND_FLOOR), Context(prec=70, rounding=ROUND_CEILING), Context(prec=70)
    numerator, denominator = Decimal(value.numerator), Decimal(value.denominator)
    lo, hi = down.divide(numerator, denominator), up.divide(numerator, denominator)
    # log10(1)=0 exactly. Decimal's next neighbour of zero is a subnormal
    # with a million-place exponent under its default context, unnecessary
    # here and prohibitively expensive to turn into a Fraction.
    return (Fraction(0) if lo == 1 else Fraction(near.log10(lo).next_minus(near)),
            Fraction(0) if hi == 1 else Fraction(near.log10(hi).next_plus(near)))


def _outward_float(value: Fraction, upper: bool) -> float:
    result = float(value)
    if (upper and Fraction(result) < value) or (not upper and Fraction(result) > value):
        result = math.nextafter(result, math.inf if upper else -math.inf)
    return result


def _bilinear_point_bracket(grid: RectilinearGrid, x: float, y: float) -> dict[str, float]:
    """Outward float readout of the exact interpolant at binary64 inputs."""
    if not math.isfinite(x) or not math.isfinite(y):
        raise ValueError("interpolation coordinates must be finite")
    if not (grid.x_axis[0] <= x <= grid.x_axis[-1] and grid.y_axis[0] <= y <= grid.y_axis[-1]):
        raise ValueError("certified interpolation cannot extrapolate")
    x0, x1, _ = grid._bracket(grid.x_axis, x, "x")
    y0, y1, _ = grid._bracket(grid.y_axis, y, "y")
    tx = (Fraction(x) - Fraction(x0)) / (Fraction(x1) - Fraction(x0))
    ty = (Fraction(y) - Fraction(y0)) / (Fraction(y1) - Fraction(y0))
    value = ((1-tx)*(1-ty)*Fraction(grid.values[x0, y0])
             + tx*(1-ty)*Fraction(grid.values[x1, y0])
             + (1-tx)*ty*Fraction(grid.values[x0, y1])
             + tx*ty*Fraction(grid.values[x1, y1]))
    return {"lower_bound": _outward_float(value, False), "upper_bound": _outward_float(value, True)}


def _ratio_profile_bracket(
    grid: RectilinearGrid,
    ratio_dm21_over_dm32: float,
    tolerance: Fraction = Fraction(1, 100_000_000),
    max_subdivisions: int = 20_000,
) -> dict[str, Any]:
    """Certified min bracket for the declared bilinear table, not its likelihood.

    Treat parsed binary64 coordinates/ordinates and the supplied ratio as exact
    rationals. The curve is x=log10(r*y/(1000*(1+r))). Every y interval maps
    into a certified x interval. On each intersected table cell a bilinear
    function has its extrema at rectangle corners, evaluated here exactly.
    Such rectangles cover the complete curve, including crossings between
    samples. Their minima are LOWER bounds; an enclosed value at a feasible
    midpoint is an UPPER bound. Best-first dyadic subdivision closes the gap.
    No differentiability bound for the experimental likelihood is asserted.
    """
    if (isinstance(ratio_dm21_over_dm32, bool)
            or not isinstance(ratio_dm21_over_dm32, (int, float))
            or not math.isfinite(ratio_dm21_over_dm32) or ratio_dm21_over_dm32 <= 0):
        raise ValueError("the normal-ordering mass-squared ratio must be finite and positive")
    if not isinstance(tolerance, Fraction) or tolerance <= 0:
        raise ValueError("the enclosure tolerance must be a positive Fraction")
    if type(max_subdivisions) is not int or max_subdivisions < 0:
        raise ValueError("max_subdivisions must be a nonnegative integer")
    if len(grid.x_axis) < 2 or len(grid.y_axis) < 2:
        raise ValueError("profiling needs at least two coordinates on each axis")
    all_values = (*grid.x_axis, *grid.y_axis, *grid.values.values())
    if any(isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v)
           for v in all_values):
        raise ValueError("profile grid coordinates and values must be finite real numbers")
    xs, ys = tuple(map(Fraction, grid.x_axis)), tuple(map(Fraction, grid.y_axis))
    if ys[0] <= 0 or any(b <= a for axis in (xs, ys) for a, b in zip(axis, axis[1:])):
        raise ValueError("profile axes must increase, with positive atmospheric coordinates")
    r = Fraction(ratio_dm21_over_dm32)
    factor = r / (1000 * (1 + r))
    table = {(Fraction(x), Fraction(y)): Fraction(z) for (x, y), z in grid.values.items()}

    def rectangle(j: int, lo: Fraction, hi: Fraction) -> tuple[Fraction, Fraction] | None:
        # The logarithm bounds are rounded outwards even for exact powers of 10.
        xlo, xhi = _log10_interval(factor * lo)[0], _log10_interval(factor * hi)[1]
        xlo, xhi = max(xs[0], xlo), min(xs[-1], xhi)
        if xlo > xhi:
            return None
        extrema = []
        first = max(0, bisect.bisect_left(xs, xlo) - 1)
        last = min(len(xs) - 2, bisect.bisect_right(xs, xhi) - 1)
        y0, y1 = ys[j], ys[j + 1]
        for i in range(first, last + 1):
            x0, x1 = xs[i], xs[i + 1]
            for x in (max(x0, xlo), min(x1, xhi)):
                for y in (lo, hi):
                    tx, ty = (x - x0) / (x1 - x0), (y - y0) / (y1 - y0)
                    extrema.append((1-tx)*(1-ty)*table[x0, y0] + tx*(1-ty)*table[x1, y0]
                                   + (1-tx)*ty*table[x0, y1] + tx*ty*table[x1, y1])
        return min(extrema), max(extrema)

    heap: list[tuple[Fraction, int, int, Fraction, Fraction]] = []
    serial, upper = 0, None

    def insert(j: int, lo: Fraction, hi: Fraction) -> None:
        nonlocal serial, upper
        box = rectangle(j, lo, hi)
        if box is None:
            return
        heapq.heappush(heap, (box[0], serial, j, lo, hi))
        serial += 1
        mid = (lo + hi) / 2
        xlo, xhi = _log10_interval(factor * mid)
        # An upper bound needs a genuinely feasible point, not a clipped one.
        if xs[0] <= xlo and xhi <= xs[-1]:
            point = rectangle(j, mid, mid)
            assert point is not None
            upper = point[1] if upper is None else min(upper, point[1])

    for j, (lo, hi) in enumerate(zip(ys, ys[1:])):
        insert(j, lo, hi)
    if not heap:
        raise ValueError("the source ratio curve does not cross the NuFIT grid")
    subdivisions = 0
    while upper is None or upper - heap[0][0] > tolerance:
        if subdivisions >= max_subdivisions:
            raise ValueError("profile enclosure did not reach its declared tolerance")
        _, _, j, lo, hi = heapq.heappop(heap)
        mid = (lo + hi) / 2
        insert(j, lo, mid)
        insert(j, mid, hi)
        subdivisions += 1
        if not heap:
            raise ValueError("the source ratio curve has no certified feasible interval")
    lower = heap[0][0]
    assert upper is not None and lower <= upper
    exported_lower, exported_upper = _outward_float(lower, False), _outward_float(upper, True)
    return {
        "lower_bound": exported_lower,
        "upper_bound": exported_upper,
        "absolute_tolerance": str(tolerance),
        "tolerance_scope": "exact_rational_interval_before_outward_binary64_export",
        "exported_width_upper": _outward_float(Fraction(exported_upper) - Fraction(exported_lower), True),
        "subdivisions": subdivisions,
        "covering_intervals": len(heap),
        "method": "exact_rational_bilinear_rectangles_and_outward_Decimal_log10",
        "scope": "minimum_of_declared_bilinear_interpolant_on_exact_binary64_input_ratio_curve",
        "experimental_likelihood_interpolation_error_bounded": False,
    }


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
    t13_t12["certified_table_value_interval"] = _bilinear_point_bracket(
        grids["T13/T12"], coordinates["sin2_theta13"], coordinates["sin2_theta12"]
    )
    t23_dcp["certified_table_value_interval"] = _bilinear_point_bracket(
        grids["T23/DCP"], coordinates["sin2_theta23"], coordinates["delta_cp_deg_wrapped"]
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
    ratio_profile["sampled_value_role"] = "uncertified_sampled_minimum_diagnostic_not_a_lower_bound"
    ratio_profile["certified_profile_bracket"] = _ratio_profile_bracket(
        grids["DMS/DMA"], coordinates["ratio_dm21_over_dm32"]
    )
    component_values = {
        "T13/T12": t13_t12["certified_table_value_interval"]["lower_bound"],
        "T23/DCP": t23_dcp["certified_table_value_interval"]["lower_bound"],
        "DMS/DMA_ratio_profile": ratio_profile["certified_profile_bracket"]["lower_bound"],
    }
    limiting_projection = max(component_values, key=component_values.__getitem__)
    lower_bound = component_values[limiting_projection]
    maximum_profile_upper = max(
        t13_t12["certified_table_value_interval"]["upper_bound"],
        t23_dcp["certified_table_value_interval"]["upper_bound"],
        ratio_profile["certified_profile_bracket"]["upper_bound"],
    )
    if lower_bound > THREE_SIGMA_TWO_DOF_DELTA_CHI2:
        comparison_status, passes = "REJECTED", False
    elif maximum_profile_upper <= THREE_SIGMA_TWO_DOF_DELTA_CHI2:
        comparison_status, passes = "COMPATIBLE_PROFILES", True
    else:
        comparison_status, passes = "UNRESOLVED_INTERVAL", None
    return {
        "source": receipt,
        "profiles": {
            "T13/T12": t13_t12,
            "T23/DCP": t23_dcp,
            "DMS/DMA_ratio_profile": ratio_profile,
        },
        "joint_fixed_candidate_delta_chi2_lower_bound": lower_bound,
        "lower_bound_reason": (
            "Under the declared bilinear table interpolation, each angular surface profiles over "
            "undisplayed parameters. The ratio curve contributes its certified lower bound, never "
            "its sampled minimum. The maximum, not the sum, is a lower bound for these profile models; "
            "interpolation error relative to the unavailable experimental likelihood is not bounded."
        ),
        "limiting_projection": limiting_projection,
        "maximum_component_profile_upper_bound": maximum_profile_upper,
        "upper_bound_scope": "maximum_of_selected_profiles_not_an_upper_bound_on_the_full_joint_likelihood",
        "published_3sigma_2d_comparison_status": comparison_status,
        "passes_published_3sigma_2d_compatibility": passes,
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


def _prior_score_projection(
    prior: dict[str, Any], coordinates: dict[str, float], candidate_sha256: str,
    manifest_sha256: str, kernel_sha256: str, candidate_content_sha256: str | None = None,
) -> dict[str, Any]:
    """Preserve a prior comparison's inputs without inventing byte custody."""
    if prior.get("artifact") != "oph_neutrino_nufit61_retrospective_profile_score":
        raise ValueError("prior artifact is not a neutrino profile score")
    old_coordinates = prior["candidate"]["coordinates"]
    if (old_coordinates.keys() != coordinates.keys()
            or any(type(old_coordinates[k]) is not type(v) or old_coordinates[k] != v
                   for k, v in coordinates.items())):
        raise ValueError("prior scored coordinate projection changed")
    if prior["nufit_release"]["source_manifest_sha256"] != manifest_sha256:
        raise ValueError("prior comparison source manifest changed")
    if prior["source_boundary"]["family_transport_kernel_sha256"] != kernel_sha256:
        raise ValueError("prior family-kernel source changed")
    old_candidate = prior["candidate"]["sha256"]
    old_content = prior["candidate"].get("content_sha256")
    return {
        "original_candidate_sha256": old_candidate,
        "replayed_candidate_sha256": candidate_sha256,
        "candidate_byte_identity": old_candidate == candidate_sha256,
        "original_candidate_content_sha256": old_content,
        "replayed_candidate_content_sha256": candidate_content_sha256,
        "candidate_content_identity": (
            old_content is not None and candidate_content_sha256 is not None
            and old_content == candidate_content_sha256
        ),
        "scored_coordinate_projection_exactly_equal": True,
        "comparison_source_manifest_unchanged": True,
        "family_kernel_unchanged": True,
        "scope": "same_scored_coordinates_and_comparison_sources_not_a_claim_of_candidate_byte_identity",
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
    parser.add_argument("--prior-score-git-object", default="",
                        help="Optional immutable full-commit:path of the earlier comparison receipt")
    parser.add_argument("--output", default=str(DEFAULT_OUTPUT))
    args = parser.parse_args()

    candidate_path = pathlib.Path(args.candidate)
    kernel_path = pathlib.Path(args.family_kernel)
    manifest_path = pathlib.Path(args.source_manifest)
    candidate = _load_json(candidate_path)
    kernel = _load_json(kernel_path)
    manifest = _load_json(manifest_path)
    coordinates = _candidate_coordinates(candidate)
    candidate_content_sha256 = _content_sha256(candidate_path)

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
        result["published_3sigma_2d_comparison_status"] == "REJECTED" for result in scores.values()
    )
    boundary = _source_boundary(kernel)
    payload = {
        "artifact": "oph_neutrino_nufit61_retrospective_profile_score",
        "generated_utc": args.generated_utc or _timestamp(),
        "scorer": {"filename": pathlib.Path(__file__).name, "sha256": _sha256(pathlib.Path(__file__))},
        "candidate": {
            "artifact": candidate.get("artifact"),
            "filename": candidate_path.name,
            "sha256": _sha256(candidate_path),
            "content_sha256": candidate_content_sha256,
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

    if args.prior_score_git_object:
        if not re.fullmatch(r"[0-9a-f]{40}:[A-Za-z0-9_./-]+", args.prior_score_git_object):
            raise ValueError("prior score must name a full commit and repository path")
        raw = subprocess.check_output(
            ["git", "cat-file", "blob", args.prior_score_git_object], cwd=ROOT
        )
        prior = _prior_score_projection(
            json.loads(raw), coordinates, _sha256(candidate_path),
            _sha256(manifest_path), _sha256(kernel_path), candidate_content_sha256
        )
        payload["prior_comparison_provenance"] = {
            "git_object": args.prior_score_git_object,
            "receipt_sha256": hashlib.sha256(raw).hexdigest(),
            **prior,
        }

    output_path = pathlib.Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_bytes((json.dumps(payload, indent=2, sort_keys=True) + "\n").encode("utf-8"))
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
