#!/usr/bin/env python3
"""Conditional edge-center arithmetic and finite-geometry packet for issue #522.

The module checks the arithmetic consequences of the conditional edge-center
reserve-generator theorem (paper/tex_fragments/SCREEN_SPECTRUM_THEOREMS.tex,
Definition ``def:oph-screen-reserve-generator`` and Theorem
``thm:oph-p-over-48-tilt``) and records finite collar combinatorics. It does not
derive the reserve generator or an operational clock from those counts.

* exact integer edge-midpoint refinement combinatorics on the icosahedral
  reference carrier (collar cells vs full cells, oriented slot counts, the
  orientation-reversal involution, all exact at every tower depth);
* the certified pixel-root interval ``P`` imported from the P_derivation
  contraction certificate, whose physical endpoint remains unpromoted;
* the conditional repair-round count ``m_rep = 24`` and the separately declared
  reserve-trace branch input ``tau_q(Z6) = P/24``. The equality of their
  denominators is not a clock derivation.

Conditional arithmetic:

* if the imported full-collar density premise is
  ``-u_full'(0) = P/24``, its interval is propagated;
* if the weighted orientation-reversal coarea identity is also supplied, then
  ``theta = P/48`` and ``n_s = 1 - P/48`` as intervals, and
  ``kappa_rep = P/(48 (P - phi))`` as an interval;
* the finite presence survival family ``u_m(1 tick) = (1 - (P/24) 2^-m)^(2^m)``
  constructed from that imported density, with its arithmetic defects computed
  and bounded at every available depth;
* the finite-transition exponent ``-log u / log b`` and the ``e`` branch as
  diagnostics only.

Equal oriented-slot counts alone do not force equal reserve weights on the two
orientations. The operational clock, source-emitted weighted half identity,
simulator receipt consumption, freeze action, and sky comparison remain open.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any

from mpmath import iv, mp, mpf

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DEFAULT_OUT = HERE / "manifests" / "edge_center_clock_certificate.json"

SCHEMA = "oph.edge_center_clock_certificate.v3"
ARTIFACT = "oph_edge_center_clock_certificate"
GITHUB_ISSUE = 522
PRECISION = 60
TOWER_DEPTH = 12

P_CERTIFICATE_PATH = (
    ROOT
    / "code"
    / "P_derivation"
    / "runtime"
    / "p_interval_contraction_certificate_2026-07-14.json"
)
P_CERTIFICATE_MODE = "thomson_structured_running"

M_REP_CERTIFICATE_PATH = (
    ROOT / "code" / "particles" / "hierarchy" / "certificates" / "R_m_rep_24_certificate.json"
)

SELECTED_PRIMARY_BRANCH = "conditional_edge_center_P_over_48"

# Defect gates (all verified with mpmath.iv outward-rounded interval arithmetic).
# Identities that hold pointwise for every representative of the P enclosure
# still produce a nonzero interval width when both sides carry the enclosure
# independently (interval dependency); the width of the P enclosure is about
# 6.8e-28, so 1e-25 is a safe representation-level gate.  A thin-representative
# check at the midpoint additionally verifies the identities to 1e-50.
EXACT_DEFECT_WIDTH_BOUND = "1e-25"
THIN_DEFECT_WIDTH_BOUND = "1e-50"

# Input hygiene: no time-dimension constant and no measured-tilt-shaped target
# may enter the computation.
FORBIDDEN_INPUT_KEY_SUBSTRINGS = (
    "time",
    "second",
    "hertz",
    "planck",
    "target",
    "measured",
    "sky",
    "likelihood",
    "fit",
    "n_s",
    "tilt",
)
MEASURED_TILT_WINDOW = ("0.95", "0.98")


class CertificateError(ValueError):
    """Raised when an input violates the contract or a defect bound fails."""


# ---------------------------------------------------------------------------
# interval helpers
# ---------------------------------------------------------------------------


def _endpoints(x: Any) -> tuple[Any, Any]:
    return mpf(x.a), mpf(x.b)


def _interval_json(x: Any) -> dict[str, str]:
    lo, hi = _endpoints(x)
    return {"lo": mp.nstr(lo, 40), "hi": mp.nstr(hi, 40)}


def _width(x: Any) -> Any:
    lo, hi = _endpoints(x)
    return hi - lo


def _contains_zero(x: Any) -> bool:
    lo, hi = _endpoints(x)
    return lo <= 0 <= hi

def _contains_value(x: Any, decimal_text: str) -> bool:
    lo, hi = _endpoints(x)
    value = mpf(decimal_text)
    return lo <= value <= hi


def _disjoint(x: Any, y: Any) -> bool:
    x_lo, x_hi = _endpoints(x)
    y_lo, y_hi = _endpoints(y)
    return x_hi < y_lo or y_hi < x_lo


def _fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def _sha256_bytes(path: Path) -> str:
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


def _payload_sha256(payload: dict[str, Any]) -> str:
    body = {key: value for key, value in payload.items() if key != "artifact_sha256"}
    encoded = json.dumps(
        body,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    ).encode("utf-8")
    return "sha256:" + hashlib.sha256(encoded).hexdigest()


def verify_artifact_self_hash(payload: dict[str, Any]) -> None:
    if payload.get("artifact_sha256") != _payload_sha256(payload):
        raise CertificateError("artifact self-hash mismatch")


# ---------------------------------------------------------------------------
# antecedent import: the certified pixel-root interval
# ---------------------------------------------------------------------------


def load_certified_p_interval(path: Path = P_CERTIFICATE_PATH) -> dict[str, Any]:
    """Load the certified P enclosure; P never enters as a measured alpha."""
    if not path.exists():
        raise CertificateError(f"missing P interval certificate: {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("artifact") != "oph_p_interval_contraction_certificate":
        raise CertificateError("unexpected artifact in P interval certificate")
    mode = payload["modes"][P_CERTIFICATE_MODE]
    enclosure = mode["certified_enclosure"]["P"]
    return {
        "lo": enclosure["lo"],
        "hi": enclosure["hi"],
        "source": path.relative_to(ROOT).as_posix(),
        "source_sha256": _sha256_bytes(path),
        "mode": P_CERTIFICATE_MODE,
        "certificate_date": payload.get("date"),
        "consumer_policy": payload.get("consumer_policy", {}),
        "claim_boundary": payload.get("claim_boundary"),
        "exact_alpha_promoted": bool(payload.get("exact_alpha_promoted", False)),
    }


def load_repair_round_invariant(path: Path = M_REP_CERTIFICATE_PATH) -> dict[str, Any]:
    """Load the derived round-count invariant m_rep = 24; fail closed."""
    if not path.exists():
        raise CertificateError(f"missing repair-round certificate: {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    result = payload.get("result", {})
    if not (
        payload.get("artifact") == "R_m_rep_24_certificate"
        and payload.get("accepted") is True
        and result.get("m_rep") == 24
        and result.get("exponent_denominator") == 48
    ):
        raise CertificateError("repair-round certificate does not certify m_rep = 24")
    return {
        "artifact": path.relative_to(ROOT).as_posix(),
        "artifact_sha256": _sha256_bytes(path),
        "status": payload.get("status"),
        "m_rep": result["m_rep"],
        "exponent_denominator": result["exponent_denominator"],
        "specialized_exponent": result.get("specialized_exponent"),
        "branch_scope": payload.get("branch_scope", {}),
    }


# ---------------------------------------------------------------------------
# input hygiene: no declared step time, no measured-tilt target
# ---------------------------------------------------------------------------


def _numeric_or_none(value: Any) -> Any:
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        return mpf(value)
    if isinstance(value, str):
        try:
            return mpf(value)
        except (ValueError, TypeError):
            return None
    return None


def reject_forbidden_inputs(inputs: dict[str, Any], path: str = "inputs") -> None:
    """Fail closed on time-dimension keys and measured-tilt-shaped values.

    No external step time is admissible until an operational binding exists,
    so any declared step time is a contract violation. Any numeric input inside
    the measured spectral-tilt window is a target injection.
    """
    window_lo = mpf(MEASURED_TILT_WINDOW[0])
    window_hi = mpf(MEASURED_TILT_WINDOW[1])
    for key, value in inputs.items():
        key_lower = str(key).lower()
        for token in FORBIDDEN_INPUT_KEY_SUBSTRINGS:
            if token in key_lower:
                raise CertificateError(
                    f"forbidden input key at {path}.{key}: contains {token!r}"
                )
        if isinstance(value, dict):
            reject_forbidden_inputs(value, path=f"{path}.{key}")
            continue
        candidates = value if isinstance(value, (list, tuple)) else [value]
        for item in candidates:
            numeric = _numeric_or_none(item)
            if numeric is not None and window_lo <= numeric <= window_hi:
                raise CertificateError(
                    f"measured-tilt-shaped input value at {path}.{key}: {item!r} "
                    f"lies inside [{MEASURED_TILT_WINDOW[0]}, {MEASURED_TILT_WINDOW[1]}]"
                )


# ---------------------------------------------------------------------------
# finite collar-count records: exact icosahedral edge-midpoint tower
# ---------------------------------------------------------------------------


def collar_tower_records(depth: int = TOWER_DEPTH) -> list[dict[str, Any]]:
    """Exact integer collar records on the icosahedral refinement tower.

    Edge-midpoint (4-to-1) refinement at depth n:
    ``F_n = 20*4^n``, ``E_n = 30*4^n``, ``V_n = 2 + 10*4^n``.  A collar cell is
    a depth-n face sharing a sub-edge with the depth-0 edge skeleton.  Two
    independent exact counting routes must agree:

    * per depth-0 face: ``3*2^n - 3`` boundary-incident small triangles
      (corner triangles counted once), all 20 faces collar at depth 0;
    * per skeleton sub-edge: two incident faces each, corner triangles
      (which meet two skeleton edges) corrected once: ``2*30*2^n - 60``.
    """
    if depth < 1:
        raise CertificateError("tower depth must be at least 1")
    records: list[dict[str, Any]] = []
    for n in range(depth + 1):
        four_n = 4**n
        two_n = 2**n
        faces = 20 * four_n
        edges = 30 * four_n
        vertices = 2 + 10 * four_n
        skeleton_sub_edges = 30 * two_n
        if n == 0:
            collar_route_a = 20
            collar_route_b = 20
        else:
            collar_route_a = 20 * (3 * two_n - 3)
            collar_route_b = 2 * skeleton_sub_edges - 60
        if collar_route_a != collar_route_b:
            raise CertificateError(f"collar counting routes disagree at depth {n}")
        euler = vertices - edges + faces
        if euler != 2:
            raise CertificateError(f"Euler characteristic failed at depth {n}")
        collar = collar_route_a
        full_cells = faces - collar
        oriented_slots = 2 * skeleton_sub_edges
        forward_slots = skeleton_sub_edges
        backward_slots = skeleton_sub_edges
        records.append(
            {
                "depth": n,
                "vertices": vertices,
                "edges": edges,
                "faces": faces,
                "euler_characteristic": euler,
                "skeleton_sub_edges": skeleton_sub_edges,
                "collar_cells": collar,
                "collar_cells_route_a_per_face": collar_route_a,
                "collar_cells_route_b_per_sub_edge": collar_route_b,
                "full_cells": full_cells,
                "collar_fraction": _fraction_text(Fraction(collar, faces)),
                "oriented_collar_slots": oriented_slots,
                "forward_oriented_slots": forward_slots,
                "backward_oriented_slots": backward_slots,
                "orientation_balance_defect": forward_slots - backward_slots,
                "orientation_reversal_involution": {
                    "fixed_point_free": True,
                    "orbit_count": skeleton_sub_edges,
                    "half_collar_slot_count": forward_slots,
                    "half_over_full_slot_ratio": _fraction_text(
                        Fraction(forward_slots, oriented_slots)
                    ),
                },
                "sub_slots_per_depth0_slot": two_n,
            }
        )
    return records


def orientation_weight_countermodel(
    tower: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """Show that equal slot counts do not force a weighted half identity.

    At every depth of the actual collar tower, N forward and N backward slots
    admit both equal positive weights and unequal positive weights. The
    counting involution is unchanged, while the forward reserve share changes
    from one half to two thirds.
    """

    records = collar_tower_records() if tower is None else tower
    if not records:
        raise CertificateError("orientation-weight countermodel needs a tower")
    schedule: list[dict[str, Any]] = []
    equal_forward_weight = Fraction(1, 1)
    equal_backward_weight = Fraction(1, 1)
    biased_forward_weight = Fraction(2, 1)
    biased_backward_weight = Fraction(1, 1)
    for record in records:
        forward_slots = record["forward_oriented_slots"]
        backward_slots = record["backward_oriented_slots"]
        if (
            type(forward_slots) is not int
            or type(backward_slots) is not int
            or forward_slots <= 0
            or forward_slots != backward_slots
        ):
            raise CertificateError(
                "orientation-weight countermodel requires equal positive "
                "forward/backward tower counts"
            )
        equal_forward_total = forward_slots * equal_forward_weight
        equal_total = (
            equal_forward_total + backward_slots * equal_backward_weight
        )
        biased_forward_total = forward_slots * biased_forward_weight
        biased_total = (
            biased_forward_total + backward_slots * biased_backward_weight
        )
        schedule.append(
            {
                "depth": record["depth"],
                "forward_slot_count": forward_slots,
                "backward_slot_count": backward_slots,
                "equal_weight_forward_reserve_share": _fraction_text(
                    equal_forward_total / equal_total
                ),
                "biased_positive_weight_forward_reserve_share": _fraction_text(
                    biased_forward_total / biased_total
                ),
            }
        )

    return {
        "tower_slot_schedule": schedule,
        "tower_slot_schedule_sha256": _payload_sha256({"schedule": schedule}),
        "all_depths_match_exact_collar_tower_counts": True,
        "slot_count_ratio_at_every_depth": "1/2",
        "equal_weight_forward_reserve_share_at_every_depth": "1/2",
        "biased_positive_weight_forward_reserve_share_at_every_depth": "2/3",
        "all_weights_strictly_positive": True,
        "same_fixed_point_free_counting_involution": True,
        "conclusion": (
            "equal orientation counts do not imply an equal weighted reserve "
            "share; the weighted coarea identity requires an additional "
            "orientation-invariance premise or source receipt"
        ),
        "exhaustive_family": (
            "at a depth with N=30*2^n slots in each orientation, every "
            "rational r in (0,1) is realized by assigning weight r/N to each "
            "forward slot and (1-r)/N to each backward slot; all weights are "
            "positive and the unchanged N+N counts have forward share r"
        ),
    }


def count_only_independence_no_go(
    tower: list[dict[str, Any]],
    weight_countermodel: dict[str, Any],
) -> dict[str, Any]:
    """Prove that the exact count interface selects no density, half, or clock.

    The exported geometry is held fixed while reserve weights, a formal
    continuous-semigroup generator coefficient, and the conversion from repair
    event number to logarithmic thickness vary independently.
    """

    geometry_projection = [
        {
            "depth": row["depth"],
            "vertices": row["vertices"],
            "edges": row["edges"],
            "faces": row["faces"],
            "collar_cells": row["collar_cells"],
            "oriented_collar_slots": row["oriented_collar_slots"],
            "forward_oriented_slots": row["forward_oriented_slots"],
            "backward_oriented_slots": row["backward_oriented_slots"],
        }
        for row in tower
    ]
    geometry_sha256 = _payload_sha256({"geometry": geometry_projection})
    expected_schedule = [
        {
            "depth": row["depth"],
            "forward_slot_count": row["forward_oriented_slots"],
            "backward_slot_count": row["backward_oriented_slots"],
            "equal_weight_forward_reserve_share": "1/2",
            "biased_positive_weight_forward_reserve_share": "2/3",
        }
        for row in tower
    ]
    if weight_countermodel["tower_slot_schedule"] != expected_schedule:
        raise CertificateError(
            "orientation-weight witness does not use the exact collar tower counts"
        )
    reserve_models = [
        {
            "name": "backward_biased_positive_weights",
            "forward_share": "1/3",
            "geometry_sha256": geometry_sha256,
        },
        {
            "name": "orientation_invariant",
            "forward_share": "1/2",
            "geometry_sha256": geometry_sha256,
        },
        {
            "name": "forward_biased_positive_weights",
            "forward_share": "2/3",
            "geometry_sha256": geometry_sha256,
        },
    ]
    generator_models = [
        {
            "name": "formal_semigroup_g_1",
            "generator_coefficient": "1/24",
            "family": "u(s)=exp(-(1/24)s)",
            "geometry_sha256": geometry_sha256,
        },
        {
            "name": "formal_semigroup_g_2",
            "generator_coefficient": "1/12",
            "family": "u(s)=exp(-(1/12)s)",
            "geometry_sha256": geometry_sha256,
        },
    ]
    clock_models = [
        {
            "name": "unit_event_scale",
            "log_thickness_per_event": "1",
            "coordinate": "s=n",
            "geometry_sha256": geometry_sha256,
        },
        {
            "name": "double_event_scale",
            "log_thickness_per_event": "2",
            "coordinate": "s=2n",
            "geometry_sha256": geometry_sha256,
        },
    ]
    return {
        "count_interface_sha256": geometry_sha256,
        "same_counts_multiple_reserve_weights": reserve_models,
        "same_counts_multiple_generator_coefficients": generator_models,
        "same_counts_multiple_event_to_log_thickness_scales": clock_models,
        "orientation_weight_family": weight_countermodel["exhaustive_family"],
        "conclusion": (
            "no deterministic selector whose complete input is the exact collar "
            "count interface can uniquely determine P/24, a weighted half "
            "factor, a reserve generator, or an event-to-log-thickness clock "
            "scale"
        ),
        "exhaustive_scope": (
            "all predicates and deterministic algorithms whose complete input "
            "is the exact integer collar-count tower"
        ),
        "negative_closure_at_count_only_interface": True,
    }


def assert_counts_select_unique_reserve_share(no_go: dict[str, Any]) -> None:
    shares = {
        row["forward_share"]
        for row in no_go["same_counts_multiple_reserve_weights"]
    }
    if len(shares) != 1:
        raise CertificateError(
            "COUNT_WEIGHT_INDEPENDENCE: identical counts admit multiple reserve shares"
        )


def assert_counts_select_unique_generator(no_go: dict[str, Any]) -> None:
    coefficients = {
        row["generator_coefficient"]
        for row in no_go["same_counts_multiple_generator_coefficients"]
    }
    if len(coefficients) != 1:
        raise CertificateError(
            "COUNT_GENERATOR_INDEPENDENCE: identical counts admit multiple generators"
        )


def assert_counts_select_unique_clock_scale(no_go: dict[str, Any]) -> None:
    scales = {
        row["log_thickness_per_event"]
        for row in no_go["same_counts_multiple_event_to_log_thickness_scales"]
    }
    if len(scales) != 1:
        raise CertificateError(
            "COUNT_CLOCK_INDEPENDENCE: identical counts admit multiple clock scales"
        )


# ---------------------------------------------------------------------------
# finite survival family and defect bounds
# ---------------------------------------------------------------------------


def survival_family_records(P: Any, depth: int = TOWER_DEPTH) -> list[dict[str, Any]]:
    """Construct a finite presence-survival family from the imported density.

    At depth m the clock tick subdivides into the exact integer ``2^m``
    sub-slots. Conditional on uniform splitting of the imported reserve-presence
    density (a branch premise, not a consequence of the collar counts), the
    finite one-tick survival is ``u_m(1) = (1 - (P/24) 2^-m)^(2^m)``.

    Bounds verified at every depth (mpmath.iv, outward rounding):

    * on-grid semigroup defect ``q^3 q^5 - q^8`` contains 0, width below the
      exact-defect gate;
    * sub-step derivative defect ``(1 - q_m) 2^m - P/24`` contains 0, width
      below the exact-defect gate;
    * limit-family derivative defect ``P/24 - (1 - e^(-(P/24) h_m))/h_m`` is
      nonnegative with sup at most ``(P/24)_hi^2 2^-m``;
    * orientation-balance survival defect ``u_half_m(1)^2 - u_m(1)`` is
      nonnegative with sup at most ``(P/24)_hi^2 2^-m``;
    * refinement defect ``e^(-P/24) - u_m(1)`` is strictly positive
      (the Poisson value is unreachable at any finite regulator) with sup at
      most ``(P/24)_hi^2 2^-m``, and ``u_m(1)`` is strictly increasing in m.
    """
    # Conditional branch arithmetic. P/24 is imported from the declared
    # reserve-trace premise; it is not inferred from the collar counts below.
    eps = P / 24
    eps_half = P / 48
    _, eps_hi = _endpoints(eps)
    exact_gate = mpf(EXACT_DEFECT_WIDTH_BOUND)
    poisson_floor = iv.exp(-eps)
    records: list[dict[str, Any]] = []
    previous_one_tick = None
    for m in range(depth + 1):
        sub_slots = 2**m
        q_full = 1 - eps / sub_slots
        q_half = 1 - eps_half / sub_slots
        u_one_tick = q_full**sub_slots
        u_half_one_tick = q_half**sub_slots

        semigroup_defect = q_full**3 * q_full**5 - q_full**8
        if not _contains_zero(semigroup_defect) or _width(semigroup_defect) > exact_gate:
            raise CertificateError(f"semigroup defect gate failed at depth {m}")

        derivative_defect = (1 - q_full) * sub_slots - eps
        if not _contains_zero(derivative_defect) or _width(derivative_defect) > exact_gate:
            raise CertificateError(f"derivative defect gate failed at depth {m}")

        shrinking_bound = eps_hi * eps_hi / sub_slots

        limit_derivative_defect = eps - (1 - iv.exp(-eps / sub_slots)) * sub_slots
        limit_lo, limit_hi = _endpoints(limit_derivative_defect)
        if limit_lo < -exact_gate or limit_hi > shrinking_bound:
            raise CertificateError(f"limit derivative defect gate failed at depth {m}")

        orientation_defect = u_half_one_tick**2 - u_one_tick
        orientation_lo, orientation_hi = _endpoints(orientation_defect)
        if orientation_lo < -exact_gate or orientation_hi > shrinking_bound:
            raise CertificateError(f"orientation-balance defect gate failed at depth {m}")

        refinement_defect = poisson_floor - u_one_tick
        refinement_lo, refinement_hi = _endpoints(refinement_defect)
        if refinement_lo <= 0 or refinement_hi > shrinking_bound:
            raise CertificateError(f"refinement defect gate failed at depth {m}")

        if previous_one_tick is not None:
            step = u_one_tick - previous_one_tick
            step_lo, _ = _endpoints(step)
            if step_lo <= 0:
                raise CertificateError(f"survival monotonicity failed at depth {m}")
        previous_one_tick = u_one_tick

        records.append(
            {
                "depth": m,
                "sub_slots_per_tick": sub_slots,
                "sub_step_thickness_ticks": _fraction_text(Fraction(1, sub_slots)),
                "one_sub_step_survival": _interval_json(q_full),
                "one_tick_survival": _interval_json(u_one_tick),
                "one_tick_survival_half_collar": _interval_json(u_half_one_tick),
                "defects": {
                    "semigroup_on_grid": {
                        "value": _interval_json(semigroup_defect),
                        "contains_zero": True,
                        "width_bound": EXACT_DEFECT_WIDTH_BOUND,
                    },
                    "sub_step_derivative": {
                        "value": _interval_json(derivative_defect),
                        "contains_zero": True,
                        "width_bound": EXACT_DEFECT_WIDTH_BOUND,
                    },
                    "limit_family_derivative": {
                        "value": _interval_json(limit_derivative_defect),
                        "nonnegative": True,
                        "bound": mp.nstr(shrinking_bound, 20),
                    },
                    "orientation_balance": {
                        "value": _interval_json(orientation_defect),
                        "nonnegative": True,
                        "bound": mp.nstr(shrinking_bound, 20),
                    },
                    "refinement_to_poisson_floor": {
                        "value": _interval_json(refinement_defect),
                        "strictly_positive": True,
                        "bound": mp.nstr(shrinking_bound, 20),
                    },
                },
            }
        )
    return records


def thin_representative_exactness(P_mid: Any, depth: int = TOWER_DEPTH) -> dict[str, Any]:
    """Verify the on-grid identities at a thin midpoint representative.

    At a single representative the semigroup and sub-step derivative defects
    collapse to representation rounding, verifying that the full-enclosure
    defect widths are inherited from the P enclosure alone.
    """
    thin_gate = mpf(THIN_DEFECT_WIDTH_BOUND)
    eps = P_mid / 24
    worst_semigroup = mpf(0)
    worst_derivative = mpf(0)
    for m in (0, depth // 2, depth):
        sub_slots = 2**m
        q = 1 - eps / sub_slots
        semigroup_defect = q**3 * q**5 - q**8
        derivative_defect = (1 - q) * sub_slots - eps
        for defect, label in (
            (semigroup_defect, "semigroup"),
            (derivative_defect, "derivative"),
        ):
            if not _contains_zero(defect) or _width(defect) > thin_gate:
                raise CertificateError(
                    f"thin-representative {label} exactness failed at depth {m}"
                )
        worst_semigroup = max(worst_semigroup, _width(semigroup_defect))
        worst_derivative = max(worst_derivative, _width(derivative_defect))
    return {
        "representative": "midpoint of the certified P enclosure",
        "depths_checked": [0, depth // 2, depth],
        "worst_semigroup_defect_width": mp.nstr(worst_semigroup, 10),
        "worst_derivative_defect_width": mp.nstr(worst_derivative, 10),
        "width_bound": THIN_DEFECT_WIDTH_BOUND,
    }


# ---------------------------------------------------------------------------
# certificate builder
# ---------------------------------------------------------------------------


def build(
    *,
    injected_inputs: dict[str, Any] | None = None,
    primary_branch: str = SELECTED_PRIMARY_BRANCH,
    tower_depth: int = TOWER_DEPTH,
) -> dict[str, Any]:
    """Build the finite edge-center generator and clock certificate."""
    if primary_branch != SELECTED_PRIMARY_BRANCH:
        raise CertificateError(
            f"primary branch must be {SELECTED_PRIMARY_BRANCH!r}; "
            f"{primary_branch!r} branches are diagnostics only"
        )

    iv.dps = PRECISION
    mp.dps = PRECISION

    p_record = load_certified_p_interval()
    m_rep_record = load_repair_round_invariant()
    structure_integers = {
        "adjoint_dimension": 12,
        "orientation_doubling": 2,
        "m_rep_round_count": 24,
        "half_collar_denominator": 48,
        "icosahedron": {"vertices": 12, "edges": 30, "faces": 20},
        "refinement_scale_ratio": 2,
    }
    inputs: dict[str, Any] = {
        "P_certified_enclosure": {"lo": p_record["lo"], "hi": p_record["hi"]},
        "P_source": p_record["source"],
        "structure_integers": structure_integers,
    }
    if injected_inputs:
        inputs = {**inputs, **injected_inputs}
    reject_forbidden_inputs(inputs)

    P = iv.mpf([p_record["lo"], p_record["hi"]])
    phi = (1 + iv.sqrt(5)) / 2
    eps = P / 24
    theta = P / 48
    n_s = 1 - theta
    kappa_rep = P / (48 * (P - phi))

    tower = collar_tower_records(tower_depth)
    weight_countermodel = orientation_weight_countermodel(tower)
    count_no_go = count_only_independence_no_go(tower, weight_countermodel)
    family = survival_family_records(P, tower_depth)
    p_lo, p_hi = _endpoints(P)
    thin_exactness = thin_representative_exactness(iv.mpf((p_lo + p_hi) / 2), tower_depth)

    # Exponential-family semigroup identity (algebraic; interval-verified).
    s_probe = iv.mpf([3, 3]) / 8
    t_probe = iv.mpf([5, 5]) / 8
    exponential_semigroup_defect = iv.exp(-eps * s_probe) * iv.exp(-eps * t_probe) - iv.exp(
        -eps * (s_probe + t_probe)
    )
    if not _contains_zero(exponential_semigroup_defect):
        raise CertificateError("exponential-family semigroup identity failed")
    if _width(exponential_semigroup_defect) > mpf(EXACT_DEFECT_WIDTH_BOUND):
        raise CertificateError("exponential-family semigroup width gate failed")

    # Conditional half-collar arithmetic: once the weighted coarea premise is
    # supplied, theta computed through that premise must agree with P/48.
    theta_via_half = eps / 2
    if _disjoint(theta, theta_via_half):
        raise CertificateError("half-collar identity arithmetic inconsistency")

    # Orientation count balance is exact at every tower depth. The separate
    # countermodel above records why this does not imply weighted balance.
    if any(record["orientation_balance_defect"] != 0 for record in tower):
        raise CertificateError("orientation balance defect nonzero on the tower")
    if any(
        record["orientation_reversal_involution"]["half_over_full_slot_ratio"] != "1/2"
        for record in tower
    ):
        raise CertificateError("half-collar slot ratio is not exactly 1/2 on the tower")

    # Diagnostics (never selected).
    presence_one_step_full = 1 - eps
    presence_one_step_half = 1 - theta
    finite_transition_exponent = -iv.log(presence_one_step_half)
    if not _disjoint(finite_transition_exponent, theta):
        raise CertificateError(
            "finite-transition exponent interval overlaps the generator interval; "
            "the diagnostic separation failed"
        )
    e_branch_theta = iv.exp(iv.mpf([1, 1])) * (P - phi)
    e_branch_n_s = 1 - e_branch_theta

    # Alternative factor-1 value. It is numerically distinct, but slot counts
    # alone do not reject it; the weighted coarea premise selects the half.
    theta_wrong_orientation = eps
    if not _disjoint(theta_wrong_orientation, theta):
        raise CertificateError("wrong-orientation branch is not separated from theta")

    # The issue body quotes kappa_rep = 2.627023712627471 at the distinct
    # canonical comparison coordinate P_C = 1.630968209403959. P_C lies
    # outside this issue's selected P_fwd enclosure. Both coordinates remain
    # typed separately, and the emitted object here is the P_fwd interval.
    p_c = iv.mpf(["1.630968209403959", "1.630968209403959"])
    kappa_rep_at_p_c = p_c / (48 * (p_c - phi))
    p_c_defect = kappa_rep_at_p_c - iv.mpf(
        ["2.627023712627471", "2.627023712627471"]
    )
    p_c_lo, p_c_hi = _endpoints(p_c_defect)
    if max(abs(p_c_lo), abs(p_c_hi)) > mpf("1e-13"):
        raise CertificateError(
            "issue-body kappa_rep value is not reproduced at the P_C comparison coordinate"
        )
    if _contains_value(kappa_rep, "2.627023712627471"):
        raise CertificateError(
            "P_C kappa_rep point unexpectedly lies inside the P_fwd certified interval"
        )

    controls = _run_controls(weight_countermodel, count_no_go)

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "artifact": ARTIFACT,
        "github_issue": GITHUB_ISSUE,
        "status": "conditional_edge_center_arithmetic_with_open_source_and_clock_gates",
        "selected_primary_branch": SELECTED_PRIMARY_BRANCH,
        "theorem_source": {
            "fragment": "paper/tex_fragments/SCREEN_SPECTRUM_THEOREMS.tex",
            "definition": "def:oph-screen-reserve-generator",
            "theorem": "thm:oph-p-over-48-tilt",
        },
        "interval_backend": {
            "library": "mpmath.iv",
            "precision_decimal_digits": PRECISION,
            "rounding": "binary interval floating-point with outward enclosures",
            "exact_arithmetic": False,
            "exact_integer_geometry": True,
            "claim": (
                "numerical intervals are rigorous enclosures at the declared "
                "precision; they are not exact rational arithmetic"
            ),
        },
        "inputs": inputs,
        "input_hygiene": {
            "forbidden_key_substrings": list(FORBIDDEN_INPUT_KEY_SUBSTRINGS),
            "measured_tilt_window_rejected": list(MEASURED_TILT_WINDOW),
            "declared_step_time_present": False,
            "time_dimension_constants_present": False,
            "measured_alpha_present": False,
        },
        "antecedents": {
            "pixel_root_interval": {
                "artifact": p_record["source"],
                "artifact_sha256": p_record["source_sha256"],
                "mode": p_record["mode"],
                "certificate_date": p_record["certificate_date"],
                "consumer_policy": p_record["consumer_policy"],
                "claim_boundary": p_record["claim_boundary"],
                "exact_alpha_promoted": p_record["exact_alpha_promoted"],
                "role": (
                    "P enters as this conditional closure-map enclosure; the "
                    "artifact does not promote it to a physical endpoint"
                ),
            },
            "repair_round_invariant": {
                "artifact": m_rep_record["artifact"],
                "artifact_sha256": m_rep_record["artifact_sha256"],
                "producer": "code/particles/hierarchy/verify_issue_343_m_rep_24.py",
                "certificate_status": m_rep_record["status"],
                "branch_scope": m_rep_record["branch_scope"],
                "m_rep": m_rep_record["m_rep"],
                "exponent_denominator": m_rep_record["exponent_denominator"],
                "specialized_exponent": m_rep_record["specialized_exponent"],
                "statement": (
                    "m_rep = 2*(8+3+1) = 24 on the declared product-adjoint, "
                    "orientation-doubling, cyclic-scheduler branch"
                ),
                "clock_role": (
                    "a conditional scheduler count only; no binding to the "
                    "collar refinement coordinate or operational clock is "
                    "proved in this packet"
                ),
            },
            "reserve_trace_branch_input": {
                "artifact": "paper/screen_microphysics_and_observer_synchronization.tex",
                "statement": (
                    "tau_q(Z6) = P/24 on the shared-edge protected-reserve "
                    "branch"
                ),
                "premises": [
                    "declared shared-cut density P/4",
                    "declared Z6 class equidistribution",
                    "scalar-reserve unbiasedness where a scalar-weighted mean is consumed",
                ],
                "derived_from_collar_counts": False,
                "denominator_match_not_used_as_a_derivation": True,
                "binding_theorem_present": False,
            },
            "presence_reading": {
                "artifact": "extra/chi_nu_susceptibility_bounds.tex",
                "statement": (
                    "the finite one-step survival is the presence value 1 - P/24; "
                    "e^(-P/24) is the depth limit, unreachable at any finite regulator"
                ),
            },
        },
        "carrier_tower": {
            "carrier": "icosahedral screen cellulation, edge-midpoint refinement",
            "depths": tower_depth,
            "records": tower,
            "derives_reserve_density": False,
            "derives_weighted_half_identity": False,
            "derives_operational_clock": False,
            "boundary": (
                "the tower certifies integer incidence and orientation counts "
                "only; reserve weights and time normalization require separate "
                "source data"
            ),
        },
        "count_only_independence_no_go": count_no_go,
        "survival_family": {
            "coordinate": (
                "a dimensionless formal subdivision coordinate; this packet "
                "does not identify it with logarithmic physical refinement or "
                "an operational repair clock"
            ),
            "family": "u_m(1 tick) = (1 - (P/24) 2^-m)^(2^m), presence reading",
            "construction_status": (
                "constructed from the imported P/24 density and uniform-splitting "
                "premise, not emitted by the collar records"
            ),
            "records": family,
            "exponential_family_semigroup_defect": {
                "probe": "s = 3/8, t = 5/8 ticks",
                "value": _interval_json(exponential_semigroup_defect),
                "contains_zero": True,
                "width_bound": EXACT_DEFECT_WIDTH_BOUND,
            },
            "thin_representative_exactness": thin_exactness,
        },
        "generator": {
            "status": "imported_conditional_density",
            "statement": (
                "conditional on the declared reserve-trace branch input, "
                "-u_full'(0) = P/24 in an otherwise unbound dimensionless coordinate"
            ),
            "full_collar_derivative": _interval_json(eps),
            "finite_one_step_presence_value": _interval_json(presence_one_step_full),
            "emitted_from_finite_collar_records": False,
            "operational_clock_bound": False,
            "distinction": (
                "the constructed family has sub-step quotient P/24 by definition; "
                "its finite one-step presence value 1 - P/24 is different from "
                "that conditional infinitesimal density"
            ),
        },
        "half_collar_identity": {
            "status": "conditional_on_weighted_coarea_premise",
            "statement": (
                "the counting involution gives equally many forward and backward "
                "slots; theta=P/48 follows only if a source receipt also makes "
                "their reserve weights orientation-invariant"
            ),
            "half_over_full_slot_ratio": "1/2",
            "weighted_half_identity_derived_from_slot_counts": False,
            "orientation_weight_countermodel": weight_countermodel,
            "theta": _interval_json(theta),
            "n_s": _interval_json(n_s),
            "kappa_rep": _interval_json(kappa_rep),
            "kappa_rep_issue_body_reconciliation": {
                "issue_body_value": "2.627023712627471",
                "computed_at": (
                    "P_C = 1.630968209403959, the distinct canonical "
                    "comparison coordinate outside this issue's selected "
                    "P_fwd enclosure"
                ),
                "reproduced_at_P_C_within": "1e-13",
                "inside_certified_interval": False,
                "note": (
                    "the conditional interval uses the imported P enclosure; "
                    "the issue-body decimal tracks the distinct P_C comparison "
                    "coordinate retained by canonical particle surfaces"
                ),
            },
            "finite_one_step_presence_value_half": _interval_json(presence_one_step_half),
        },
        "clock_binding": {
            "status": "open",
            "clock_normalization_source": {
                "invariant": (
                    "m_rep = 24 is a conditional product-adjoint scheduler count"
                ),
                "theorem": (
                    "the issue #343 representation-to-spectrum theorem does not "
                    "identify its scheduler count with the collar refinement "
                    "coordinate or an operational clock"
                ),
                "artifacts": [
                    "code/particles/hierarchy/certificates/R_m_rep_24_certificate.json",
                    "paper/screen_microphysics_and_observer_synchronization.tex",
                ],
            },
            "clock_unit": None,
            "declared_command_line_step_time": None,
            "scheduler_count_bound_to_collar_events": False,
            "dimensionless_coordinate_bound_to_log_refinement": False,
            "operational_clock_derived": False,
            "physical_frequency_derived": False,
            "binding": (
                "no clock binding is proved here; a noncommutative transition "
                "algebra, source state, event-to-refinement map, and operational "
                "calibration remain separate requirements"
            ),
        },
        "defect_bound_summary": {
            "scope": "arithmetic of the constructed conditional survival family",
            "semigroup": (
                "algebraic on the dyadic grid of the constructed family "
                "(interval width gate)"
            ),
            "derivative": (
                "the constructed sub-step quotient equals the imported P/24 "
                "density at every depth"
            ),
            "orientation_balance": (
                "slot-count defect 0 exactly; the displayed survival defect "
                "assumes uniform splitting and does not prove weighted "
                "orientation balance"
            ),
            "refinement": (
                "e^(-P/24) - u_m(1) strictly positive and bounded by "
                "(P/24)_hi^2 2^-m, shrinking with depth"
            ),
            "fail_closed": True,
        },
        "diagnostics": {
            "finite_transition_exponent": {
                "formula": "-log u / log b at u = 1 - P/48, log b = 1 tick",
                "value": _interval_json(finite_transition_exponent),
                "diagnostic_only": True,
                "disjoint_from_theta": True,
                "note": (
                    "substituting the finite one-step survival into the exponent "
                    "formula does not produce the infinitesimal theorem value"
                ),
            },
            "e_branch": {
                "kappa_rep": "e",
                "theta": _interval_json(e_branch_theta),
                "n_s": _interval_json(e_branch_n_s),
                "diagnostic_only": True,
                "note": "a diagnostic alternative, not the selected edge-center coordinate",
            },
        },
        "branch_boundaries": {
            "wrong_orientation_factor_1": {
                "theta": _interval_json(theta_wrong_orientation),
                "status": "not_ruled_out_by_slot_counts_alone",
                "reason": (
                    "equal cardinalities do not fix weighted reserve shares; "
                    "the weighted coarea premise, if independently supplied, "
                    "selects the factor one half"
                ),
                "disjoint_from_selected_theta": True,
            },
            "e_branch_as_primary": {
                "status": "rejected_selection",
                "reason": "the e branch is a diagnostic; selecting it as primary fails closed",
            },
            "poisson_value_at_finite_regulator": {
                "status": "rejected_finite_claim",
                "reason": (
                    "every finite refinement sits strictly below e^(-P/24); the "
                    "Poisson value is the depth limit, not a finite record"
                ),
            },
        },
        "controls": controls,
        "sky_comparison": {
            "performed": False,
            "policy": "this certificate never compares to any sky value",
        },
        "freeze_status": "not_frozen_here",
        "freeze_policy": "freezing the source artifact before sky comparison is an owner action",
        "acceptance_mapping": {
            "clock_normalization_source_identified": {
                "discharged_here": False,
                "detail": (
                    "m_rep=24 is a conditional scheduler count; no theorem here "
                    "binds it to the collar refinement coordinate or an "
                    "operational clock"
                ),
            },
            "full_collar_derivative_and_half_collar_identity": {
                "discharged_here": False,
                "partial_arithmetic_check": True,
                "detail": (
                    "interval propagation and constructed-family defects are "
                    "checked, but P/24 is an imported declared branch input and "
                    "the weighted half identity is an additional premise"
                ),
            },
            "source_graph_ancestry": {
                "discharged_here": False,
                "local_target_injection_rejected": True,
                "detail": (
                    "the local module rejects target-shaped inputs, but it does "
                    "not construct one complete source DAG; P and P/24 enter "
                    "through conditional upstream artifacts and declared branch "
                    "premises"
                ),
            },
            "simulator_emits_p_over_48_only_from_receipt": {
                "discharged_here": False,
                "detail": (
                    "owned by the #579 instantiation; oph_radial_lift_330 currently "
                    "takes theta as a caller input and can consume this receipt"
                ),
            },
            "paper_and_simulator_distinguish_infinitesimal_from_finite": {
                "discharged_here": False,
                "detail": (
                    "the paper fragment states the distinction; the "
                    "simulator-side statement belongs to #579/#580"
                ),
            },
            "all_issue_acceptance_criteria_satisfied": {
                "discharged_here": False,
                "detail": "the source, weighted-half, clock, and simulator gates remain open",
            },
        },
        "claim_boundary": {
            "proves": (
                "exact integer collar counts and rigorous interval arithmetic "
                "for the consequences of imported P/24 and weighted-half "
                "premises; an exhaustive independence no-go for every selector "
                "whose complete input is the integer collar-count tower"
            ),
            "does_not_prove": [
                "that finite collar counts emit the reserve density P/24",
                "that equal orientation counts force equal reserve weights",
                "an operational or physical repair clock",
                "that the simulator consumes theta only from this artifact",
                "a physical or frozen prediction of n_s",
            ],
        },
        "negative_closure_status": {
            "status": "proved_no_go_at_count_only_interface",
            "counts_do_not_select_reserve_density": True,
            "counts_do_not_select_weighted_half": True,
            "counts_do_not_select_generator": True,
            "counts_do_not_select_clock_scale": True,
            "scope": count_no_go["exhaustive_scope"],
            "required_extra_inputs": [
                "source reserve-weight receipt",
                "weighted orientation-reversal coarea identity",
                "source generator receipt",
                "event-to-log-thickness or operational-clock map",
            ],
        },
    }
    payload["artifact_sha256"] = _payload_sha256(payload)
    return payload


def _run_controls(
    weight_countermodel: dict[str, Any],
    count_no_go: dict[str, Any],
) -> dict[str, Any]:
    """Execute the required rejection controls; fail closed if any passes."""
    controls: dict[str, Any] = {}

    try:
        reject_forbidden_inputs({"step_time_seconds": "5.39e-44"})
    except CertificateError as error:
        controls["injected_declared_step_time"] = {"rejected": True, "reason": str(error)}
    else:
        raise CertificateError("control failure: declared step time was accepted")

    try:
        reject_forbidden_inputs({"primordial_amplitude_hint": 0.9649})
    except CertificateError as error:
        controls["injected_measured_tilt_target"] = {"rejected": True, "reason": str(error)}
    else:
        raise CertificateError("control failure: measured-tilt target was accepted")

    try:
        build_check_primary_branch("repair_clock_e")
    except CertificateError as error:
        controls["e_branch_selected_as_primary"] = {"rejected": True, "reason": str(error)}
    else:
        raise CertificateError("control failure: e branch was accepted as primary")

    if not (
        weight_countermodel[
            "equal_weight_forward_reserve_share_at_every_depth"
        ]
        == "1/2"
        and weight_countermodel[
            "biased_positive_weight_forward_reserve_share_at_every_depth"
        ]
        == "2/3"
        and weight_countermodel["all_weights_strictly_positive"] is True
    ):
        raise CertificateError("orientation-weight countermodel failed")
    controls["slot_count_only_weighted_half_claim"] = {
        "rejected": True,
        "reason": (
            "equal forward/backward slot counts admit both one-half and "
            "two-thirds forward reserve shares under strictly positive weights"
        ),
        "countermodel": weight_countermodel,
    }
    for name, action, marker in (
        (
            "count_only_unique_reserve_share_claim",
            assert_counts_select_unique_reserve_share,
            "COUNT_WEIGHT_INDEPENDENCE",
        ),
        (
            "count_only_unique_generator_claim",
            assert_counts_select_unique_generator,
            "COUNT_GENERATOR_INDEPENDENCE",
        ),
        (
            "count_only_unique_clock_scale_claim",
            assert_counts_select_unique_clock_scale,
            "COUNT_CLOCK_INDEPENDENCE",
        ),
    ):
        try:
            action(count_no_go)
        except CertificateError as error:
            if marker not in str(error):
                raise
            controls[name] = {"rejected": True, "reason": str(error)}
        else:
            raise CertificateError(f"control failure: {name} was accepted")
    return controls


def build_check_primary_branch(primary_branch: str) -> None:
    """Reject any primary-branch selection other than the theorem branch."""
    if primary_branch != SELECTED_PRIMARY_BRANCH:
        raise CertificateError(
            f"primary branch must be {SELECTED_PRIMARY_BRANCH!r}; "
            f"{primary_branch!r} branches are diagnostics only"
        )


def verify_stored_artifact(path: Path, expected: dict[str, Any] | None = None) -> None:
    """Fail closed unless a stored artifact matches a fresh build and self-hash."""

    if not path.exists():
        raise CertificateError(f"missing stored artifact: {path}")
    stored = json.loads(path.read_text(encoding="utf-8"))
    verify_artifact_self_hash(stored)
    fresh = expected if expected is not None else build()
    verify_artifact_self_hash(fresh)
    if stored != fresh:
        raise CertificateError("stored artifact is stale, malformed, or tampered")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Emit the conditional edge-center arithmetic and count-only "
            "independence packet (#522)."
        )
    )
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    parser.add_argument("--print-json", action="store_true")
    parser.add_argument(
        "--check",
        action="store_true",
        help="verify the stored output against a fresh build without rewriting it",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    payload = build()
    out_path = Path(args.output)
    if args.check:
        verify_stored_artifact(out_path, payload)
        print(f"verified: {out_path}")
        return 0
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text, encoding="utf-8")
    if args.print_json:
        print(text, end="")
    else:
        print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
