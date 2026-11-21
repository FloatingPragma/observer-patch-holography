#!/usr/bin/env python3
"""Build the cosmology postdiction ledger.

Every row compares a corpus-side value, read from a typed corpus input or from
a pinned parent receipt, with a public measurement typed once in
``public_inputs.json``.  All comparisons are on seen data.  The ledger
promotes nothing: no row is a frozen prediction, a score, or evidence for or
against OPH.

Outputs, written deterministically (sorted keys, two-space indent, ASCII, LF,
no timestamps, no absolute paths):

* ``runtime/cosmology_postdiction_ledger.json``
* ``COSMOLOGY_POSTDICTION_LEDGER.md``, rendered from the JSON.

Exact rational arithmetic on the typed decimal digits is used wherever the
theory value is rational; float arithmetic in a fixed operation order is used
for the irrational constants (golden ratio, e, pi, square roots).  The
independent verifier replays both without importing this module.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
from decimal import Decimal
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import strict_json  # noqa: E402

SCHEMA = "oph.cosmology.postdiction_ledger.v1"
PUBLIC_SCHEMA = "oph.cosmology.postdiction_ledger.public_inputs.v1"
CORPUS_SCHEMA = "oph.cosmology.postdiction_ledger.corpus_inputs.v1"
CLASSIFICATION = (
    "retrospective postdiction ledger; every comparison is on seen data; "
    "no row is a frozen prediction or a score"
)
PACKAGE = "code/cosmology/postdiction_ledger"
PUBLIC_INPUTS_PATH = HERE / "public_inputs.json"
CORPUS_INPUTS_PATH = HERE / "corpus_inputs.json"
RECEIPT_NAME = "runtime/cosmology_postdiction_ledger.json"
MARKDOWN_NAME = "COSMOLOGY_POSTDICTION_LEDGER.md"
PINNED_CODE = (
    f"{PACKAGE}/build_cosmology_postdiction_ledger.py",
    f"{PACKAGE}/verify_cosmology_postdiction_ledger_independent.py",
    f"{PACKAGE}/test_cosmology_postdiction_ledger.py",
)
ALLOWED_CLASSES = (
    "conditional_theorem_postdiction",
    "shared_baseline",
    "closure_candidate_display",
    "fitted_comparison_value",
    "not_evaluable",
)
ALLOWED_VERDICTS = (
    "consistent",
    "tension",
    "exceeds_three_sigma_diagnostic",
    "exceeds_stated_confidence_bound",
    "not_evaluable",
)
VERDICT_RULE = {
    "sigma_rows": (
        "|sigma| < 2 consistent; 2 <= |sigma| < 3 tension; "
        "|sigma| >= 3 exceeds_three_sigma_diagnostic"
    ),
    "bound_rows": (
        "consistent exactly when the theory value lies strictly below the "
        "quoted upper bound; otherwise exceeds_stated_confidence_bound"
    ),
    "rows_without_comparison_contract": "not_evaluable",
}
# The frozen FZ-01 capacity band (claims/frozen_prediction_register.json, FZ-01
# decision policy); a registered target, rendered beside the display candidates.
FZ01_CAPACITY_BAND = (3.30836e122, 3.31137e122)
SIGMA_SIGN_CONVENTION = (
    "theory minus measurement, divided by the quoted one-sigma uncertainty"
)
CLASS_DESCRIPTIONS = {
    "conditional_theorem_postdiction": (
        "value follows from a corpus theorem under named premises and was "
        "exposed after the data"
    ),
    "shared_baseline": (
        "identical to the standard LambdaCDM or single-field expectation, so "
        "agreement supports nothing"
    ),
    "closure_candidate_display": "target-informed declared map",
    "fitted_comparison_value": "fitted comparison value with no source value",
    "not_evaluable": "no comparison contract",
}
ROW_FIELDS = (
    "row_id",
    "section",
    "quantity",
    "class",
    "oph_value",
    "oph_value_display",
    "oph_value_exact_rational",
    "oph_premises",
    "measurement",
    "measurement_display",
    "dataset",
    "citation",
    "public_input_ids",
    "corpus_source_ids",
    "sigma_distance",
    "sigma_distance_kind",
    "verdict",
    "discriminates_from_baseline",
    "prospective_data",
    "epistemic_note",
)
NS_DATASETS = (
    "planck2018_vi_ns",
    "spt3g_d1_planck_ns",
    "cmb_spa_ns",
    "spt3g_d1_alone_ns",
    "act_dr6_p_act_ns",
    "act_dr6_p_act_lb_ns",
)
DESI_COMBINATIONS = (
    ("DESI_DR2_BAO+CMB", "desi_dr2_bao_cmb"),
    ("DESI_DR2_BAO+CMB+PantheonPlus", "desi_dr2_bao_cmb_pantheonplus"),
    ("DESI_DR2_BAO+CMB+Union3", "desi_dr2_bao_cmb_union3"),
    ("DESI_DR2_BAO+CMB+DESY5", "desi_dr2_bao_cmb_desy5"),
)
DESI_CITATION = (
    "DESI Collaboration (2025), DESI DR2 Results II, arXiv:2503.14738; "
    "official Cobaya chains, "
    "https://data.desi.lbl.gov/public/papers/y3/bao-cosmo-params/"
)


class LedgerError(RuntimeError):
    """Raised when an input, pin, anchor, or cross-check fails."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise LedgerError(message)


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def canonical_bytes(payload: Any) -> bytes:
    return json.dumps(
        payload, ensure_ascii=True, separators=(",", ":"), sort_keys=True
    ).encode("utf-8")


def render_bytes(payload: Any) -> bytes:
    return (
        json.dumps(payload, ensure_ascii=True, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")


def load_exact_bytes(payload: bytes) -> Any:
    """Strictly decode JSON, keeping every decimal literal exact."""

    text = payload.decode("utf-8")
    strict_json.loads(text)
    return json.loads(text, parse_float=Decimal)


def frac(value: Any) -> Fraction:
    require(
        isinstance(value, (int, Decimal)) and not isinstance(value, bool),
        f"exact number required, got {type(value).__name__}",
    )
    return Fraction(value)


def normalize_whitespace(text: str) -> str:
    return re.sub(r"\s+", " ", text)


def rational_string(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def verdict_from_sigma(sigma: Fraction | float) -> str:
    magnitude = abs(sigma)
    if magnitude < 2:
        return "consistent"
    if magnitude < 3:
        return "tension"
    return "exceeds_three_sigma_diagnostic"


def verdict_from_bound(theory: Fraction | float, bound: Fraction | float) -> str:
    return "consistent" if theory < bound else "exceeds_stated_confidence_bound"


def signed_sigma(theory: Fraction, value: Fraction, uncertainty: Fraction) -> Fraction:
    require(uncertainty > 0, "positive uncertainty required")
    return (theory - value) / uncertainty


def fmt_num(value: Any) -> str:
    return f"{float(value):g}"


def fmt_sigma(sigma: float | None) -> str:
    return "n/a" if sigma is None else f"{sigma:+.2f}"


def fmt_a0(value: float) -> str:
    return f"{value / 1e-10:.4f}e-10"


def check_anchors(typed: dict[str, Any], root: Path) -> list[dict[str, Any]]:
    checks = []
    for key in typed:
        entry = typed[key]
        source = entry["source_file"]
        anchor = entry["anchor_text"]
        if source is None:
            require(anchor is None, f"{key}: anchor without source file")
            continue
        require(not Path(source).is_absolute(), f"{key}: absolute source path")
        text = (root / source).read_text(encoding="utf-8")
        require(
            normalize_whitespace(anchor) in normalize_whitespace(text),
            f"{key}: anchor text absent from {source}",
        )
        checks.append(
            {"corpus_input_id": key, "source_file": source, "anchor_verified": True}
        )
    return checks


def load_parents(spec: dict[str, Any], root: Path) -> dict[str, Any]:
    parents = {}
    for key in spec:
        entry = spec[key]
        path = root / entry["path"]
        payload = path.read_bytes()
        digest = sha256_bytes(payload)
        require(digest == entry["sha256"], f"{key}: parent receipt sha256 mismatch")
        data = load_exact_bytes(payload)
        require(data["schema"] == entry["schema"], f"{key}: parent schema mismatch")
        parents[key] = data
    return parents


def derive_constants(
    items: dict[str, Any], typed: dict[str, Any], parents: dict[str, Any]
) -> dict[str, Any]:
    c = float(items["si_speed_of_light"]["value"])
    mpc = float(items["iau2015_megaparsec"]["value"])
    lp = float(items["codata2022_planck_length"]["value"])
    h0 = float(items["planck2018_vi_h0"]["value"])
    omega_lambda = float(items["planck2018_vi_omega_lambda"]["value"])
    desi_constants = parents["official_desi_dr2_fz13_retrospective"][
        "base_lcdm_capacity_display"
    ]["constants"]
    require(
        float(desi_constants["Planck_length_m_CODATA_2022_central"]) == lp
        and float(desi_constants["speed_of_light_m_s_exact"]) == c
        and float(desi_constants["Mpc_in_m"]) == mpc,
        "constants differ between public inputs and the DESI parent receipt",
    )

    h0_si = h0 * 1000.0 / mpc
    lambda_planck = 3.0 * omega_lambda * (h0_si / c) ** 2
    lambda_lp2_planck = lambda_planck * lp**2
    n_planck = 3.0 * math.pi / lambda_lp2_planck
    a_ds = c**2 * math.sqrt(lambda_planck / 3.0)
    c_h0 = c * h0_si
    flagship_lambda = typed["flagship_lambda_planck_central_display"]
    digits = int(flagship_lambda["significant_digits"])
    require(
        float(f"{lambda_planck:.{digits - 1}e}") == float(flagship_lambda["value"]),
        "Planck-central Lambda disagrees with the flagship display",
    )

    p_c = frac(typed["P_C"]["value"])
    p_fwd = frac(typed["P_fwd"]["value"])
    ns_p_c = 1 - p_c / 48
    ns_p_fwd = 1 - p_fwd / 48
    phi = (1.0 + math.sqrt(5.0)) / 2.0
    ns_clock = 1.0 - math.e * (float(p_c) - phi)
    edge_display = typed["edge_center_paper_display"]
    require(
        round(float(ns_p_c), int(edge_display["decimals"]))
        == float(edge_display["value"]),
        "edge-center n_s disagrees with the inflation paper display",
    )
    clock_display = typed["clock_branch_paper_display"]
    require(
        round(ns_clock, int(clock_display["decimals"])) == float(clock_display["value"]),
        "clock-branch n_s disagrees with the inflation paper display",
    )

    desi_display = parents["official_desi_dr2_fz13_retrospective"][
        "base_lcdm_capacity_display"
    ]["combined"]
    desi_lambda_lp2_mean = float(desi_display["Lambda_lP2"]["weighted_mean"])
    desi_lambda_lp2_std = float(desi_display["Lambda_lP2"]["weighted_std"])
    n_desi = 3.0 * math.pi / desi_lambda_lp2_mean

    return {
        "speed_of_light_m_s": c,
        "megaparsec_m": mpc,
        "planck_length_m": lp,
        "planck2018_H0_km_s_Mpc": h0,
        "planck2018_Omega_Lambda": omega_lambda,
        "planck2018_H0_s_minus_1": h0_si,
        "Lambda_planck_central_m_minus_2": lambda_planck,
        "Lambda_lP2_planck_central": lambda_lp2_planck,
        "capacity_coordinate_planck_central": n_planck,
        "a_dS_m_s2": a_ds,
        "c_H0_m_s2": c_h0,
        "P_C": float(p_c),
        "P_fwd": float(p_fwd),
        "golden_ratio": phi,
        "euler_number": math.e,
        "ns_edge_center_P_C": float(ns_p_c),
        "ns_edge_center_P_fwd": float(ns_p_fwd),
        "ns_edge_center_branch_spread": float(ns_p_fwd - ns_p_c),
        "ns_clock_branch_P_C": ns_clock,
        "Lambda_lP2_desi_dr2_display_mean": desi_lambda_lp2_mean,
        "Lambda_lP2_desi_dr2_display_std": desi_lambda_lp2_std,
        "capacity_coordinate_desi_dr2_display": n_desi,
        "_exact": {"P_C": p_c, "P_fwd": p_fwd, "ns_P_C": ns_p_c, "ns_P_fwd": ns_p_fwd},
    }


def make_row(**fields: Any) -> dict[str, Any]:
    extras = fields.pop("extras", None)
    require(set(fields) == set(ROW_FIELDS), f"row fields mismatch: {sorted(set(fields) ^ set(ROW_FIELDS))}")
    require(fields["class"] in ALLOWED_CLASSES, "unknown row class")
    require(fields["verdict"] in ALLOWED_VERDICTS, "unknown verdict")
    require(isinstance(fields["discriminates_from_baseline"], bool), "discriminates flag")
    row = dict(fields)
    if extras:
        row["extras"] = extras
    return row


def gaussian_row(
    *,
    row_id: str,
    section: str,
    quantity: str,
    row_class: str,
    theory: Fraction | float,
    theory_display: str,
    premises: list[str],
    item_id: str,
    item: dict[str, Any],
    corpus_ids: list[str],
    discriminates: bool,
    prospective: str,
    note: str,
    extras: dict[str, Any] | None = None,
) -> dict[str, Any]:
    value = frac(item["value"])
    uncertainty = frac(item["uncertainty"])
    if isinstance(theory, Fraction):
        sigma_exact = signed_sigma(theory, value, uncertainty)
        sigma = float(sigma_exact)
        verdict = verdict_from_sigma(sigma_exact)
        exact = rational_string(theory)
        oph_value: float = float(theory)
    else:
        sigma = (theory - float(value)) / float(uncertainty)
        verdict = verdict_from_sigma(sigma)
        exact = None
        oph_value = theory
    unit = item.get("unit")
    return make_row(
        row_id=row_id,
        section=section,
        quantity=quantity,
        **{"class": row_class},
        oph_value=oph_value,
        oph_value_display=theory_display,
        oph_value_exact_rational=exact,
        oph_premises=premises,
        measurement={
            "kind": "gaussian",
            "value": float(value),
            "sigma": float(uncertainty),
            "unit": unit,
        },
        measurement_display=f"{fmt_num(value)} +/- {fmt_num(uncertainty)}"
        + (f" {unit}" if unit else ""),
        dataset=item["dataset"],
        citation=item["citation"],
        public_input_ids=[item_id],
        corpus_source_ids=corpus_ids,
        sigma_distance=sigma,
        sigma_distance_kind="signed_one_dimensional",
        verdict=verdict,
        discriminates_from_baseline=discriminates,
        prospective_data=prospective,
        epistemic_note=note,
        extras=extras,
    )


def bound_row(
    *,
    row_id: str,
    section: str,
    quantity: str,
    row_class: str,
    theory: Fraction,
    theory_display: str,
    premises: list[str],
    item_id: str,
    item: dict[str, Any],
    corpus_ids: list[str],
    discriminates: bool,
    prospective: str,
    note: str,
    extras: dict[str, Any] | None = None,
) -> dict[str, Any]:
    bound = frac(item["upper_bound"])
    confidence = frac(item["confidence_level"])
    measurement: dict[str, Any] = {
        "kind": "upper_bound",
        "upper_bound": float(bound),
        "confidence_level": float(confidence),
        "unit": item.get("unit"),
    }
    if "k_pivot_mpc_inverse" in item:
        measurement["k_pivot_mpc_inverse"] = float(item["k_pivot_mpc_inverse"])
    if "sigma_r" in item:
        measurement["sigma_r"] = float(item["sigma_r"])
    display = f"< {fmt_num(bound)} ({fmt_num(confidence * 100)} percent CL)"
    if "k_pivot_mpc_inverse" in item:
        display += f" at k = {fmt_num(item['k_pivot_mpc_inverse'])} Mpc^-1"
    return make_row(
        row_id=row_id,
        section=section,
        quantity=quantity,
        **{"class": row_class},
        oph_value=float(theory),
        oph_value_display=theory_display,
        oph_value_exact_rational=rational_string(theory),
        oph_premises=premises,
        measurement=measurement,
        measurement_display=display,
        dataset=item["dataset"],
        citation=item["citation"],
        public_input_ids=[item_id],
        corpus_source_ids=corpus_ids,
        sigma_distance=None,
        sigma_distance_kind="none_upper_bound",
        verdict=verdict_from_bound(theory, bound),
        discriminates_from_baseline=discriminates,
        prospective_data=prospective,
        epistemic_note=note,
        extras=extras,
    )


def primordial_rows(items: dict[str, Any], consts: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    exact = consts["_exact"]
    spread = consts["ns_edge_center_branch_spread"]
    edge_note = (
        "Conditional theorem: n_s = 1 - P/48 follows from the declared "
        "full-collar generator density P/24 and the orientation-half identity; "
        "the source generator receipt is not supplied. Both P values were "
        f"exposed after the data. The P_fwd minus P_C branch spread in n_s is {spread:.1e}, "
        "far below every quoted uncertainty. Seen data; the row is a "
        "postdiction, not a frozen prediction."
    )
    for branch, p_value, ns_value, p_label in (
        ("P_C", exact["P_C"], exact["ns_P_C"], "the CODATA-located comparison pixel"),
        ("P_fwd", exact["P_fwd"], exact["ns_P_fwd"], "the source closure-map root"),
    ):
        premises = [
            "full-collar reserve-generator density P/24 (declared branch input)",
            "orientation-half identity giving the half-collar density P/48",
            "continuous dilation cocycle on the source-facing half collar",
            f"P = {branch} = {float(p_value)!r} ({p_label})",
        ]
        for item_id in NS_DATASETS:
            rows.append(
                gaussian_row(
                    row_id=f"ns_edge_center_{branch}_{item_id}",
                    section="primordial",
                    quantity=f"n_s = 1 - {branch}/48 (edge-center branch)",
                    row_class="conditional_theorem_postdiction",
                    theory=ns_value,
                    theory_display=f"{float(ns_value):.6f}",
                    premises=premises,
                    item_id=item_id,
                    item=items[item_id],
                    corpus_ids=[branch, "edge_center_tilt_theorem", "edge_center_paper_display"],
                    discriminates=True,
                    prospective="Simons Observatory",
                    note=edge_note,
                )
            )
    clock_premises = [
        "clock-branch coordinate theta = e (P - varphi)",
        f"P = P_C = {consts['P_C']!r}",
        "e taken as a separate diagnostic hypothesis, not a source-derived coordinate",
    ]
    clock_note = (
        "diagnostic alternative branch; no selection between branches is made "
        "by data. The inflation paper displays this branch value as 0.96484. "
        "Float arithmetic: 1 - e (P_C - varphi) with varphi = (1 + sqrt 5)/2."
    )
    for item_id in NS_DATASETS:
        rows.append(
            gaussian_row(
                row_id=f"ns_clock_branch_P_C_{item_id}",
                section="primordial",
                quantity="n_s = 1 - e (P_C - varphi) (clock branch, diagnostic)",
                row_class="conditional_theorem_postdiction",
                theory=consts["ns_clock_branch_P_C"],
                theory_display=f"{consts['ns_clock_branch_P_C']:.6f}",
                premises=clock_premises,
                item_id=item_id,
                item=items[item_id],
                corpus_ids=["P_C", "golden_ratio", "euler_number", "clock_branch_diagnostic_alternative", "clock_branch_paper_display"],
                discriminates=True,
                prospective="Simons Observatory",
                note=clock_note,
            )
        )

    running_note = (
        "Zero running is exact for the fixed dilation exponent (pure power "
        "law), conditional on that exponent. Generic slow-roll "
        "running is second order in the slow-roll parameters, of order 1e-3, "
        "so a running measurement at that level separates the branches. Seen data."
    )
    for item_id in ("planck2018_running", "act_dr6_p_act_lb_running"):
        rows.append(
            gaussian_row(
                row_id=f"running_zero_{item_id}",
                section="primordial",
                quantity="dn_s/dlnk",
                row_class="conditional_theorem_postdiction",
                theory=Fraction(0),
                theory_display="0 (exact)",
                premises=[
                    "covariance-survival cocycle with a fixed dilation exponent",
                    "physical dilation intertwiner transferring the cocycle to the curvature covariance",
                    "homogeneous radial source family: exact power law",
                ],
                item_id=item_id,
                item=items[item_id],
                corpus_ids=["zero_running_dilation_exponent", "homogeneous_radial_source_family_theorem"],
                discriminates=True,
                prospective="Simons Observatory",
                note=running_note,
            )
        )

    tensor_note = (
        "Zero linear primordial tensor amplitude requires zero initial tensor "
        "amplitude and velocity, no tensor forcing, and positive scalar power. "
        "A scalar source alone leaves homogeneous tensor modes unconstrained. "
        "These extra initial conditions are branch premises, not source-derived "
        "results. Slow roll has no universal positive lower bound on r; an upper "
        "limit tests only specified alternatives above it. Observable comparisons "
        "require separate treatment of scalar-induced tensors, lensing and "
        "foregrounds. Bound rows carry no sigma distance. Seen data."
    )
    for item_id in ("bk18_r_0p05", "tristram2022_r", "planck2018_x_r_0p002_bk15"):
        rows.append(
            bound_row(
                row_id=f"tensor_zero_{item_id}",
                section="primordial",
                quantity=items[item_id]["quantity"] + " (tensor-to-scalar ratio)",
                row_class="conditional_theorem_postdiction",
                theory=Fraction(0),
                theory_display="0 (exact)",
                premises=[
                    "linear tensor evolution on the declared background",
                    "zero initial tensor amplitude and velocity",
                    "no tensor forcing",
                    "positive scalar power at the comparison pivot",
                ],
                item_id=item_id,
                item=items[item_id],
                corpus_ids=["zero_tensor_scalar_screen_field"],
                discriminates=False,
                prospective="LiteBIRD (tests specified tensor alternatives above the achieved bound)",
                note=tensor_note,
            )
        )

    iso_note = (
        "Zero source-side isocurvature follows for a Gaussian single-clock "
        "source. Adiabatic initial conditions are the single-field "
        "baseline as well, so a bound satisfied by zero supports nothing. "
        "Seen data."
    )
    for item_id in (
        "planck2018_x_cdi_axion_i_100beta_klow",
        "planck2018_x_cdi_axion_i_100beta_kmid",
        "planck2018_x_cdi_axion_i_100beta_khigh",
        "planck2018_x_cdi_general_100beta_klow",
    ):
        rows.append(
            bound_row(
                row_id=f"isocurvature_zero_{item_id}",
                section="primordial",
                quantity=items[item_id]["quantity"],
                row_class="shared_baseline",
                theory=Fraction(0),
                theory_display="0 (exact)",
                premises=["Gaussian single-clock source", "one common release clock"],
                item_id=item_id,
                item=items[item_id],
                corpus_ids=["zero_isocurvature_single_clock"],
                discriminates=False,
                prospective="LiteBIRD (bound tightening only; zero is the shared baseline)",
                note=iso_note,
            )
        )

    fnl_note = (
        "The Gaussian screen covariance gives zero primordial non-Gaussianity at "
        "the source; the standard nonlinear transfer of order unity is not "
        "computed here. Gaussian initial conditions are the single-field "
        "baseline as well, so agreement supports nothing. Seen data."
    )
    for item_id in ("planck2018_ix_fnl_local", "planck2018_ix_fnl_equil", "planck2018_ix_fnl_ortho"):
        rows.append(
            gaussian_row(
                row_id=f"fnl_zero_{item_id}",
                section="primordial",
                quantity=items[item_id]["quantity"],
                row_class="shared_baseline",
                theory=Fraction(0),
                theory_display="0 at the source (transfer not computed)",
                premises=["screen covariance theorem: the entropy-maximizing law at fixed energy is Gaussian", "nonlinear transfer of order unity not computed"],
                item_id=item_id,
                item=items[item_id],
                corpus_ids=["gaussian_maxent_release"],
                discriminates=False,
                prospective="none scoring; large-scale-structure f_NL programmes tighten a shared-baseline bound",
                note=fnl_note,
            )
        )

    curvature_note = (
        "Zero curvature follows from the small-loop holonomy identification of "
        "a flat branch, conditional on that identification. Flatness is the "
        "LambdaCDM baseline as well, so agreement supports nothing. The DESI "
        "DR2+CMB row sits above two sigma under the fixed rule and is reported "
        "as it stands. Seen data."
    )
    for item_id in ("planck2018_vi_omegak", "planck2018_vi_omegak_bao", "desi_dr2_cmb_omegak"):
        rows.append(
            gaussian_row(
                row_id=f"curvature_zero_{item_id}",
                section="primordial",
                quantity="Omega_K",
                row_class="shared_baseline",
                theory=Fraction(0),
                theory_display="0 (exact)",
                premises=[
                    "vanishing area-normalized small-loop holonomy",
                    "clock-slice flat FLRW branch identification",
                ],
                item_id=item_id,
                item=items[item_id],
                corpus_ids=["flat_branch_small_loop_holonomy"],
                discriminates=False,
                prospective="DESI DR3; Euclid",
                note=curvature_note,
            )
        )
    return rows


def dark_energy_rows(parents: dict[str, Any]) -> list[dict[str, Any]]:
    desi = parents["official_desi_dr2_fz13_retrospective"]
    rows = []
    note = (
        "Gaussian moment diagnostic on the official chains: the Mahalanobis "
        "distance of the point (-1, 0) from the posterior mean under the "
        "sampled w0-wa covariance, converted to a two-sided normal "
        "equivalent; it is not the collaboration's delta chi-squared or "
        "evidence. Seen data; direction-neutral. The point (-1, 0) is also the "
        "LambdaCDM null, so consistency would support nothing; tension rows are "
        "reported as they stand. The posterior fractions are empirical masses "
        "under the DESI model, priors, and likelihoods."
    )
    premises = [
        "density, continuity, and closed-sector premises of the capacity law",
        "fixed record capacity: d ln N / d ln a = 0",
        "exposure map w(a) = -1 + (1/3) d ln N / d ln a (Lean FixedCapacityWLaw)",
    ]
    for dataset_key, slug in DESI_COMBINATIONS:
        block = desi["datasets"][dataset_key]
        combined = block["combined"]
        diagnostic = block["fixed_capacity_point_gaussian_diagnostic"]
        w0_mean = frac(combined["w0_mean"])
        w0_std = frac(combined["w0_std"])
        wa_mean = frac(combined["wa_mean"])
        wa_std = frac(combined["wa_std"])
        sigma_equivalent = float(diagnostic["two_sided_normal_sigma_equivalent"])
        rows.append(
            make_row(
                row_id=f"w0wa_fixed_capacity_{slug}",
                section="background_dark_energy",
                quantity="(w0, wa) fixed-capacity point",
                **{"class": "conditional_theorem_postdiction"},
                oph_value={"w0": -1.0, "wa": 0.0},
                oph_value_display="(w0, wa) = (-1, 0) exact",
                oph_value_exact_rational=None,
                oph_premises=premises,
                measurement={
                    "kind": "two_parameter_gaussian_moments",
                    "w0_mean": float(w0_mean),
                    "w0_std": float(w0_std),
                    "wa_mean": float(wa_mean),
                    "wa_std": float(wa_std),
                    "w0_wa_covariance": float(combined["w0_wa_covariance"]),
                    "w0_wa_correlation": float(combined["w0_wa_correlation"]),
                },
                measurement_display=(
                    f"w0 = {float(w0_mean):.3f} +/- {float(w0_std):.3f}, "
                    f"wa = {float(wa_mean):.3f} +/- {float(wa_std):.3f}"
                ),
                dataset=f"{dataset_key} (official DESI DR2 base_w_wa chains)",
                citation=DESI_CITATION,
                public_input_ids=[],
                corpus_source_ids=["fixed_capacity_cpl_point", "official_desi_dr2_fz13_retrospective"],
                sigma_distance=sigma_equivalent,
                sigma_distance_kind="two_dof_gaussian_mahalanobis_two_sided_normal_equivalent",
                verdict=verdict_from_sigma(sigma_equivalent),
                discriminates_from_baseline=True,
                prospective_data="DESI DR3; Euclid",
                epistemic_note=note,
                extras={
                    "mahalanobis_squared": float(diagnostic["mahalanobis_squared"]),
                    "chi2_2dof_survival": float(diagnostic["chi2_2dof_survival"]),
                    "component_sigma_distances": {
                        "w0": float(signed_sigma(Fraction(-1), w0_mean, w0_std)),
                        "wa": float(signed_sigma(Fraction(0), wa_mean, wa_std)),
                    },
                    "posterior_mass_w_ge_minus_one_for_0_le_z_le_2": float(
                        combined["posterior_mass_w_ge_minus_one_for_0_le_z_le_2"]
                    ),
                    "posterior_mass_w0_gt_minus_one": float(combined["posterior_mass_w0_gt_minus_one"]),
                    "posterior_mass_wa_nonnegative": float(combined["posterior_mass_wa_nonnegative"]),
                },
            )
        )
    return rows


def capacity_rows(
    items: dict[str, Any], typed: dict[str, Any], parents: dict[str, Any], consts: dict[str, Any]
) -> list[dict[str, Any]]:
    rows = []
    lp = consts["planck_length_m"]
    lambda_planck = consts["Lambda_planck_central_m_minus_2"]
    n_planck = consts["capacity_coordinate_planck_central"]
    desi_mean = consts["Lambda_lP2_desi_dr2_display_mean"]
    desi_std = consts["Lambda_lP2_desi_dr2_display_std"]
    n_desi = consts["capacity_coordinate_desi_dr2_display"]
    desi_display = parents["official_desi_dr2_fz13_retrospective"]["base_lcdm_capacity_display"]
    quantiles = {
        key: float(value)
        for key, value in desi_display["combined"]["Lambda_lP2"]["weighted_step_cdf_quantiles"].items()
    }
    desi_h0 = float(desi_display["combined"]["H0_km_s_Mpc"]["weighted_mean"])
    desi_omega_lambda = float(desi_display["combined"]["OmegaLambda"]["weighted_mean"])
    percent_display = [float(v) for v in typed["capacity_lambda_percent_display"]["values"]]
    premises = [
        "capacity read as de Sitter horizon entropy in nats: Lambda = 3 pi/(N l_P^2)",
        "candidate N from the closure hypothesis (declared, target-informed map)",
        "calibration import of hbar, G, and c through l_P",
    ]
    note = (
        "Target-informed declared map: the capacity coordinate N = 3 pi/(Lambda "
        "l_P^2) is read from the data and both candidates were exposed after it, "
        "so the residual carries no predictive weight. Large N explains the small "
        "dimensional constant, not the closeness of the residual: a relative "
        "capacity residual delta gives a relative constant residual -delta/(1+delta). "
        "The DESI display is a base-LambdaCDM posterior display whose centrals "
        f"(H0 = {desi_h0:.2f}, Omega_Lambda = {desi_omega_lambda:.4f}) differ from the "
        "Planck 2018 centrals; that shift moves the coordinate by about four "
        "percent, larger than the candidate spread. Seen data."
    )
    for index, candidate_id in enumerate(("capacity_candidate_a", "capacity_candidate_b")):
        n_candidate = float(typed[candidate_id]["value"])
        lambda_candidate = 3.0 * math.pi / (n_candidate * lp**2)
        lambda_lp2_candidate = 3.0 * math.pi / n_candidate
        pct_lambda_planck = 100.0 * (lambda_candidate - lambda_planck) / lambda_planck
        pct_capacity_planck = 100.0 * (n_candidate - n_planck) / n_planck
        require(
            round(pct_lambda_planck, 2) == percent_display[index],
            f"{candidate_id}: percent residual disagrees with the paper display",
        )
        rows.append(
            make_row(
                row_id=f"lambda_{candidate_id}_planck2018_central",
                section="capacity_lambda",
                quantity=f"Lambda from N = {n_candidate:.4e} (candidate {'a' if index == 0 else 'b'})",
                **{"class": "closure_candidate_display"},
                oph_value=lambda_candidate,
                oph_value_display=f"{lambda_candidate:.4e} m^-2",
                oph_value_exact_rational=None,
                oph_premises=premises,
                measurement={
                    "kind": "central_value_display",
                    "value": lambda_planck,
                    "sigma": None,
                    "unit": "m^-2",
                    "inputs": {
                        "H0_km_s_Mpc": consts["planck2018_H0_km_s_Mpc"],
                        "Omega_Lambda": consts["planck2018_Omega_Lambda"],
                    },
                    "capacity_coordinate": n_planck,
                },
                measurement_display=f"{lambda_planck:.4e} m^-2 (central; no uncertainty attached)",
                dataset="Planck 2018 TT,TE,EE+lowE+lensing base-LambdaCDM centrals",
                citation=items["planck2018_vi_h0"]["citation"],
                public_input_ids=["planck2018_vi_h0", "planck2018_vi_omega_lambda", "codata2022_planck_length", "si_speed_of_light", "iau2015_megaparsec"],
                corpus_source_ids=[candidate_id, "capacity_lambda_dictionary", "flagship_lambda_planck_central_display", "capacity_lambda_percent_display"],
                sigma_distance=None,
                sigma_distance_kind="none_central_value_display",
                verdict="not_evaluable",
                discriminates_from_baseline=False,
                prospective_data="DESI DR3; Euclid (background posterior displays only)",
                epistemic_note=note + " Central-value display without an uncertainty; percent residuals only.",
                extras={
                    "percent_residual_in_lambda": pct_lambda_planck,
                    "percent_residual_in_capacity": pct_capacity_planck,
                    "planck_chain_capacity_coordinate_display": float(
                        typed["planck_chain_capacity_coordinate_display"]["value"]
                    ),
                },
            )
        )
        sigma = (lambda_lp2_candidate - desi_mean) / desi_std
        rows.append(
            make_row(
                row_id=f"lambda_lp2_{candidate_id}_desi_dr2_display",
                section="capacity_lambda",
                quantity=f"Lambda l_P^2 = 3 pi/N, N = {n_candidate:.4e} (candidate {'a' if index == 0 else 'b'})",
                **{"class": "closure_candidate_display"},
                oph_value=lambda_lp2_candidate,
                oph_value_display=f"{lambda_lp2_candidate:.5e}",
                oph_value_exact_rational=None,
                oph_premises=premises,
                measurement={
                    "kind": "gaussian",
                    "value": desi_mean,
                    "sigma": desi_std,
                    "unit": None,
                    "weighted_step_cdf_quantiles": quantiles,
                    "capacity_coordinate": n_desi,
                },
                measurement_display=f"{desi_mean:.5e} +/- {desi_std:.5e} (weighted mean and std)",
                dataset="DESI DR2 BAO + CMB flat base-LambdaCDM chains, sample-level Lambda l_P^2 display",
                citation=DESI_CITATION,
                public_input_ids=[],
                corpus_source_ids=[candidate_id, "capacity_lambda_dictionary", "official_desi_dr2_fz13_retrospective"],
                sigma_distance=None,
                sigma_distance_kind="none_not_evaluable",
                verdict="not_evaluable",
                discriminates_from_baseline=False,
                prospective_data="DESI DR3; Euclid (background posterior displays only)",
                epistemic_note=note
                + " A target-informed display carries no verdict; the standardized"
                " distance is display information only.",
                extras={
                    "percent_residual_in_lambda": 100.0 * (lambda_lp2_candidate - desi_mean) / desi_mean,
                    "percent_residual_in_capacity": 100.0 * (n_candidate - n_desi) / n_desi,
                    "standardized_distance_display": sigma,
                },
            )
        )
    return rows


def a0_row(
    *,
    row_id: str,
    quantity: str,
    a0: float,
    measurement: dict[str, Any],
    measurement_display: str,
    dataset: str,
    citation: str,
    public_input_ids: list[str],
    corpus_ids: list[str],
    consts: dict[str, Any],
    note: str,
) -> dict[str, Any]:
    return make_row(
        row_id=row_id,
        section="dark_sector",
        quantity=quantity,
        **{"class": "fitted_comparison_value"},
        oph_value=None,
        oph_value_display="open (dictionary a0 = G n^2 c exact; value not derived)",
        oph_value_exact_rational=None,
        oph_premises=["none: fitted comparison value with no source value"],
        measurement=measurement,
        measurement_display=measurement_display,
        dataset=dataset,
        citation=citation,
        public_input_ids=public_input_ids,
        corpus_source_ids=corpus_ids,
        sigma_distance=None,
        sigma_distance_kind="none_no_source_value",
        verdict="not_evaluable",
        discriminates_from_baseline=False,
        prospective_data="none: no source value to score",
        epistemic_note=note,
        extras={
            "a_dS_m_s2": consts["a_dS_m_s2"],
            "xi_a0_over_a_dS": a0 / consts["a_dS_m_s2"],
            "a0_over_c_H0": a0 / consts["c_H0_m_s2"],
        },
    )


def dark_sector_rows(
    items: dict[str, Any], typed: dict[str, Any], parents: dict[str, Any], consts: dict[str, Any]
) -> list[dict[str, Any]]:
    rows = []
    sparc = parents["sparc_full_rar_calibration"]
    deep = parents["sparc_deep_regime_diagnostic"]
    a0_note = (
        "The OPH dictionary a0 = G n^2 c (collar density n, per-cut variance "
        "slope c) is exact and its value is open; no source value exists to "
        "compare. The ratio xi = a0/a_dS with a_dS = c^2 sqrt(Lambda/3) at the "
        "Planck-central Lambda is the number a source derivation must produce. "
        "The coincidence a0 ~ c H0 is Milgrom's (1983). The fits are same-data "
        "calibration diagnostics, not OPH evidence."
    )
    fits = sparc["fits"]
    paper_display = typed["dark_matter_paper_full_rar_a0_display"]
    for upsilon in ("0.5", "0.4", "0.6"):
        fit = fits[f"disk_mass_to_light_{upsilon}"]
        a0 = float(fit["a0_si_m_per_s2"])
        std = float(fit["galaxy_bootstrap_std_si_m_per_s2"])
        if upsilon == "0.5":
            decimals = int(paper_display["decimals"])
            require(
                round(a0 / 1e-10, decimals) == round(float(paper_display["value"]) / 1e-10, decimals)
                and round(std / 1e-10, decimals) == round(float(paper_display["uncertainty"]) / 1e-10, decimals),
                "full-RAR a0 disagrees with the dark-matter paper display",
            )
        rows.append(
            a0_row(
                row_id=f"a0_sparc_full_rar_upsilon_disk_{upsilon.replace('.', 'p')}",
                quantity=f"a0, SPARC full-RAR fit (Upsilon_disk = {upsilon}, Upsilon_bulge = 0.7)",
                a0=a0,
                measurement={
                    "kind": "fitted_value",
                    "value": a0,
                    "sigma": std,
                    "sigma_kind": "galaxy bootstrap standard deviation",
                    "unit": "m s^-2",
                    "n_galaxies": int(fit["n_galaxies"]),
                    "n_points": int(fit["n_points"]),
                    "scatter_dex": float(fit["scatter_dex"]),
                },
                measurement_display=f"{fmt_a0(a0)} +/- {fmt_a0(std)} m s^-2 (galaxy bootstrap)",
                dataset="SPARC (CDS J/AJ/152/157), unweighted point-level log-residual fit",
                citation="parent receipt code/cosmology/rar_deep_regime/receipts/sparc_full_rar_calibration.json",
                public_input_ids=[],
                corpus_ids=["sparc_full_rar_calibration", "a0_dictionary"] + (["dark_matter_paper_full_rar_a0_display"] if upsilon == "0.5" else []),
                consts=consts,
                note=a0_note,
            )
        )

    mcgaugh = items["mcgaugh2016_a0"]
    literature = sparc["reference_literature"]
    require(
        frac(literature["reported_a0_si_m_per_s2"]) == frac(mcgaugh["value"])
        and frac(literature["reported_random_uncertainty_si_m_per_s2"]) == frac(mcgaugh["uncertainty"])
        and frac(literature["reported_systematic_uncertainty_si_m_per_s2"]) == frac(mcgaugh["systematic_uncertainty"]),
        "McGaugh 2016 values differ between public inputs and the SPARC parent receipt",
    )
    a0_published = float(mcgaugh["value"])
    rows.append(
        a0_row(
            row_id="a0_mcgaugh2016_published",
            quantity="a0, published radial acceleration relation fit",
            a0=a0_published,
            measurement={
                "kind": "published_value",
                "value": a0_published,
                "sigma": float(mcgaugh["uncertainty"]),
                "systematic_uncertainty": float(mcgaugh["systematic_uncertainty"]),
                "unit": "m s^-2",
            },
            measurement_display=(
                f"{fmt_a0(a0_published)} +/- {fmt_a0(float(mcgaugh['uncertainty']))} (random) "
                f"+/- {fmt_a0(float(mcgaugh['systematic_uncertainty']))} (systematic) m s^-2"
            ),
            dataset=mcgaugh["dataset"],
            citation=mcgaugh["citation"],
            public_input_ids=["mcgaugh2016_a0"],
            corpus_ids=["a0_dictionary"],
            consts=consts,
            note=a0_note,
        )
    )

    sweep = deep["deep_radial_acceleration_fraction_sweep"]
    by_fraction = {str(entry["deep_fraction_of_reference_a0"]): entry for entry in sweep}
    require(set(by_fraction) == {"0.3", "0.1", "0.03"}, "unexpected deep-regime fraction sweep")
    deep_display = typed["dark_matter_paper_deep_regime_a0_display"]
    for index, fraction_key in enumerate(("0.3", "0.1", "0.03")):
        entry = by_fraction[fraction_key]
        a0 = float(entry["a0_fixed_exponent_m_s2"])
        interval = [float(v) for v in entry["a0_galaxy_cluster_bootstrap_95pct_m_s2"]]
        decimals = int(deep_display["decimals"])
        require(
            round(a0 / 1e-10, decimals) == round(float(deep_display["values"][index]) / 1e-10, decimals),
            f"deep-regime a0 at f = {fraction_key} disagrees with the dark-matter paper display",
        )
        rows.append(
            a0_row(
                row_id=f"a0_sparc_deep_regime_fixed_exponent_f_{fraction_key.replace('.', 'p')}",
                quantity=f"a0, deep-regime fixed-exponent fit (g_bar < {fraction_key} x 1.2e-10 m s^-2)",
                a0=a0,
                measurement={
                    "kind": "fitted_value",
                    "value": a0,
                    "sigma": None,
                    "interval_95_galaxy_cluster_bootstrap": interval,
                    "unit": "m s^-2",
                    "n_galaxies": int(entry["n_galaxies"]),
                    "n_points": int(entry["n_points"]),
                },
                measurement_display=(
                    f"{fmt_a0(a0)} m s^-2, 95 percent bootstrap interval "
                    f"[{fmt_a0(interval[0])}, {fmt_a0(interval[1])}]"
                ),
                dataset="SPARC deep-regime subset, total model g_obs = g_bar + sqrt(g_bar a0)",
                citation="parent receipt code/cosmology/rar_deep_regime/receipts/sparc_deep_regime_diagnostic.json",
                public_input_ids=[],
                corpus_ids=["sparc_deep_regime_diagnostic", "a0_dictionary", "dark_matter_paper_deep_regime_a0_display"],
                consts=consts,
                note=a0_note,
            )
        )

    btfr = deep["baryonic_tully_fisher"]
    require(frac(btfr["law_exponent"]) == 4, "parent BTFR law exponent is not 4")
    require(frac(typed["btfr_exponent_four"]["value"]) == 4, "typed BTFR exponent is not 4")
    exponent = frac(btfr["free_exponent_M_of_V"])
    exponent_sd = frac(btfr["free_exponent_bootstrap_sd"])
    btfr_display = typed["dark_matter_paper_btfr_exponent_display"]
    require(
        round(float(exponent), int(btfr_display["decimals"])) == float(btfr_display["value"]),
        "BTFR exponent disagrees with the dark-matter paper display",
    )
    sigma_exact = signed_sigma(Fraction(4), exponent, exponent_sd)
    rows.append(
        make_row(
            row_id="btfr_exponent_four_sparc_ols",
            section="dark_sector",
            quantity="baryonic Tully-Fisher exponent in M_b proportional to v_f^x",
            **{"class": "conditional_theorem_postdiction"},
            oph_value=4.0,
            oph_value_display="4 (exact)",
            oph_value_exact_rational="4/1",
            oph_premises=[
                "declared deep-regime scale covariance",
                "quadrature source composition",
                "characterization theorem M_A(r) = r sqrt(M_b a0/G), hence v^4 = G M_b a0",
            ],
            measurement={
                "kind": "gaussian",
                "value": float(exponent),
                "sigma": float(exponent_sd),
                "sigma_kind": "galaxy bootstrap standard deviation of the reciprocal OLS slope",
                "unit": None,
                "n_galaxies": int(btfr["n_galaxies"]),
                "paper_display": {"value": float(btfr_display["value"]), "uncertainty": float(btfr_display["uncertainty"])},
            },
            measurement_display=f"{float(exponent):.3f} +/- {float(exponent_sd):.3f} (paper display 3.75 +/- 0.10)",
            dataset="SPARC, 123 galaxies with positive catalogued V_flat, Upsilon_disk = 0.5, M_b = 0.5 L_[3.6] + 1.33 M_HI",
            citation="parent receipt code/cosmology/rar_deep_regime/receipts/sparc_deep_regime_diagnostic.json",
            public_input_ids=[],
            corpus_source_ids=["btfr_exponent_four", "sparc_deep_regime_diagnostic", "dark_matter_paper_btfr_exponent_display"],
            sigma_distance=float(sigma_exact),
            sigma_distance_kind="signed_one_dimensional",
            verdict=verdict_from_sigma(sigma_exact),
            discriminates_from_baseline=False,
            prospective_data="none: a calibrated errors-in-variables Tully-Fisher likelihood is required before scoring",
            epistemic_note=(
                "The exponent 4 is a conditional characterization theorem (declared "
                "deep-regime scale covariance and quadrature source composition); "
                "neither premise nor a0 is source-derived. The fitted exponent is the "
                "reciprocal of an ordinary least-squares slope of log v_f on log M_b; "
                "ordinary least squares with error in the regressor attenuates that "
                "slope, so the estimator is biased and this is not a calibrated test. "
                "The paper displays 3.75 +/- 0.10; the receipt bootstrap standard "
                "deviation of the exponent is 0.093 and the receipt values are used here."
            ),
        )
    )
    return rows


def not_evaluable_rows(items: dict[str, Any]) -> list[dict[str, Any]]:
    rows = [
        make_row(
            row_id="ringdown_integer_k_comb_fz14",
            section="not_evaluable",
            quantity="black-hole ringdown integer-k comb (FZ-14)",
            **{"class": "not_evaluable"},
            oph_value=None,
            oph_value_display="no comparison contract",
            oph_value_exact_rational=None,
            oph_premises=["none: registered pending owner freeze"],
            measurement={"kind": "none"},
            measurement_display="none",
            dataset="no event likelihood evaluated",
            citation="docs/FROZEN_PREDICTION_LADDER.md, row FZ-14",
            public_input_ids=[],
            corpus_source_ids=["ringdown_fz14_registration"],
            sigma_distance=None,
            sigma_distance_kind="none_not_evaluable",
            verdict="not_evaluable",
            discriminates_from_baseline=False,
            prospective_data="post-anchoring gravitational-wave events after the FZ-14 owner freeze",
            epistemic_note=(
                "FZ-14 is registered pending owner freeze; comparison access is "
                "forbidden until the complete contract is anchored, and no event "
                "likelihood has been evaluated. No comparison contract exists for "
                "this ledger."
            ),
        )
    ]
    neutrino_note = (
        "The corpus supplies no source-derived neutrino mass sum or family "
        "census eligible for comparison (the FZ-01 erratum). The public value is "
        "recorded for completeness; no comparison is made."
    )
    mnu = items["desi_dr2_sum_mnu"]
    rows.append(
        make_row(
            row_id="sum_mnu_desi_dr2_not_evaluable",
            section="not_evaluable",
            quantity="sum of neutrino masses",
            **{"class": "not_evaluable"},
            oph_value=None,
            oph_value_display="none (no corpus value)",
            oph_value_exact_rational=None,
            oph_premises=["none"],
            measurement={
                "kind": "upper_bound",
                "upper_bound": float(mnu["upper_bound"]),
                "confidence_level": float(mnu["confidence_level"]),
                "unit": mnu["unit"],
            },
            measurement_display=f"< {fmt_num(mnu['upper_bound'])} {mnu['unit']} ({fmt_num(frac(mnu['confidence_level']) * 100)} percent CL)",
            dataset=mnu["dataset"],
            citation=mnu["citation"],
            public_input_ids=["desi_dr2_sum_mnu"],
            corpus_source_ids=[],
            sigma_distance=None,
            sigma_distance_kind="none_not_evaluable",
            verdict="not_evaluable",
            discriminates_from_baseline=False,
            prospective_data="none: no corpus value",
            epistemic_note=neutrino_note,
        )
    )
    for item_id in ("planck2018_vi_bao_neff", "act_dr6_p_act_lb_neff"):
        item = items[item_id]
        rows.append(
            make_row(
                row_id=f"neff_{item_id}_not_evaluable",
                section="not_evaluable",
                quantity="N_eff",
                **{"class": "not_evaluable"},
                oph_value=None,
                oph_value_display="none (no corpus value)",
                oph_value_exact_rational=None,
                oph_premises=["none"],
                measurement={
                    "kind": "gaussian",
                    "value": float(item["value"]),
                    "sigma": float(item["uncertainty"]),
                    "unit": None,
                },
                measurement_display=f"{fmt_num(item['value'])} +/- {fmt_num(item['uncertainty'])}",
                dataset=item["dataset"],
                citation=item["citation"],
                public_input_ids=[item_id],
                corpus_source_ids=[],
                sigma_distance=None,
                sigma_distance_kind="none_not_evaluable",
                verdict="not_evaluable",
                discriminates_from_baseline=False,
                prospective_data="none: no corpus value",
                epistemic_note=neutrino_note,
            )
        )
    return rows


def conditional_tilt_attachment(parent: dict[str, Any]) -> dict[str, Any]:
    """Project the pinned overlapping-data diagnostic without changing a row verdict."""
    require(parent["scope"] == "fixed_background_ACT_only_likelihood_slice",
            "ACT likelihood scope mismatch")
    for key in ("p_value", "FZ15_verdict"):
        require(parent[key] is None, f"ACT likelihood cannot assign {key}")
    for key in ("full_cosmological_profile", "native_primordial_generator_identified"):
        require(parent[key] is False, f"ACT likelihood cannot promote {key}")
    fixed, free = (parent["fits"][key] for key in ("fixed_tilt", "free_tilt"))
    delta = float(parent["delta_chi2"])
    values = [float(fit[key]) for fit in (fixed, free) for key in ("ns", "chi2")]
    require(all(math.isfinite(value) for value in values + [delta]),
            "nonfinite ACT likelihood diagnostic")
    require(abs(float(fixed["chi2"]) - float(free["chi2"]) - delta) < 1e-10,
            "ACT likelihood difference mismatch")
    return {
        "artifact_ref": "code/cosmology/tilt_posterior/likelihood_slice.json",
        "scope": parent["scope"],
        "background": {key: float(value) for key, value in parent["background"].items()},
        "calibration_A_act": float(parent["calibration_A_act"]),
        "fixed_tilt": {key: float(fixed[key]) for key in ("ns", "chi2")},
        "free_tilt": {key: float(free[key]) for key in ("ns", "chi2")},
        "delta_chi2": delta,
        "new_independent_data": False,
        "p_value": None,
        "FZ15_verdict": None,
        "full_cosmological_profile": False,
        "native_primordial_generator_identified": False,
    }


def build_receipt(public_bytes: bytes, corpus_bytes: bytes, root: Path = ROOT) -> dict[str, Any]:
    public = load_exact_bytes(public_bytes)
    corpus = load_exact_bytes(corpus_bytes)
    require(public["schema"] == PUBLIC_SCHEMA, "public inputs schema mismatch")
    require(corpus["schema"] == CORPUS_SCHEMA, "corpus inputs schema mismatch")
    items = public["items"]
    typed = corpus["typed"]
    anchor_checks = check_anchors(typed, root)
    parents = load_parents(corpus["parent_receipts"], root)
    consts = derive_constants(items, typed, parents)

    sections = [
        {"id": "primordial", "title": "Primordial", "rows": primordial_rows(items, consts)},
        {
            "id": "background_dark_energy",
            "title": "Background and dark energy",
            "rows": dark_energy_rows(parents),
        },
        {
            "id": "capacity_lambda",
            "title": "Capacity and the cosmological constant",
            "rows": capacity_rows(items, typed, parents, consts),
        },
        {"id": "dark_sector", "title": "Dark sector", "rows": dark_sector_rows(items, typed, parents, consts)},
        {"id": "not_evaluable", "title": "Not evaluable", "rows": not_evaluable_rows(items)},
    ]
    attachment = conditional_tilt_attachment(parents["act_dr6_conditional_likelihood"])
    require(attachment["fixed_tilt"]["ns"] == consts["ns_edge_center_P_C"],
            "ACT likelihood uses a different fixed tilt")
    tilt_row = next(row for row in sections[0]["rows"]
                    if row["row_id"] == "ns_edge_center_P_C_act_dr6_p_act_lb_ns")
    tilt_row["extras"] = {"act_only_conditional_likelihood": attachment}
    tilt_row["epistemic_note"] += (
        " The attached ACT-only likelihood slice shares ACT observations with this "
        "combined posterior. Its fixed-background comparison fits amplitude and "
        "polarization calibration; it changes no uncertainty, verdict, row count, "
        "or frozen prediction and supplies no independent dataset."
    )
    row_ids = [row["row_id"] for section in sections for row in section["rows"]]
    require(len(row_ids) == len(set(row_ids)), "duplicate row id")
    for section in sections:
        for row in section["rows"]:
            require(row["section"] == section["id"], "row placed in the wrong section")
            for item_id in row["public_input_ids"]:
                require(item_id in items, f"unknown public input {item_id}")
            for corpus_id in row["corpus_source_ids"]:
                require(
                    corpus_id in typed or corpus_id in corpus["parent_receipts"],
                    f"unknown corpus source {corpus_id}",
                )

    pins = {
        f"{PACKAGE}/public_inputs.json": sha256_bytes(public_bytes),
        f"{PACKAGE}/corpus_inputs.json": sha256_bytes(corpus_bytes),
    }
    for relative in PINNED_CODE:
        pins[relative] = sha256_bytes((root / relative).read_bytes())
    parent_pins = {}
    for key in corpus["parent_receipts"]:
        entry = corpus["parent_receipts"][key]
        pins[entry["path"]] = entry["sha256"]
        parent_pins[key] = {"path": entry["path"], "sha256": entry["sha256"], "schema": entry["schema"]}

    derived = {key: value for key, value in consts.items() if not key.startswith("_")}
    return {
        "schema": SCHEMA,
        "classification": CLASSIFICATION,
        "promotes_nothing": True,
        "seen_data": True,
        "frozen_prediction_rows": 0,
        "allowed_classes": list(ALLOWED_CLASSES),
        "class_descriptions": CLASS_DESCRIPTIONS,
        "allowed_verdicts": list(ALLOWED_VERDICTS),
        "verdict_rule": VERDICT_RULE,
        "sigma_sign_convention": SIGMA_SIGN_CONVENTION,
        "public_inputs_schema": PUBLIC_SCHEMA,
        "corpus_inputs_schema": CORPUS_SCHEMA,
        "pins": pins,
        "parent_receipts": parent_pins,
        "corpus_anchor_checks": anchor_checks,
        "derived_constants": derived,
        "sections": sections,
        "row_count": len(row_ids),
        "rows_sha256": sha256_bytes(canonical_bytes(sections)),
    }


def _cell(text: str) -> str:
    return str(text).replace("|", "\\|").replace("\n", " ")


def render_markdown(receipt: dict[str, Any], receipt_sha256: str) -> str:
    lines = [
        "# Cosmology postdiction ledger",
        "",
        f"Generated deterministically by `{PACKAGE}/build_cosmology_postdiction_ledger.py`; "
        f"the JSON artifact is `{PACKAGE}/{RECEIPT_NAME}` (sha256 `{receipt_sha256}`, "
        f"rows digest `{receipt['rows_sha256']}`).",
        "",
        "This ledger promotes nothing. Every comparison is on seen data; no row is a "
        "frozen prediction, a score, or evidence for or against OPH. Every number is "
        "read mechanically from `public_inputs.json`, `corpus_inputs.json`, or a pinned "
        "parent receipt. Sigma distances are signed as theory minus measurement over the "
        "quoted one-sigma uncertainty. The fixed verdict rule is: |sigma| < 2 consistent; "
        "2 <= |sigma| < 3 tension; |sigma| >= 3 exceeds_three_sigma_diagnostic; an "
        "upper-bound row is consistent exactly when the theory value lies below the bound; "
        "rows without a comparison contract are not_evaluable. The (w0, wa) rows use the "
        "two-sided normal equivalent of a two-degree-of-freedom Gaussian Mahalanobis distance.",
        "",
        "Row classes:",
        "",
    ]
    for name in receipt["allowed_classes"]:
        lines.append(f"- `{name}`: {receipt['class_descriptions'][name]}")
    lines += [
        "",
        "`discriminates` marks the specific n_s value, zero running, and w = -1 "
        "against specified nonzero alternatives. It does not assert separation "
        "from every member of a model family. The tensor upper limits cannot "
        "distinguish exact zero from arbitrarily small slow-roll signals.",
        "",
        "The [ordered ACT DR6 posterior audit](../tilt_posterior/README.md) "
        "reanalyses the official P-ACT-LB CAMB samples underlying this ledger's "
        "rounded ACT joint summary. It retains both exposed tilt branches, "
        "original multiplicity weights and per-chain burn-in, and reports exact "
        "empirical posterior tail positions and chain sensitivity. It supplies "
        "no additional independent dataset or FZ-15 verdict. The "
        "[public-data correspondence audit](../../../docs/research/PUBLIC_DATA_PREDICTION_AUDIT.md) "
        "distinguishes these conditional cosmological restrictions from a "
        "derivation using native observer records. The scalar tilt is a free "
        "primordial input of the baseline model; fixing it is a conditional "
        "parameter restriction, not a demonstrated change of the baseline dynamics.",
        "",
        "DOE and NSF ended their support for the CMB-S4 project on 9 July 2025 "
        "([CMB-S4](https://cmb-s4.org/)), so the prospective-data column names no "
        "CMB-S4 product. FZ-01's frozen scorer order, which names a CMB-S4 baseline "
        "as its terminal pass, is unchanged; that pass has no decision data unless "
        "such a product is released.",
        "",
    ]
    header = (
        "| quantity | OPH value | premises | measurement (dataset) | sigma | verdict | "
        "discriminates | prospective data |"
    )
    rule = "| --- | --- | --- | --- | --- | --- | --- | --- |"
    for section in receipt["sections"]:
        lines += [f"## {section['title']}", "", header, rule]
        for row in section["rows"]:
            lines.append(
                "| "
                + " | ".join(
                    [
                        _cell(row["quantity"]),
                        _cell(row["oph_value_display"]),
                        _cell("; ".join(row["oph_premises"])),
                        _cell(f"{row['measurement_display']} ({row['dataset']})"),
                        fmt_sigma(row["sigma_distance"]),
                        row["verdict"],
                        "yes" if row["discriminates_from_baseline"] else "no",
                        _cell(row["prospective_data"]),
                    ]
                )
                + " |"
            )
        lines.append("")
        notes = section_notes(section)
        if notes:
            lines += notes + [""]
    lines += [
        "## Reproduce",
        "",
        "```bash",
        f"python3 {PACKAGE}/build_cosmology_postdiction_ledger.py --check",
        f"python3 {PACKAGE}/verify_cosmology_postdiction_ledger_independent.py",
        "```",
        "",
    ]
    return "\n".join(lines)


def section_notes(section: dict[str, Any]) -> list[str]:
    notes: list[str] = []
    if section["id"] == "primordial":
        row = next(row for row in section["rows"]
                   if row["row_id"] == "ns_edge_center_P_C_act_dr6_p_act_lb_ns")
        diagnostic = row["extras"]["act_only_conditional_likelihood"]
        notes.append(
            "Attached to the P_C ACT joint-posterior row: the "
            "[conditional ACT-only likelihood](../tilt_posterior/README.md#conditional-act-only-likelihood-slice) "
            f"gives chi^2 {diagnostic['fixed_tilt']['chi2']:.9f} for the fixed tilt "
            f"and {diagnostic['free_tilt']['chi2']:.9f} for free tilt "
            f"(fitted n_s {diagnostic['free_tilt']['ns']:.9f}); "
            f"the difference is {diagnostic['delta_chi2']:.9f}. "
            "The background and absolute calibration are fixed; amplitude and "
            "polarization calibration are fitted under both branches. These ACT "
            "observations overlap the combined posterior. This attachment adds no "
            "independent dataset or row, changes no quoted uncertainty or verdict, "
            "and supplies neither a full cosmological/nuisance profile nor a p-value "
            "or FZ-15 adjudication. The pinned parent retains the complete data, "
            "covariance and supplied radiation-transfer assumptions."
        )
    elif section["id"] == "background_dark_energy":
        notes.append("Chain diagnostics read from the parent DESI receipt:")
        notes.append("")
        for row in section["rows"]:
            extras = row["extras"]
            notes.append(
                f"- {row['dataset'].split(' (')[0]}: Mahalanobis squared "
                f"{extras['mahalanobis_squared']:.3f}; component sigma w0 "
                f"{extras['component_sigma_distances']['w0']:+.2f}, wa "
                f"{extras['component_sigma_distances']['wa']:+.2f}; posterior mass with "
                f"w(a) >= -1 on 0 <= z <= 2: {extras['posterior_mass_w_ge_minus_one_for_0_le_z_le_2']:.2e}; "
                f"mass with w0 > -1: {extras['posterior_mass_w0_gt_minus_one']:.4f}; "
                f"mass with wa >= 0: {extras['posterior_mass_wa_nonnegative']:.2e}."
            )
    elif section["id"] == "capacity_lambda":
        notes.append("Percent residuals, candidate minus comparison over comparison:")
        notes.append("")
        for row in section["rows"]:
            extras = row["extras"]
            line = (
                f"- {row['quantity']} against {row['dataset']}: "
                f"{extras['percent_residual_in_lambda']:+.2f} percent in the constant, "
                f"{extras['percent_residual_in_capacity']:+.2f} percent in the capacity coordinate"
            )
            if "standardized_distance_display" in extras:
                line += (
                    f"; standardized display distance {extras['standardized_distance_display']:+.2f}"
                    " (no verdict)"
                )
            notes.append(line + ".")
        low, high = FZ01_CAPACITY_BAND
        notes.append("")
        notes.append(
            f"The frozen FZ-01 capacity band is N in [{low:.5e}, {high:.5e}], the"
            " registered target; the display candidates above are not that target."
            " Its registered policy scores ln N against the band edges on the frozen"
            " Planck 2018 posterior and issues no verdict until sigma_lnN <= 0.01, so"
            " no FZ-01 verdict exists at today's precision"
            " ([ladder](../../../docs/FROZEN_PREDICTION_LADDER.md#fz-01))."
        )
    elif section["id"] == "dark_sector":
        notes.append(
            "Dimensionless ratios with a_dS = c^2 sqrt(Lambda/3) at the Planck-central "
            "Lambda and c H0 at H0 = 67.36 km/s/Mpc:"
        )
        notes.append("")
        for row in section["rows"]:
            if row["class"] != "fitted_comparison_value":
                continue
            extras = row["extras"]
            notes.append(
                f"- {row['quantity']}: xi = a0/a_dS = {extras['xi_a0_over_a_dS']:.4f}; "
                f"a0/(c H0) = {extras['a0_over_c_H0']:.4f}."
            )
    return notes


def build_outputs(root: Path = ROOT) -> tuple[bytes, bytes]:
    receipt = build_receipt(
        (root / PACKAGE / "public_inputs.json").read_bytes(),
        (root / PACKAGE / "corpus_inputs.json").read_bytes(),
        root,
    )
    receipt_bytes = render_bytes(receipt)
    markdown = render_markdown(receipt, sha256_bytes(receipt_bytes)).encode("utf-8")
    return receipt_bytes, markdown


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output-dir", type=Path, default=HERE)
    parser.add_argument(
        "--check",
        action="store_true",
        help="require the committed receipt and Markdown to be byte-identical to a rebuild",
    )
    args = parser.parse_args(argv)
    receipt_bytes, markdown = build_outputs()
    receipt_path = args.output_dir / RECEIPT_NAME
    markdown_path = args.output_dir / MARKDOWN_NAME
    if args.check:
        current = (
            receipt_path.is_file()
            and markdown_path.is_file()
            and receipt_path.read_bytes() == receipt_bytes
            and markdown_path.read_bytes() == markdown
        )
        if not current:
            print("FAIL: committed cosmology postdiction ledger is stale")
            return 1
        print(f"PASS: {RECEIPT_NAME} and {MARKDOWN_NAME} are current")
        return 0
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.write_bytes(receipt_bytes)
    markdown_path.write_bytes(markdown)
    print(f"wrote {receipt_path.name} ({sha256_bytes(receipt_bytes)}) and {markdown_path.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
