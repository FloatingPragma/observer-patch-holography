#!/usr/bin/env python3
"""Independent verifier for the cosmology postdiction ledger receipt.

This script imports neither the producer nor any of its helpers.  It re-reads
the two input files and the pinned parent receipts, recomputes every derived
number (exact rationals on the typed decimal digits where the theory value is
rational, floats in a fixed operation order otherwise), checks every row
field against that recomputation, checks every pin and every corpus anchor,
scans all prose for the banned progress vocabulary, and exits 0 only when the
committed receipt and its Markdown rendering agree with the replay.

What is not proved here.  Replay equality certifies that the committed
receipt is exactly the declared computation on the committed inputs.  It does
not make any row a prediction, a score, or evidence for or against OPH.
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
from statistics import NormalDist
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import strict_json  # noqa: E402

PACKAGE = "code/cosmology/postdiction_ledger"
RECEIPT_PATH = HERE / "runtime" / "cosmology_postdiction_ledger.json"
MARKDOWN_PATH = HERE / "COSMOLOGY_POSTDICTION_LEDGER.md"
PUBLIC_PATH = HERE / "public_inputs.json"
CORPUS_PATH = HERE / "corpus_inputs.json"

SCHEMA = "oph.cosmology.postdiction_ledger.v1"
PUBLIC_SCHEMA = "oph.cosmology.postdiction_ledger.public_inputs.v1"
CORPUS_SCHEMA = "oph.cosmology.postdiction_ledger.corpus_inputs.v1"
CLASSIFICATION = (
    "retrospective postdiction ledger; every comparison is on seen data; "
    "no row is a frozen prediction or a score"
)
CLASSES = {
    "conditional_theorem_postdiction",
    "shared_baseline",
    "closure_candidate_display",
    "fitted_comparison_value",
    "not_evaluable",
}
VERDICTS = {
    "consistent",
    "tension",
    "exceeds_three_sigma_diagnostic",
    "exceeds_stated_confidence_bound",
    "not_evaluable",
}
TOP_LEVEL_KEYS = {
    "schema",
    "classification",
    "promotes_nothing",
    "seen_data",
    "frozen_prediction_rows",
    "allowed_classes",
    "class_descriptions",
    "allowed_verdicts",
    "verdict_rule",
    "sigma_sign_convention",
    "public_inputs_schema",
    "corpus_inputs_schema",
    "pins",
    "parent_receipts",
    "corpus_anchor_checks",
    "derived_constants",
    "sections",
    "row_count",
    "rows_sha256",
}
ROW_KEYS = {
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
}
PINNED_CODE = (
    f"{PACKAGE}/build_cosmology_postdiction_ledger.py",
    f"{PACKAGE}/verify_cosmology_postdiction_ledger_independent.py",
    f"{PACKAGE}/test_cosmology_postdiction_ledger.py",
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
BANNED_PROSE = [
    (re.compile(r"\bnow\b", re.IGNORECASE), "now"),
    (re.compile(r"\balready\b", re.IGNORECASE), "already"),
    (re.compile(r"\bstill\b", re.IGNORECASE), "still"),
    (re.compile(r"\bcurrently\b", re.IGNORECASE), "currently"),
    (re.compile(r"\brecently\b", re.IGNORECASE), "recently"),
    (re.compile(r"\bpreviously\b", re.IGNORECASE), "previously"),
    (re.compile(r"\bno longer\b", re.IGNORECASE), "no longer"),
    (re.compile(r"\bnot yet\b", re.IGNORECASE), "not yet"),
    (re.compile(r"\bremains open\b", re.IGNORECASE), "remains open"),
    (re.compile(r"\bso far\b", re.IGNORECASE), "so far"),
    (re.compile(r"\bfuture work\b", re.IGNORECASE), "future work"),
    (re.compile(r"\bnext step\b", re.IGNORECASE), "next step"),
    (re.compile(r"—|---"), "em dash"),
    (re.compile(r"\bnot\b[^.!?\n]{0,120}\bbut\b", re.IGNORECASE), "not-X-but-Y"),
    (re.compile(r"\bhonest", re.IGNORECASE), "h-word"),
]
DEVELOPER_PATH = re.compile(r"/Users/[^/]+/|/home/[^/]+/|[A-Za-z]:[\\/]Users[\\/]")
TIMESTAMP = re.compile(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}")


class VerificationError(ValueError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def load_exact(path: Path) -> Any:
    text = path.read_text(encoding="utf-8")
    strict_json.loads(text)
    return json.loads(text, parse_float=Decimal)


def fr(value: Any) -> Fraction:
    require(isinstance(value, (int, Decimal)) and not isinstance(value, bool), "exact number required")
    return Fraction(value)


def fl(value: Any) -> float:
    require(isinstance(value, (int, Decimal)) and not isinstance(value, bool), "number required")
    return float(value)


def norm(text: str) -> str:
    return re.sub(r"\s+", " ", text)


def rule(sigma: Fraction | float) -> str:
    magnitude = abs(sigma)
    if magnitude < 2:
        return "consistent"
    if magnitude < 3:
        return "tension"
    return "exceeds_three_sigma_diagnostic"


def same_float(actual: Any, expected: float, label: str) -> None:
    require(isinstance(actual, (int, Decimal)) and not isinstance(actual, bool), f"{label}: number required")
    require(float(actual) == expected, f"{label}: {float(actual)!r} != {expected!r}")


def scan_prose(value: Any, where: str, problems: list[str]) -> None:
    if isinstance(value, dict):
        for key, item in value.items():
            scan_prose(item, f"{where}/{key}", problems)
    elif isinstance(value, list):
        for index, item in enumerate(value):
            scan_prose(item, f"{where}[{index}]", problems)
    elif isinstance(value, str):
        for pattern, label in BANNED_PROSE:
            if pattern.search(value):
                problems.append(f"{where}: banned prose ({label})")
        if DEVELOPER_PATH.search(value):
            problems.append(f"{where}: developer home path")
        if TIMESTAMP.search(value):
            problems.append(f"{where}: timestamp")


class Expected:
    """Recompute every hard row field from the inputs and parent receipts."""

    def __init__(self, public: dict[str, Any], corpus: dict[str, Any], parents: dict[str, Any]) -> None:
        self.items = public["items"]
        self.typed = corpus["typed"]
        self.parents = parents
        self.rows: dict[str, dict[str, Any]] = {}
        self.order: list[tuple[str, str]] = []
        self.constants()
        self.primordial()
        self.dark_energy()
        self.capacity()
        self.dark_sector()
        self.not_evaluable()

    def add(self, section: str, row: dict[str, Any]) -> None:
        row["section"] = section
        require(row["row_id"] not in self.rows, "duplicate expected row")
        self.rows[row["row_id"]] = row
        self.order.append((section, row["row_id"]))

    def constants(self) -> None:
        items, typed = self.items, self.typed
        c = fl(items["si_speed_of_light"]["value"])
        mpc = fl(items["iau2015_megaparsec"]["value"])
        lp = fl(items["codata2022_planck_length"]["value"])
        h0 = fl(items["planck2018_vi_h0"]["value"])
        omega_lambda = fl(items["planck2018_vi_omega_lambda"]["value"])
        h0_si = h0 * 1000.0 / mpc
        lambda_planck = 3.0 * omega_lambda * (h0_si / c) ** 2
        lambda_lp2_planck = lambda_planck * lp**2
        p_c = fr(typed["P_C"]["value"])
        p_fwd = fr(typed["P_fwd"]["value"])
        phi = (1.0 + math.sqrt(5.0)) / 2.0
        desi = self.parents["official_desi_dr2_fz13_retrospective"]["base_lcdm_capacity_display"]["combined"]
        desi_mean = fl(desi["Lambda_lP2"]["weighted_mean"])
        desi_std = fl(desi["Lambda_lP2"]["weighted_std"])
        self.exact = {"P_C": p_c, "P_fwd": p_fwd, "ns_P_C": 1 - p_c / 48, "ns_P_fwd": 1 - p_fwd / 48}
        self.derived = {
            "speed_of_light_m_s": c,
            "megaparsec_m": mpc,
            "planck_length_m": lp,
            "planck2018_H0_km_s_Mpc": h0,
            "planck2018_Omega_Lambda": omega_lambda,
            "planck2018_H0_s_minus_1": h0_si,
            "Lambda_planck_central_m_minus_2": lambda_planck,
            "Lambda_lP2_planck_central": lambda_lp2_planck,
            "capacity_coordinate_planck_central": 3.0 * math.pi / lambda_lp2_planck,
            "a_dS_m_s2": c**2 * math.sqrt(lambda_planck / 3.0),
            "c_H0_m_s2": c * h0_si,
            "P_C": float(p_c),
            "P_fwd": float(p_fwd),
            "golden_ratio": phi,
            "euler_number": math.e,
            "ns_edge_center_P_C": float(self.exact["ns_P_C"]),
            "ns_edge_center_P_fwd": float(self.exact["ns_P_fwd"]),
            "ns_edge_center_branch_spread": float(self.exact["ns_P_fwd"] - self.exact["ns_P_C"]),
            "ns_clock_branch_P_C": 1.0 - math.e * (float(p_c) - phi),
            "Lambda_lP2_desi_dr2_display_mean": desi_mean,
            "Lambda_lP2_desi_dr2_display_std": desi_std,
            "capacity_coordinate_desi_dr2_display": 3.0 * math.pi / desi_mean,
        }
        display = typed["flagship_lambda_planck_central_display"]
        rounded = float("%.*e" % (int(display["significant_digits"]) - 1, lambda_planck))
        require(rounded == fl(display["value"]), "flagship Lambda display cross-check")
        require(round(float(self.exact["ns_P_C"]), 5) == fl(typed["edge_center_paper_display"]["value"]), "n_s paper display")
        require(round(self.derived["ns_clock_branch_P_C"], 5) == fl(typed["clock_branch_paper_display"]["value"]), "clock display")

    def gaussian(self, row_id: str, cls: str, theory: Fraction | float, item_id: str, discriminates: bool, prospective: str) -> dict[str, Any]:
        item = self.items[item_id]
        value, unc = fr(item["value"]), fr(item["uncertainty"])
        if isinstance(theory, Fraction):
            sigma_exact = (theory - value) / unc
            sigma, verdict = float(sigma_exact), rule(sigma_exact)
            oph, exact = float(theory), f"{theory.numerator}/{theory.denominator}"
        else:
            sigma = (theory - float(value)) / float(unc)
            verdict, oph, exact = rule(sigma), theory, None
        return {
            "row_id": row_id,
            "class": cls,
            "oph_value": oph,
            "oph_value_exact_rational": exact,
            "measurement": {"kind": "gaussian", "value": float(value), "sigma": float(unc)},
            "public_input_ids": [item_id],
            "sigma_distance": sigma,
            "sigma_distance_kind": "signed_one_dimensional",
            "verdict": verdict,
            "discriminates_from_baseline": discriminates,
            "prospective_data": prospective,
            "dataset": item["dataset"],
            "citation": item["citation"],
        }

    def bound(self, row_id: str, cls: str, item_id: str, discriminates: bool, prospective: str) -> dict[str, Any]:
        item = self.items[item_id]
        bound = fr(item["upper_bound"])
        return {
            "row_id": row_id,
            "class": cls,
            "oph_value": 0.0,
            "oph_value_exact_rational": "0/1",
            "measurement": {
                "kind": "upper_bound",
                "upper_bound": float(bound),
                "confidence_level": float(fr(item["confidence_level"])),
            },
            "public_input_ids": [item_id],
            "sigma_distance": None,
            "sigma_distance_kind": "none_upper_bound",
            "verdict": "consistent" if Fraction(0) < bound else "exceeds_stated_confidence_bound",
            "discriminates_from_baseline": discriminates,
            "prospective_data": prospective,
            "dataset": item["dataset"],
            "citation": item["citation"],
        }

    def primordial(self) -> None:
        cmb = "Simons Observatory"
        for branch in ("P_C", "P_fwd"):
            for item_id in NS_DATASETS:
                self.add("primordial", self.gaussian(
                    f"ns_edge_center_{branch}_{item_id}", "conditional_theorem_postdiction",
                    self.exact[f"ns_{branch}"], item_id, True, cmb))
        parent = self.parents["act_dr6_conditional_likelihood"]
        require(parent["scope"] == "fixed_background_ACT_only_likelihood_slice", "ACT slice scope")
        require(parent["p_value"] is None and parent["FZ15_verdict"] is None,
                "ACT slice cannot assign a significance or verdict")
        require(parent["full_cosmological_profile"] is False
                and parent["native_primordial_generator_identified"] is False,
                "ACT slice cannot promote physical or profiling scope")
        fixed, free = parent["fits"]["fixed_tilt"], parent["fits"]["free_tilt"]
        same_float(fixed["ns"], float(self.exact["ns_P_C"]), "ACT fixed branch")
        same_float(parent["delta_chi2"], fl(fixed["chi2"]) - fl(free["chi2"]), "ACT delta chi2")
        diagnostic = {
            "artifact_ref": "code/cosmology/tilt_posterior/likelihood_slice.json",
            "scope": parent["scope"],
            "background": {key: fl(value) for key, value in parent["background"].items()},
            "calibration_A_act": fl(parent["calibration_A_act"]),
            "fixed_tilt": {key: fl(fixed[key]) for key in ("ns", "chi2")},
            "free_tilt": {key: fl(free[key]) for key in ("ns", "chi2")},
            "delta_chi2": fl(parent["delta_chi2"]),
            "new_independent_data": False,
            "p_value": None,
            "FZ15_verdict": None,
            "full_cosmological_profile": False,
            "native_primordial_generator_identified": False,
        }
        self.rows["ns_edge_center_P_C_act_dr6_p_act_lb_ns"]["extras"] = {
            "act_only_conditional_likelihood": diagnostic}
        for item_id in NS_DATASETS:
            self.add("primordial", self.gaussian(
                f"ns_clock_branch_P_C_{item_id}", "conditional_theorem_postdiction",
                self.derived["ns_clock_branch_P_C"], item_id, True, cmb))
        for item_id in ("planck2018_running", "act_dr6_p_act_lb_running"):
            self.add("primordial", self.gaussian(
                f"running_zero_{item_id}", "conditional_theorem_postdiction", Fraction(0), item_id, True, cmb))
        for item_id in ("bk18_r_0p05", "tristram2022_r", "planck2018_x_r_0p002_bk15"):
            row = self.bound(
                f"tensor_zero_{item_id}", "conditional_theorem_postdiction", item_id, False,
                "LiteBIRD (tests specified tensor alternatives above the achieved bound)")
            # Bind the scientific premises, not only the zero and upper bound.
            # No source term in a second-order tensor equation does not fix
            # its two homogeneous initial data.
            row["oph_premises"] = [
                "linear tensor evolution on the declared background",
                "zero initial tensor amplitude and velocity",
                "no tensor forcing",
                "positive scalar power at the comparison pivot",
            ]
            row["epistemic_note"] = (
                "Zero linear primordial tensor amplitude requires zero initial tensor "
                "amplitude and velocity, no tensor forcing, and positive scalar power. "
                "A scalar source alone leaves homogeneous tensor modes unconstrained. "
                "These extra initial conditions are branch premises, not source-derived "
                "results. Slow roll has no universal positive lower bound on r; an upper "
                "limit tests only specified alternatives above it. Observable comparisons "
                "require separate treatment of scalar-induced tensors, lensing and "
                "foregrounds. Bound rows carry no sigma distance. Seen data."
            )
            self.add("primordial", row)
        for item_id in (
            "planck2018_x_cdi_axion_i_100beta_klow",
            "planck2018_x_cdi_axion_i_100beta_kmid",
            "planck2018_x_cdi_axion_i_100beta_khigh",
            "planck2018_x_cdi_general_100beta_klow",
        ):
            self.add("primordial", self.bound(
                f"isocurvature_zero_{item_id}", "shared_baseline", item_id, False,
                "LiteBIRD (bound tightening only; zero is the shared baseline)"))
        for item_id in ("planck2018_ix_fnl_local", "planck2018_ix_fnl_equil", "planck2018_ix_fnl_ortho"):
            self.add("primordial", self.gaussian(
                f"fnl_zero_{item_id}", "shared_baseline", Fraction(0), item_id, False,
                "none scoring; large-scale-structure f_NL programmes tighten a shared-baseline bound"))
        for item_id in ("planck2018_vi_omegak", "planck2018_vi_omegak_bao", "desi_dr2_cmb_omegak"):
            self.add("primordial", self.gaussian(
                f"curvature_zero_{item_id}", "shared_baseline", Fraction(0), item_id, False, "DESI DR3; Euclid"))

    def dark_energy(self) -> None:
        desi = self.parents["official_desi_dr2_fz13_retrospective"]
        for dataset_key, slug in DESI_COMBINATIONS:
            block = desi["datasets"][dataset_key]
            combined = block["combined"]
            diagnostic = block["fixed_capacity_point_gaussian_diagnostic"]
            w0_mean, w0_std = fr(combined["w0_mean"]), fr(combined["w0_std"])
            wa_mean, wa_std = fr(combined["wa_mean"]), fr(combined["wa_std"])
            cov = fl(combined["w0_wa_covariance"])
            # Cross-check the parent's Gaussian diagnostic from its own moments.
            v0, va = float(w0_std) ** 2, float(wa_std) ** 2
            det = v0 * va - cov * cov
            d0, da = float(w0_mean) + 1.0, float(wa_mean)
            q = (va * d0 * d0 - 2.0 * cov * d0 * da + v0 * da * da) / det
            survival = math.exp(-0.5 * q)
            sigma_equivalent = NormalDist().inv_cdf(1.0 - survival / 2.0)
            require(abs(q - fl(diagnostic["mahalanobis_squared"])) <= 1e-9 * q, f"{dataset_key}: Mahalanobis replay")
            require(
                abs(sigma_equivalent - fl(diagnostic["two_sided_normal_sigma_equivalent"])) <= 1e-9 * sigma_equivalent,
                f"{dataset_key}: sigma-equivalent replay",
            )
            sigma = fl(diagnostic["two_sided_normal_sigma_equivalent"])
            self.add("background_dark_energy", {
                "row_id": f"w0wa_fixed_capacity_{slug}",
                "class": "conditional_theorem_postdiction",
                "oph_value": {"w0": -1.0, "wa": 0.0},
                "oph_value_exact_rational": None,
                "measurement": {
                    "kind": "two_parameter_gaussian_moments",
                    "w0_mean": float(w0_mean),
                    "w0_std": float(w0_std),
                    "wa_mean": float(wa_mean),
                    "wa_std": float(wa_std),
                    "w0_wa_covariance": cov,
                    "w0_wa_correlation": fl(combined["w0_wa_correlation"]),
                },
                "public_input_ids": [],
                "sigma_distance": sigma,
                "sigma_distance_kind": "two_dof_gaussian_mahalanobis_two_sided_normal_equivalent",
                "verdict": rule(sigma),
                "discriminates_from_baseline": True,
                "prospective_data": "DESI DR3; Euclid",
                "extras": {
                    "mahalanobis_squared": fl(diagnostic["mahalanobis_squared"]),
                    "chi2_2dof_survival": fl(diagnostic["chi2_2dof_survival"]),
                    "component_sigma_distances": {
                        "w0": float((Fraction(-1) - w0_mean) / w0_std),
                        "wa": float((Fraction(0) - wa_mean) / wa_std),
                    },
                    "posterior_mass_w_ge_minus_one_for_0_le_z_le_2": fl(combined["posterior_mass_w_ge_minus_one_for_0_le_z_le_2"]),
                    "posterior_mass_w0_gt_minus_one": fl(combined["posterior_mass_w0_gt_minus_one"]),
                    "posterior_mass_wa_nonnegative": fl(combined["posterior_mass_wa_nonnegative"]),
                },
            })

    def capacity(self) -> None:
        d = self.derived
        lp = d["planck_length_m"]
        percent_display = [fl(v) for v in self.typed["capacity_lambda_percent_display"]["values"]]
        for index, candidate_id in enumerate(("capacity_candidate_a", "capacity_candidate_b")):
            n = fl(self.typed[candidate_id]["value"])
            lambda_candidate = 3.0 * math.pi / (n * lp**2)
            lambda_lp2_candidate = 3.0 * math.pi / n
            pct_lambda = 100.0 * (lambda_candidate - d["Lambda_planck_central_m_minus_2"]) / d["Lambda_planck_central_m_minus_2"]
            require(round(pct_lambda, 2) == percent_display[index], f"{candidate_id}: percent display")
            self.add("capacity_lambda", {
                "row_id": f"lambda_{candidate_id}_planck2018_central",
                "class": "closure_candidate_display",
                "oph_value": lambda_candidate,
                "oph_value_exact_rational": None,
                "measurement": {
                    "kind": "central_value_display",
                    "value": d["Lambda_planck_central_m_minus_2"],
                    "sigma": None,
                    "capacity_coordinate": d["capacity_coordinate_planck_central"],
                },
                "public_input_ids": ["planck2018_vi_h0", "planck2018_vi_omega_lambda", "codata2022_planck_length", "si_speed_of_light", "iau2015_megaparsec"],
                "sigma_distance": None,
                "sigma_distance_kind": "none_central_value_display",
                "verdict": "not_evaluable",
                "discriminates_from_baseline": False,
                "prospective_data": "DESI DR3; Euclid (background posterior displays only)",
                "extras": {
                    "percent_residual_in_lambda": pct_lambda,
                    "percent_residual_in_capacity": 100.0 * (n - d["capacity_coordinate_planck_central"]) / d["capacity_coordinate_planck_central"],
                    "planck_chain_capacity_coordinate_display": fl(self.typed["planck_chain_capacity_coordinate_display"]["value"]),
                },
            })
            mean, std = d["Lambda_lP2_desi_dr2_display_mean"], d["Lambda_lP2_desi_dr2_display_std"]
            sigma = (lambda_lp2_candidate - mean) / std
            self.add("capacity_lambda", {
                "row_id": f"lambda_lp2_{candidate_id}_desi_dr2_display",
                "class": "closure_candidate_display",
                "oph_value": lambda_lp2_candidate,
                "oph_value_exact_rational": None,
                "measurement": {
                    "kind": "gaussian",
                    "value": mean,
                    "sigma": std,
                    "capacity_coordinate": d["capacity_coordinate_desi_dr2_display"],
                },
                "public_input_ids": [],
                "sigma_distance": None,
                "sigma_distance_kind": "none_not_evaluable",
                "verdict": "not_evaluable",
                "discriminates_from_baseline": False,
                "prospective_data": "DESI DR3; Euclid (background posterior displays only)",
                "extras": {
                    "percent_residual_in_lambda": 100.0 * (lambda_lp2_candidate - mean) / mean,
                    "percent_residual_in_capacity": 100.0 * (n - d["capacity_coordinate_desi_dr2_display"]) / d["capacity_coordinate_desi_dr2_display"],
                    "standardized_distance_display": sigma,
                },
            })

    def a0(self, row_id: str, a0: float, measurement: dict[str, Any], public_ids: list[str]) -> dict[str, Any]:
        d = self.derived
        return {
            "row_id": row_id,
            "class": "fitted_comparison_value",
            "oph_value": None,
            "oph_value_exact_rational": None,
            "measurement": measurement,
            "public_input_ids": public_ids,
            "sigma_distance": None,
            "sigma_distance_kind": "none_no_source_value",
            "verdict": "not_evaluable",
            "discriminates_from_baseline": False,
            "prospective_data": "none: no source value to score",
            "extras": {
                "a_dS_m_s2": d["a_dS_m_s2"],
                "xi_a0_over_a_dS": a0 / d["a_dS_m_s2"],
                "a0_over_c_H0": a0 / d["c_H0_m_s2"],
            },
        }

    def dark_sector(self) -> None:
        sparc = self.parents["sparc_full_rar_calibration"]
        deep = self.parents["sparc_deep_regime_diagnostic"]
        for upsilon in ("0.5", "0.4", "0.6"):
            fit = sparc["fits"][f"disk_mass_to_light_{upsilon}"]
            a0 = fl(fit["a0_si_m_per_s2"])
            self.add("dark_sector", self.a0(
                f"a0_sparc_full_rar_upsilon_disk_{upsilon.replace('.', 'p')}", a0,
                {"kind": "fitted_value", "value": a0, "sigma": fl(fit["galaxy_bootstrap_std_si_m_per_s2"]),
                 "n_galaxies": int(fit["n_galaxies"]), "n_points": int(fit["n_points"]), "scatter_dex": fl(fit["scatter_dex"])},
                []))
        paper = self.typed["dark_matter_paper_full_rar_a0_display"]
        primary = sparc["fits"]["disk_mass_to_light_0.5"]
        require(round(fl(primary["a0_si_m_per_s2"]) / 1e-10, 4) == round(fl(paper["value"]) / 1e-10, 4), "full-RAR paper display")
        mcgaugh = self.items["mcgaugh2016_a0"]
        literature = sparc["reference_literature"]
        require(fr(literature["reported_a0_si_m_per_s2"]) == fr(mcgaugh["value"]), "McGaugh value cross-check")
        require(fr(literature["reported_random_uncertainty_si_m_per_s2"]) == fr(mcgaugh["uncertainty"]), "McGaugh random")
        require(fr(literature["reported_systematic_uncertainty_si_m_per_s2"]) == fr(mcgaugh["systematic_uncertainty"]), "McGaugh systematic")
        a0_pub = fl(mcgaugh["value"])
        self.add("dark_sector", self.a0(
            "a0_mcgaugh2016_published", a0_pub,
            {"kind": "published_value", "value": a0_pub, "sigma": fl(mcgaugh["uncertainty"]),
             "systematic_uncertainty": fl(mcgaugh["systematic_uncertainty"])},
            ["mcgaugh2016_a0"]))
        sweep = {str(entry["deep_fraction_of_reference_a0"]): entry for entry in deep["deep_radial_acceleration_fraction_sweep"]}
        require(set(sweep) == {"0.3", "0.1", "0.03"}, "fraction sweep keys")
        displays = [fl(v) for v in self.typed["dark_matter_paper_deep_regime_a0_display"]["values"]]
        for index, key in enumerate(("0.3", "0.1", "0.03")):
            entry = sweep[key]
            a0 = fl(entry["a0_fixed_exponent_m_s2"])
            require(round(a0 / 1e-10, 3) == round(displays[index] / 1e-10, 3), f"deep-regime display f = {key}")
            self.add("dark_sector", self.a0(
                f"a0_sparc_deep_regime_fixed_exponent_f_{key.replace('.', 'p')}", a0,
                {"kind": "fitted_value", "value": a0, "sigma": None,
                 "interval_95_galaxy_cluster_bootstrap": [fl(v) for v in entry["a0_galaxy_cluster_bootstrap_95pct_m_s2"]],
                 "n_galaxies": int(entry["n_galaxies"]), "n_points": int(entry["n_points"])},
                []))
        btfr = deep["baryonic_tully_fisher"]
        require(fr(btfr["law_exponent"]) == 4 and fr(self.typed["btfr_exponent_four"]["value"]) == 4, "BTFR exponent four")
        exponent, sd = fr(btfr["free_exponent_M_of_V"]), fr(btfr["free_exponent_bootstrap_sd"])
        require(round(float(exponent), 2) == fl(self.typed["dark_matter_paper_btfr_exponent_display"]["value"]), "BTFR paper display")
        sigma_exact = (Fraction(4) - exponent) / sd
        self.add("dark_sector", {
            "row_id": "btfr_exponent_four_sparc_ols",
            "class": "conditional_theorem_postdiction",
            "oph_value": 4.0,
            "oph_value_exact_rational": "4/1",
            "measurement": {"kind": "gaussian", "value": float(exponent), "sigma": float(sd), "n_galaxies": int(btfr["n_galaxies"])},
            "public_input_ids": [],
            "sigma_distance": float(sigma_exact),
            "sigma_distance_kind": "signed_one_dimensional",
            "verdict": rule(sigma_exact),
            "discriminates_from_baseline": False,
            "prospective_data": "none: a calibrated errors-in-variables Tully-Fisher likelihood is required before scoring",
        })

    def not_evaluable(self) -> None:
        self.add("not_evaluable", {
            "row_id": "ringdown_integer_k_comb_fz14",
            "class": "not_evaluable",
            "oph_value": None,
            "oph_value_exact_rational": None,
            "measurement": {"kind": "none"},
            "public_input_ids": [],
            "sigma_distance": None,
            "sigma_distance_kind": "none_not_evaluable",
            "verdict": "not_evaluable",
            "discriminates_from_baseline": False,
            "prospective_data": "post-anchoring gravitational-wave events after the FZ-14 owner freeze",
        })
        mnu = self.items["desi_dr2_sum_mnu"]
        self.add("not_evaluable", {
            "row_id": "sum_mnu_desi_dr2_not_evaluable",
            "class": "not_evaluable",
            "oph_value": None,
            "oph_value_exact_rational": None,
            "measurement": {"kind": "upper_bound", "upper_bound": fl(mnu["upper_bound"]), "confidence_level": fl(mnu["confidence_level"])},
            "public_input_ids": ["desi_dr2_sum_mnu"],
            "sigma_distance": None,
            "sigma_distance_kind": "none_not_evaluable",
            "verdict": "not_evaluable",
            "discriminates_from_baseline": False,
            "prospective_data": "none: no corpus value",
            "dataset": mnu["dataset"],
            "citation": mnu["citation"],
        })
        for item_id in ("planck2018_vi_bao_neff", "act_dr6_p_act_lb_neff"):
            item = self.items[item_id]
            self.add("not_evaluable", {
                "row_id": f"neff_{item_id}_not_evaluable",
                "class": "not_evaluable",
                "oph_value": None,
                "oph_value_exact_rational": None,
                "measurement": {"kind": "gaussian", "value": fl(item["value"]), "sigma": fl(item["uncertainty"])},
                "public_input_ids": [item_id],
                "sigma_distance": None,
                "sigma_distance_kind": "none_not_evaluable",
                "verdict": "not_evaluable",
                "discriminates_from_baseline": False,
                "prospective_data": "none: no corpus value",
                "dataset": item["dataset"],
                "citation": item["citation"],
            })


def compare_measurement(actual: dict[str, Any], expected: dict[str, Any], label: str) -> None:
    require(isinstance(actual, dict), f"{label}: measurement object")
    for key, value in expected.items():
        require(key in actual, f"{label}: measurement field {key} absent")
        if value is None:
            require(actual[key] is None, f"{label}: {key} expected null")
        elif isinstance(value, float):
            same_float(actual[key], value, f"{label}: {key}")
        elif isinstance(value, list):
            require(len(actual[key]) == len(value), f"{label}: {key} length")
            for a, b in zip(actual[key], value):
                same_float(a, b, f"{label}: {key} entry")
        else:
            require(actual[key] == value, f"{label}: {key} value")


def compare_extras(actual: Any, expected: Any, label: str) -> None:
    if isinstance(expected, dict):
        require(isinstance(actual, dict) and actual.keys() == expected.keys(), f"{label}: extras keys")
        for key in expected:
            compare_extras(actual[key], expected[key], f"{label}/{key}")
    elif isinstance(expected, float):
        same_float(actual, expected, label)
    else:
        require(actual == expected, f"{label}: value")


def verify_row(actual: dict[str, Any], expected: dict[str, Any], section_id: str) -> None:
    label = expected["row_id"]
    keys = set(actual)
    require(ROW_KEYS <= keys and keys <= ROW_KEYS | {"extras"}, f"{label}: row keys {sorted(keys ^ ROW_KEYS)}")
    require(actual["section"] == section_id == expected["section"], f"{label}: section")
    require(isinstance(actual["class"], str) and actual["class"] in CLASSES, f"{label}: class")
    require(actual["class"] == expected["class"], f"{label}: class value")
    expected_oph = expected["oph_value"]
    if expected_oph is None:
        require(actual["oph_value"] is None, f"{label}: oph_value null")
    elif isinstance(expected_oph, dict):
        require(isinstance(actual["oph_value"], dict) and actual["oph_value"].keys() == expected_oph.keys(), f"{label}: oph_value keys")
        for key in expected_oph:
            same_float(actual["oph_value"][key], expected_oph[key], f"{label}: oph_value {key}")
    else:
        same_float(actual["oph_value"], expected_oph, f"{label}: oph_value")
    require(actual["oph_value_exact_rational"] == expected["oph_value_exact_rational"], f"{label}: exact rational")
    if expected["oph_value_exact_rational"] is not None:
        exact = Fraction(expected["oph_value_exact_rational"])
        require(float(exact) == float(actual["oph_value"]), f"{label}: exact rational and float disagree")
    compare_measurement(actual["measurement"], expected["measurement"], label)
    require(actual["public_input_ids"] == expected["public_input_ids"], f"{label}: public input ids")
    if expected["sigma_distance"] is None:
        require(actual["sigma_distance"] is None, f"{label}: sigma expected null")
    else:
        same_float(actual["sigma_distance"], expected["sigma_distance"], f"{label}: sigma")
        if actual["sigma_distance_kind"] == "signed_one_dimensional":
            difference = float(actual["oph_value"]) - float(actual["measurement"]["value"])
            require(
                (difference > 0) == (expected["sigma_distance"] > 0) and (difference == 0) == (expected["sigma_distance"] == 0),
                f"{label}: sigma sign convention",
            )
        require(actual["verdict"] == rule(expected["sigma_distance"]), f"{label}: verdict rule")
    require(actual["sigma_distance_kind"] == expected["sigma_distance_kind"], f"{label}: sigma kind")
    require(actual["verdict"] in VERDICTS and actual["verdict"] == expected["verdict"], f"{label}: verdict")
    if actual["measurement"].get("kind") == "upper_bound" and actual["class"] != "not_evaluable":
        below = float(actual["oph_value"]) < float(actual["measurement"]["upper_bound"])
        require(actual["verdict"] == ("consistent" if below else "exceeds_stated_confidence_bound"), f"{label}: bound rule")
    require(actual["discriminates_from_baseline"] is expected["discriminates_from_baseline"], f"{label}: discriminates")
    require(actual["prospective_data"] == expected["prospective_data"], f"{label}: prospective data")
    for key in ("dataset", "citation", "oph_premises", "epistemic_note"):
        if key in expected:
            require(actual[key] == expected[key], f"{label}: {key}")
    for key in ("quantity", "oph_value_display", "measurement_display", "dataset", "citation", "epistemic_note", "prospective_data"):
        require(isinstance(actual[key], str) and actual[key].strip(), f"{label}: {key} prose")
    require(isinstance(actual["oph_premises"], list) and actual["oph_premises"] and all(isinstance(p, str) and p for p in actual["oph_premises"]), f"{label}: premises")
    require(isinstance(actual["corpus_source_ids"], list), f"{label}: corpus source ids")
    if "extras" in expected:
        require("extras" in actual, f"{label}: extras absent")
        compare_extras(actual["extras"], expected["extras"], f"{label}: extras")
    if actual["row_id"].startswith("ns_clock_branch"):
        require(
            "diagnostic alternative branch; no selection between branches is made by data" in actual["epistemic_note"],
            f"{label}: clock-branch note",
        )
    if actual["class"] == "shared_baseline":
        require(actual["discriminates_from_baseline"] is False, f"{label}: shared baseline discriminates")


def verify(receipt_path: Path = RECEIPT_PATH, markdown_path: Path = MARKDOWN_PATH) -> dict[str, Any]:
    receipt_bytes = receipt_path.read_bytes()
    require(b"\r" not in receipt_bytes and receipt_bytes.endswith(b"\n"), "receipt line endings")
    require(all(byte < 128 for byte in receipt_bytes), "receipt is not ASCII")
    receipt = load_exact(receipt_path)
    require(set(receipt) == TOP_LEVEL_KEYS, f"top-level keys {sorted(set(receipt) ^ TOP_LEVEL_KEYS)}")
    require(receipt["schema"] == SCHEMA, "schema")
    require(receipt["classification"] == CLASSIFICATION, "classification")
    require(receipt["promotes_nothing"] is True and receipt["seen_data"] is True, "promotion flags")
    require(receipt["frozen_prediction_rows"] == 0, "frozen prediction rows")
    require(set(receipt["allowed_classes"]) == CLASSES and set(receipt["allowed_verdicts"]) == VERDICTS, "allowed sets")
    require(set(receipt["class_descriptions"]) == CLASSES, "class descriptions")
    require(receipt["public_inputs_schema"] == PUBLIC_SCHEMA and receipt["corpus_inputs_schema"] == CORPUS_SCHEMA, "input schemas")

    public_bytes, corpus_bytes = PUBLIC_PATH.read_bytes(), CORPUS_PATH.read_bytes()
    public, corpus = load_exact(PUBLIC_PATH), load_exact(CORPUS_PATH)
    require(public["schema"] == PUBLIC_SCHEMA and corpus["schema"] == CORPUS_SCHEMA, "input file schemas")

    expected_pins = {
        f"{PACKAGE}/public_inputs.json": sha256_bytes(public_bytes),
        f"{PACKAGE}/corpus_inputs.json": sha256_bytes(corpus_bytes),
    }
    for relative in PINNED_CODE:
        expected_pins[relative] = sha256_bytes((ROOT / relative).read_bytes())
    parents: dict[str, Any] = {}
    expected_parent_pins = {}
    for key, entry in corpus["parent_receipts"].items():
        path = ROOT / entry["path"]
        payload = path.read_bytes()
        require(sha256_bytes(payload) == entry["sha256"], f"{key}: parent receipt sha256")
        data = load_exact(path)
        require(data["schema"] == entry["schema"], f"{key}: parent schema")
        parents[key] = data
        expected_pins[entry["path"]] = entry["sha256"]
        expected_parent_pins[key] = {"path": entry["path"], "sha256": entry["sha256"], "schema": entry["schema"]}
    require(receipt["pins"] == expected_pins, "pins")
    require(receipt["parent_receipts"] == expected_parent_pins, "parent receipt pins")

    expected_anchors = []
    for key, entry in corpus["typed"].items():
        if entry["source_file"] is None:
            require(entry["anchor_text"] is None, f"{key}: anchor without source")
            continue
        text = (ROOT / entry["source_file"]).read_text(encoding="utf-8")
        require(norm(entry["anchor_text"]) in norm(text), f"{key}: anchor text absent")
        expected_anchors.append({"corpus_input_id": key, "source_file": entry["source_file"], "anchor_verified": True})
    require(receipt["corpus_anchor_checks"] == expected_anchors, "corpus anchor checks")

    expected = Expected(public, corpus, parents)
    require(set(receipt["derived_constants"]) == set(expected.derived), "derived constant keys")
    for key, value in expected.derived.items():
        same_float(receipt["derived_constants"][key], value, f"derived constant {key}")

    sections = receipt["sections"]
    require([s["id"] for s in sections] == ["primordial", "background_dark_energy", "capacity_lambda", "dark_sector", "not_evaluable"], "section order")
    actual_order = [(s["id"], r["row_id"]) for s in sections for r in s["rows"]]
    require(actual_order == expected.order, "row order or row set")
    for section in sections:
        require(set(section) == {"id", "title", "rows"}, "section keys")
        for row in section["rows"]:
            verify_row(row, expected.rows[row["row_id"]], section["id"])
    require(receipt["row_count"] == len(actual_order), "row count")
    canonical = json.dumps(json.loads(json.dumps(sections, default=float)), ensure_ascii=True, separators=(",", ":"), sort_keys=True).encode("utf-8")
    require(receipt["rows_sha256"] == sha256_bytes(canonical), "rows digest")

    problems: list[str] = []
    scan_prose(json.loads(receipt_bytes.decode("utf-8")), "receipt", problems)
    markdown = markdown_path.read_text(encoding="utf-8")
    require(markdown.startswith("# Cosmology postdiction ledger\n"), "markdown title")
    require(sha256_bytes(receipt_bytes) in markdown and receipt["rows_sha256"] in markdown, "markdown digest lines")
    require("promotes nothing" in markdown and "## Reproduce" in markdown, "markdown header")
    for section in sections:
        require(f"## {section['title']}" in markdown, f"markdown section {section['id']}")
        for row in section["rows"]:
            require(row["verdict"] in markdown, f"markdown verdict {row['row_id']}")
    prose_only = "\n".join(
        line for line in markdown.splitlines() if not re.fullmatch(r"\|(\s*-{3,}\s*\|)+\s*", line)
    )
    scan_prose(prose_only, "markdown", problems)
    require(not problems, "; ".join(problems[:10]))

    verdicts: dict[str, int] = {}
    for _section, row_id in actual_order:
        verdict = expected.rows[row_id]["verdict"]
        verdicts[verdict] = verdicts.get(verdict, 0) + 1
    return {"schema": SCHEMA, "rows": len(actual_order), "verdict_counts": verdicts, "promotes_nothing": True}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--receipt", type=Path, default=RECEIPT_PATH)
    parser.add_argument("--markdown", type=Path, default=MARKDOWN_PATH)
    args = parser.parse_args()
    try:
        summary = verify(args.receipt, args.markdown)
    except (VerificationError, KeyError, TypeError, ValueError, OSError) as exc:
        print(f"FAIL: {exc}")
        return 1
    print("PASS: " + json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
