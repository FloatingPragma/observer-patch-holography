#!/usr/bin/env python3
"""Audit FZ-13 CPL subsets and a base-LCDM display on official DESI chains.

This is a retrospective comparison.  DESI DR2 was public before this audit
and before any future FZ-13 freeze, so the output can diagnose conditional
branches but cannot score a prediction or confirm OPH.

The input files are the collaboration-produced Cobaya chains documented at
https://data.desi.lbl.gov/public/papers/y3/bao-cosmo-params/README.html.
Their hashes below are copied from the collaboration's v1.0 SHA-256 manifest.
No OPH code generated or refit these posterior samples.

The base-LambdaCDM companion calculation evaluates ``Lambda*l_P^2`` on each
weighted ``(H0, omegal)`` sample.  It therefore uses the chain covariance
directly and does not assume a correlation coefficient.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from statistics import NormalDist
from typing import Iterable


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
SCRIPT_PATH = Path(__file__).resolve()
SOURCE_COBAYA_ROOT = (
    "https://data.desi.lbl.gov/public/papers/y3/bao-cosmo-params/cobaya"
)
SOURCE_MANIFEST = (
    "https://data.desi.lbl.gov/public/papers/y3/bao-cosmo-params/"
    "dr2_vac_dr2_bao-cosmo-params_v1.0.sha256sum"
)
SOURCE_MANIFEST_SHA256 = (
    "df78872aa8b2d3473a9e8de78f498180efd7cbcbeb18211ce4787fac52067ee5"
)

# Central constants used only to form the dimensionless retrospective display
# Lambda*l_P^2 = 3 Omega_Lambda (H0/c)^2 l_P^2 sample by sample.  The exact SI
# speed of light and the stated central conversion values are recorded in the
# receipt; their tiny uncertainties are not folded into the cosmology posterior.
SPEED_OF_LIGHT_M_S = 299_792_458.0
MPC_IN_M = 3.0856775814913673e22
PLANCK_LENGTH_M = 1.616255e-35
PLANCK_LENGTH_STANDARD_UNCERTAINTY_M = 0.000018e-35
CODATA_2022_PLANCK_LENGTH_SOURCE = "https://physics.nist.gov/cgi-bin/cuu/Value?plkl"

_CMB = (
    "planck2018-lowl-TT-clik_planck2018-lowl-EE-clik_"
    "planck-NPIPE-highl-CamSpec-TTTEEE_planck-act-dr6-lensing"
)


DATASETS = {
    "DESI_DR2_BAO+CMB": {
        "slug": "cmb",
        "model": "base_w_wa",
        "directory": f"desi-bao-all_{_CMB}",
        "sha256": [
            "c228de7bbaec19ddb22eec25c3dd7c40ef218976c020a7b55dd1c78dc3a638c5",
            "01bb30f43b3207d8575cc16354159fceff2ad4deffacc87964b0aef7f8e8ee44",
            "db12623c8c69c03ad219b28bbc517416fc941cbb1e717054344f30b4e43a4adc",
            "339312d9b7d3027147c38433fb4335cb9a12585776c4dfffc592c7c43af9a5ff",
        ],
    },
    "DESI_DR2_BAO+CMB+PantheonPlus": {
        "slug": "pantheon",
        "model": "base_w_wa",
        "directory": f"desi-bao-all_pantheonplus_{_CMB}",
        "sha256": [
            "db81d299d59051ae5e1e8f67952320e08a1644a4a1f8d19267705aee32798877",
            "442390266f6a3bfe88e9c7cd8e29d1e7cff47fc8582f3a3af80a6b40b9f83939",
            "1c92edf523df34784633f6852f7564f3d953203ae939d5d2e8cd4c3fd6f580ae",
            "b17dc02689c3a30dd34a96d91ac35180006f17d10b82331396ef55ce999594af",
        ],
    },
    "DESI_DR2_BAO+CMB+Union3": {
        "slug": "union3",
        "model": "base_w_wa",
        "directory": f"desi-bao-all_union3_{_CMB}",
        "sha256": [
            "f95d091203f615e1654ea9f4fea332db08fa66c59a5fe188b4915e67e1d73b11",
            "23aac58326257d207b8ee3699d21ed908e46238e42041d8c2ccf12fd85c8b839",
            "9e0998c5685e7552b94fd45e985dd137604f34692bebbc6fb922a1588669ab6a",
            "7710ababc2ab0c7fb5f3f67c8b790409703bce18707a2f99a1cd3fea531b4c20",
        ],
    },
    "DESI_DR2_BAO+CMB+DESY5": {
        "slug": "desy5",
        "model": "base_w_wa",
        "directory": f"desi-bao-all_desy5sn_{_CMB}",
        "sha256": [
            "8c783ebf283a205b7f569ce36a6694a32646ced5e28fd3b4683742733f6165e0",
            "cd4f2ff3a66aa92aceecd47c8452520c2de09d887f231fc0df033cd68597a886",
            "f7f8dbf28ff23d371e0b987b938d321adb920f1a09e97f746b3c0bb2d6247a01",
            "4a3607867d34890832f431c4d61947e5572a72ff1a407141f4b8cb8d3d7ec769",
        ],
    },
}

BASE_LCDM_NAME = "DESI_DR2_BAO+CMB_base_LCDM"
BASE_LCDM_DATASET = {
    "slug": "lcdm",
    "model": "base",
    "directory": f"desi-bao-all_{_CMB}",
    "sha256": [
        "00f3766f7a7b6370d21323886cd72869087b2b1346a04d729c8f3bc9e65ef698",
        "33b154eebdf4e9dca3b8f02ed2680120879d35c10b32fef42261a490104e1dc1",
        "d4717e7e5a13de851c86f24c87213faccef2b5f8747900274ab509d9dfa40aa2",
        "c827cd767a4864ca28aa15c902bda32004e803050d4be330e25aefddd78b5c36",
    ],
}

DOWNLOAD_SPECS = {**DATASETS, BASE_LCDM_NAME: BASE_LCDM_DATASET}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def producer_metadata() -> dict[str, str]:
    """Bind a receipt to the exact postprocessor that emitted it."""

    return {
        "script_path": str(SCRIPT_PATH.relative_to(REPO_ROOT)),
        "script_sha256": sha256(SCRIPT_PATH),
    }


def lambda_lp2_from_base_lcdm_sample(h0_km_s_mpc: float, omega_lambda: float) -> float:
    """Return ``Lambda*l_P^2`` for one flat base-LCDM posterior sample.

    The chain's derived ``omegal`` column is used directly; no replacement by
    ``1-omegam`` and no fitted H0--density correlation is introduced.
    """

    if not (math.isfinite(h0_km_s_mpc) and math.isfinite(omega_lambda)):
        raise ValueError("H0 and OmegaLambda must be finite")
    if h0_km_s_mpc <= 0 or not 0.0 < omega_lambda < 1.0:
        raise ValueError("H0 must be positive and OmegaLambda must lie in (0,1)")
    h0_s = h0_km_s_mpc * 1_000.0 / MPC_IN_M
    return 3.0 * omega_lambda * (h0_s / SPEED_OF_LIGHT_M_S) ** 2 * PLANCK_LENGTH_M**2


def weighted_quantiles(
    samples: list[tuple[float, float]], probabilities: tuple[float, ...]
) -> dict[str, float]:
    """Return step-CDF quantiles of positive-weight samples."""

    if not samples:
        raise ValueError("cannot compute quantiles of an empty sample")
    if any(not 0.0 <= probability <= 1.0 for probability in probabilities):
        raise ValueError("quantile probabilities must lie in [0,1]")
    if tuple(sorted(probabilities)) != probabilities:
        raise ValueError("quantile probabilities must be sorted")
    ordered = sorted(samples)
    if any(
        not (math.isfinite(value) and math.isfinite(weight)) or weight <= 0
        for value, weight in ordered
    ):
        raise ValueError("quantile samples require finite values and positive weights")
    total = sum(weight for _value, weight in ordered)
    result: dict[str, float] = {}
    cumulative = 0.0
    probability_index = 0
    for value, weight in ordered:
        cumulative += weight
        while (
            probability_index < len(probabilities)
            and cumulative >= probabilities[probability_index] * total
        ):
            probability = probabilities[probability_index]
            result[f"{probability:.3f}"] = value
            probability_index += 1
    if probability_index != len(probabilities):
        raise RuntimeError("weighted quantile traversal ended early")
    return result


@dataclass
class BaseLCDMAccumulator:
    raw_rows: int = 0
    sum_weight: float = 0.0
    sum_weight_sq: float = 0.0
    sum_h0: float = 0.0
    sum_h0_sq: float = 0.0
    sum_omega_lambda: float = 0.0
    sum_omega_lambda_sq: float = 0.0
    sum_h0_omega_lambda: float = 0.0
    sum_lambda_lp2: float = 0.0
    sum_lambda_lp2_sq: float = 0.0
    lambda_lp2_samples: list[tuple[float, float]] | None = None

    def __post_init__(self) -> None:
        if self.lambda_lp2_samples is None:
            self.lambda_lp2_samples = []

    def add(self, weight: float, h0: float, omega_lambda: float) -> None:
        if not math.isfinite(weight) or weight <= 0:
            raise ValueError("chain weights must be finite and positive")
        lambda_lp2 = lambda_lp2_from_base_lcdm_sample(h0, omega_lambda)
        self.raw_rows += 1
        self.sum_weight += weight
        self.sum_weight_sq += weight * weight
        self.sum_h0 += weight * h0
        self.sum_h0_sq += weight * h0 * h0
        self.sum_omega_lambda += weight * omega_lambda
        self.sum_omega_lambda_sq += weight * omega_lambda * omega_lambda
        self.sum_h0_omega_lambda += weight * h0 * omega_lambda
        self.sum_lambda_lp2 += weight * lambda_lp2
        self.sum_lambda_lp2_sq += weight * lambda_lp2 * lambda_lp2
        assert self.lambda_lp2_samples is not None
        self.lambda_lp2_samples.append((lambda_lp2, weight))

    def merge(self, other: "BaseLCDMAccumulator") -> None:
        self.raw_rows += other.raw_rows
        self.sum_weight += other.sum_weight
        self.sum_weight_sq += other.sum_weight_sq
        self.sum_h0 += other.sum_h0
        self.sum_h0_sq += other.sum_h0_sq
        self.sum_omega_lambda += other.sum_omega_lambda
        self.sum_omega_lambda_sq += other.sum_omega_lambda_sq
        self.sum_h0_omega_lambda += other.sum_h0_omega_lambda
        self.sum_lambda_lp2 += other.sum_lambda_lp2
        self.sum_lambda_lp2_sq += other.sum_lambda_lp2_sq
        assert self.lambda_lp2_samples is not None
        assert other.lambda_lp2_samples is not None
        self.lambda_lp2_samples.extend(other.lambda_lp2_samples)

    def summary(self) -> dict[str, object]:
        if self.sum_weight <= 0:
            raise ValueError("empty chain")
        mean_h0 = self.sum_h0 / self.sum_weight
        mean_omega_lambda = self.sum_omega_lambda / self.sum_weight
        mean_lambda_lp2 = self.sum_lambda_lp2 / self.sum_weight
        var_h0 = max(0.0, self.sum_h0_sq / self.sum_weight - mean_h0**2)
        var_omega_lambda = max(
            0.0,
            self.sum_omega_lambda_sq / self.sum_weight - mean_omega_lambda**2,
        )
        var_lambda_lp2 = max(
            0.0,
            self.sum_lambda_lp2_sq / self.sum_weight - mean_lambda_lp2**2,
        )
        covariance = (
            self.sum_h0_omega_lambda / self.sum_weight - mean_h0 * mean_omega_lambda
        )
        correlation = (
            covariance / math.sqrt(var_h0 * var_omega_lambda)
            if var_h0 > 0 and var_omega_lambda > 0
            else None
        )
        assert self.lambda_lp2_samples is not None
        return {
            "raw_rows": self.raw_rows,
            "expanded_posterior_weight": self.sum_weight,
            "weight_concentration_ess_not_autocorrelation_corrected": (
                self.sum_weight**2 / self.sum_weight_sq
            ),
            "H0_km_s_Mpc": {
                "weighted_mean": mean_h0,
                "weighted_std": math.sqrt(var_h0),
            },
            "OmegaLambda": {
                "weighted_mean": mean_omega_lambda,
                "weighted_std": math.sqrt(var_omega_lambda),
            },
            "H0_OmegaLambda_weighted_correlation": correlation,
            "Lambda_lP2": {
                "weighted_mean": mean_lambda_lp2,
                "weighted_std": math.sqrt(var_lambda_lp2),
                "fractional_std_about_weighted_mean": (
                    math.sqrt(var_lambda_lp2) / mean_lambda_lp2
                ),
                "weighted_step_cdf_quantiles": weighted_quantiles(
                    self.lambda_lp2_samples, (0.025, 0.16, 0.5, 0.84, 0.975)
                ),
            },
        }


@dataclass
class Accumulator:
    raw_rows: int = 0
    sum_weight: float = 0.0
    sum_weight_sq: float = 0.0
    sum_w0: float = 0.0
    sum_wa: float = 0.0
    sum_w0_sq: float = 0.0
    sum_wa_sq: float = 0.0
    sum_w0_wa: float = 0.0
    monotone_weight: float = 0.0
    monotone_rows: int = 0
    w0_gt_neg_one_weight: float = 0.0
    wa_nonneg_weight: float = 0.0

    def add(self, weight: float, w0: float, wa: float) -> None:
        if not (math.isfinite(weight) and math.isfinite(w0) and math.isfinite(wa)):
            raise ValueError("chain contains a non-finite weight, w, or wa")
        if weight <= 0:
            raise ValueError("chain weights must be positive")
        self.raw_rows += 1
        self.sum_weight += weight
        self.sum_weight_sq += weight * weight
        self.sum_w0 += weight * w0
        self.sum_wa += weight * wa
        self.sum_w0_sq += weight * w0 * w0
        self.sum_wa_sq += weight * wa * wa
        self.sum_w0_wa += weight * w0 * wa

        # CPL is affine in a: w(a)=w0+wa(1-a).  On a in [1/3,1],
        # w(a)>=-1 iff the inequality holds at both endpoints.
        monotone = w0 >= -1.0 and w0 + (2.0 / 3.0) * wa >= -1.0
        if monotone:
            self.monotone_weight += weight
            self.monotone_rows += 1
        if w0 > -1.0:
            self.w0_gt_neg_one_weight += weight
        if wa >= 0.0:
            self.wa_nonneg_weight += weight

    def merge(self, other: "Accumulator") -> None:
        for field in self.__dataclass_fields__:
            setattr(self, field, getattr(self, field) + getattr(other, field))

    def summary(self) -> dict[str, float | int | None]:
        if self.sum_weight <= 0:
            raise ValueError("empty chain")
        mean_w0 = self.sum_w0 / self.sum_weight
        mean_wa = self.sum_wa / self.sum_weight
        var_w0 = max(0.0, self.sum_w0_sq / self.sum_weight - mean_w0**2)
        var_wa = max(0.0, self.sum_wa_sq / self.sum_weight - mean_wa**2)
        cov = self.sum_w0_wa / self.sum_weight - mean_w0 * mean_wa
        monotone_mass = self.monotone_weight / self.sum_weight
        correlation = (
            cov / math.sqrt(var_w0 * var_wa) if var_w0 > 0 and var_wa > 0 else None
        )
        return {
            "raw_rows": self.raw_rows,
            "expanded_posterior_weight": self.sum_weight,
            "weight_concentration_ess_not_autocorrelation_corrected": (
                self.sum_weight**2 / self.sum_weight_sq
            ),
            "w0_mean": mean_w0,
            "w0_std": math.sqrt(var_w0),
            "wa_mean": mean_wa,
            "wa_std": math.sqrt(var_wa),
            "w0_wa_covariance": cov,
            "w0_wa_correlation": correlation,
            "posterior_mass_w_ge_minus_one_for_0_le_z_le_2": monotone_mass,
            "posterior_mass_capacity_loss_somewhere_for_0_le_z_le_2": (
                1.0 - monotone_mass
            ),
            "posterior_mass_w0_gt_minus_one": (
                self.w0_gt_neg_one_weight / self.sum_weight
            ),
            "posterior_mass_wa_nonnegative": self.wa_nonneg_weight / self.sum_weight,
            "raw_rows_in_monotone_subset": self.monotone_rows,
        }


def read_chain(path: Path) -> Accumulator:
    header: list[str] | None = None
    out = Accumulator()
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            if line.startswith("#"):
                if header is None:
                    header = line.lstrip("#").split()
                continue
            if header is None:
                raise ValueError(f"{path}:{line_number}: missing header")
            values = line.split()
            if len(values) != len(header):
                raise ValueError(
                    f"{path}:{line_number}: {len(values)} values for {len(header)} columns"
                )
            row = dict(zip(header, values, strict=True))
            try:
                out.add(float(row["weight"]), float(row["w"]), float(row["wa"]))
            except KeyError as exc:
                raise ValueError(
                    f"{path}: required column absent: {exc.args[0]}"
                ) from exc
    return out


def read_base_lcdm_chain(path: Path) -> BaseLCDMAccumulator:
    header: list[str] | None = None
    out = BaseLCDMAccumulator()
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            if line.startswith("#"):
                if header is None:
                    header = line.lstrip("#").split()
                continue
            if header is None:
                raise ValueError(f"{path}:{line_number}: missing header")
            values = line.split()
            if len(values) != len(header):
                raise ValueError(
                    f"{path}:{line_number}: {len(values)} values for {len(header)} columns"
                )
            row = dict(zip(header, values, strict=True))
            try:
                out.add(
                    float(row["weight"]),
                    float(row["H0"]),
                    float(row["omegal"]),
                )
            except KeyError as exc:
                raise ValueError(
                    f"{path}: required column absent: {exc.args[0]}"
                ) from exc
    return out


def source_directory(spec: dict[str, object]) -> str:
    return f"{SOURCE_COBAYA_ROOT}/{spec['model']}/{spec['directory']}/"


def gaussian_fixed_point_diagnostic(
    summary: dict[str, float | int],
) -> dict[str, float | str]:
    """Return a labelled Gaussian moment diagnostic, never a likelihood verdict."""

    v0 = float(summary["w0_std"]) ** 2
    va = float(summary["wa_std"]) ** 2
    cov = float(summary["w0_wa_covariance"])
    det = v0 * va - cov * cov
    if det <= 0:
        raise ValueError("non-positive posterior covariance determinant")
    d0 = float(summary["w0_mean"]) + 1.0
    da = float(summary["wa_mean"])
    q = (va * d0 * d0 - 2.0 * cov * d0 * da + v0 * da * da) / det
    survival_chi2_2dof = math.exp(-0.5 * q)
    sigma_two_sided = NormalDist().inv_cdf(1.0 - survival_chi2_2dof / 2.0)
    return {
        "classification": "Gaussian moment summary; not official delta-chi2 or evidence",
        "mahalanobis_squared": q,
        "chi2_2dof_survival": survival_chi2_2dof,
        "two_sided_normal_sigma_equivalent": sigma_two_sided,
    }


def audit_dataset(data_dir: Path, spec: dict[str, object]) -> dict[str, object]:
    combined = Accumulator()
    chains: list[dict[str, object]] = []
    expected_hashes = list(spec["sha256"])
    for index, expected in enumerate(expected_hashes, start=1):
        path = data_dir / f"{spec['slug']}_chain.{index}.txt"
        actual = sha256(path)
        if actual != expected:
            raise ValueError(f"SHA-256 mismatch for {path}: {actual} != {expected}")
        accumulator = read_chain(path)
        summary = accumulator.summary()
        chains.append(
            {
                "chain": index,
                "file": path.name,
                "bytes": path.stat().st_size,
                "sha256": actual,
                **summary,
            }
        )
        combined.merge(accumulator)

    summary = combined.summary()
    masses = [
        float(chain["posterior_mass_w_ge_minus_one_for_0_le_z_le_2"])
        for chain in chains
    ]
    raw_tail_rows = int(summary["raw_rows_in_monotone_subset"])
    return {
        "source_directory": source_directory(spec),
        "chains": chains,
        "combined": summary,
        "chain_range_for_monotone_subset_mass": [min(masses), max(masses)],
        "rare_tail_resolution": {
            "classification": (
                "resolved_in_all_four_chains"
                if min(int(chain["raw_rows_in_monotone_subset"]) for chain in chains)
                > 0
                else "at_least_one_chain_has_no_raw_tail_row"
            ),
            "combined_raw_tail_rows": raw_tail_rows,
            "warning": (
                "The weighted chain fraction is an empirical posterior mass under the DESI "
                "model, priors, and likelihoods. It is not a frequentist exclusion, a direct "
                "capacity measurement, or an OPH prediction score."
            ),
        },
        "fixed_capacity_point_gaussian_diagnostic": gaussian_fixed_point_diagnostic(
            summary
        ),
    }


def audit_base_lcdm_dataset(
    data_dir: Path, spec: dict[str, object] = BASE_LCDM_DATASET
) -> dict[str, object]:
    """Compute the nonlinear Lambda*l_P^2 posterior from official samples."""

    combined = BaseLCDMAccumulator()
    chains: list[dict[str, object]] = []
    for index, expected in enumerate(spec["sha256"], start=1):
        path = data_dir / f"{spec['slug']}_chain.{index}.txt"
        actual = sha256(path)
        if actual != expected:
            raise ValueError(f"SHA-256 mismatch for {path}: {actual} != {expected}")
        accumulator = read_base_lcdm_chain(path)
        chains.append(
            {
                "chain": index,
                "file": path.name,
                "bytes": path.stat().st_size,
                "sha256": actual,
                **accumulator.summary(),
            }
        )
        combined.merge(accumulator)
    return {
        "source_directory": source_directory(spec),
        "model_scope": (
            "flat base-LambdaCDM posterior under the official DESI DR2 BAO+CMB "
            "model, priors, and likelihoods"
        ),
        "sample_level_formula": (
            "Lambda*l_P^2 = 3*omegal*(H0*1000/Mpc_in_m/c)^2*l_P^2"
        ),
        "constants": {
            "speed_of_light_m_s_exact": SPEED_OF_LIGHT_M_S,
            "Mpc_in_m": MPC_IN_M,
            "Planck_length_m_CODATA_2022_central": PLANCK_LENGTH_M,
            "Planck_length_standard_uncertainty_m": (
                PLANCK_LENGTH_STANDARD_UNCERTAINTY_M
            ),
            "Planck_length_source": CODATA_2022_PLANCK_LENGTH_SOURCE,
            "constants_uncertainty_propagated": False,
        },
        "chains": chains,
        "combined": combined.summary(),
        "classification": (
            "sample-level retrospective display posterior; not an OPH prediction, "
            "model evidence, or model-independent measurement"
        ),
    }


def build_receipt(
    data_dir: Path, selected: Iterable[str] | None = None
) -> dict[str, object]:
    names = list(DATASETS) if selected is None else list(selected)
    unknown = sorted(set(names) - set(DATASETS))
    if unknown:
        raise ValueError(f"unknown datasets: {unknown}")
    return {
        "schema": "oph.official_desi_dr2_fz13_retrospective.v2",
        "producer": producer_metadata(),
        "source": {
            "publisher": "DESI Collaboration / DESI Data",
            "documentation": (
                "https://data.desi.lbl.gov/public/papers/y3/"
                "bao-cosmo-params/README.html"
            ),
            "official_sha256_manifest": SOURCE_MANIFEST,
            "official_sha256_manifest_sha256": SOURCE_MANIFEST_SHA256,
            "release_class": "official collaboration posterior chains",
        },
        "epistemic_status": {
            "retrospective_seen_data": True,
            "frozen_prediction_score": False,
            "oph_confirmation": False,
            "direct_capacity_measurement": False,
            "conditional_model_test_only": True,
            "notes": (
                "The CPL inequalities test only the conditional FZ-13 map after its density, "
                "closed-sector, and capacity-history premises. Agreement with (-1,0) is shared "
                "with LambdaCDM and earns no OPH-specific support."
            ),
        },
        "subset_definition": {
            "cpl": "w(a)=w0+wa(1-a)",
            "redshift_range": "0 <= z <= 2, equivalently 1/3 <= a <= 1",
            "monotone_capacity_condition": (
                "w(a)>=-1 on the full interval; because CPL is affine, this is exactly "
                "w0>=-1 and w0+(2/3)wa>=-1"
            ),
        },
        "datasets": {name: audit_dataset(data_dir, DATASETS[name]) for name in names},
        "base_lcdm_capacity_display": audit_base_lcdm_dataset(data_dir),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--dataset", action="append", choices=sorted(DATASETS))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    receipt = build_receipt(args.data_dir, args.dataset)
    rendered = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
