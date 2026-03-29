#!/usr/bin/env python3
"""Build the perturbative QCD Lambda bridge used by the hadron lane.

Chain role: connect the live D10/QCD calibration data and local quark masses to
the `Lambda_MSbar^(3)` descendant that seeds the hadron ensemble family.

Mathematics: perturbative MS-bar beta running with threshold-style descent from
the live strong-coupling input.

OPH-derived inputs: the D10 calibration core and the local quark singular-value
masses from the active flavor chain.

Output: the `Lambda_MSbar_3_gev` bridge consumed by the unquenched hadron
producer.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_D10 = ROOT / "particles" / "runs" / "calibration" / "d10_ew_observable_family.json"
DEFAULT_FORWARD = ROOT / "particles" / "runs" / "flavor" / "forward_yukawas.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "qcd" / "lambda_msbar_descendant.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def beta_coeffs_msbar(n_f: int) -> tuple[float, float, float, float]:
    """Return (beta0, beta1, beta2, beta3) for SU(3) QCD in MS-bar, a = alpha_s/(4*pi)."""
    nf = float(n_f)
    zeta3 = 1.2020569031595942
    beta0 = 11.0 - (2.0 / 3.0) * nf
    beta1 = 102.0 - (38.0 / 3.0) * nf
    beta2 = 2857.0 / 2.0 - (5033.0 / 18.0) * nf + (325.0 / 54.0) * (nf ** 2)
    beta3 = (
        149753.0 / 6.0
        + 3564.0 * zeta3
        + (-1078361.0 / 162.0 - (6508.0 / 27.0) * zeta3) * nf
        + (50065.0 / 162.0 + (6472.0 / 81.0) * zeta3) * (nf ** 2)
        + (1093.0 / 729.0) * (nf ** 3)
    )
    return (beta0, beta1, beta2, beta3)


def beta_a(a_value: float, n_f: int, loops: int = 4) -> float:
    if a_value <= 0.0:
        raise ValueError("a must be > 0")
    beta0, beta1, beta2, beta3 = beta_coeffs_msbar(n_f)
    out = -beta0 * a_value * a_value
    if loops >= 2:
        out += -beta1 * a_value ** 3
    if loops >= 3:
        out += -beta2 * a_value ** 4
    if loops >= 4:
        out += -beta3 * a_value ** 5
    return out


def _simpson(f: Callable[[float], float], start: float, stop: float, panels: int) -> float:
    if panels % 2:
        raise ValueError("panels must be even")
    h = (stop - start) / panels
    acc = f(start) + f(stop)
    for idx in range(1, panels):
        x_value = start + idx * h
        acc += f(x_value) * (4 if idx % 2 else 2)
    return acc * h / 3.0


def lambda_msbar_from_alpha(mu_gev: float, alpha: float, n_f: int, loops: int = 4, panels: int = 20000) -> float:
    if mu_gev <= 0.0:
        raise ValueError("mu_gev must be > 0")
    if alpha <= 0.0:
        raise ValueError("alpha must be > 0")
    if not 1 <= loops <= 4:
        raise ValueError("loops must be 1..4")
    if panels % 2:
        panels += 1

    a_value = alpha / (4.0 * math.pi)
    beta0, beta1, _, _ = beta_coeffs_msbar(n_f)
    running = 1.0 / (beta0 * a_value)
    if loops >= 2:
        running += (beta1 / (beta0 * beta0)) * math.log(beta0 * a_value)
    if loops >= 3:
        eps = max(1.0e-8, a_value * 1.0e-6)

        def integrand(x_value: float) -> float:
            return (
                1.0 / beta_a(x_value, n_f=n_f, loops=loops)
                + 1.0 / (beta0 * x_value * x_value)
                - (beta1 / (beta0 * beta0)) / x_value
            )

        running += _simpson(integrand, eps, a_value, panels)
    return mu_gev * math.exp(-0.5 * running)


def alpha_from_lambda_msbar(mu_gev: float, lambda_msbar_gev: float, n_f: int, loops: int = 4) -> float:
    if mu_gev <= 0.0 or lambda_msbar_gev <= 0.0:
        raise ValueError("mu_gev and lambda_msbar_gev must be > 0")
    lo = 1.0e-4
    hi = 1.0
    f_lo = lambda_msbar_from_alpha(mu_gev, lo, n_f=n_f, loops=loops) - lambda_msbar_gev
    f_hi = lambda_msbar_from_alpha(mu_gev, hi, n_f=n_f, loops=loops) - lambda_msbar_gev
    attempts = 0
    while f_lo * f_hi > 0.0 and attempts < 30:
        hi *= 1.5
        f_hi = lambda_msbar_from_alpha(mu_gev, hi, n_f=n_f, loops=loops) - lambda_msbar_gev
        attempts += 1
    if f_lo * f_hi > 0.0:
        raise RuntimeError("could not bracket alpha_s for Lambda inversion")
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        f_mid = lambda_msbar_from_alpha(mu_gev, mid, n_f=n_f, loops=loops) - lambda_msbar_gev
        if f_lo * f_mid > 0.0:
            lo = mid
            f_lo = f_mid
        else:
            hi = mid
            f_hi = f_mid
    return 0.5 * (lo + hi)


def build_artifact(d10: dict[str, object], forward: dict[str, object], loops: int) -> dict[str, object]:
    core = dict(d10["core_source"])
    mz_run = float(core["mz_run"])
    alpha3_mz = float(core["alpha3_mz"])
    singular_values_u = [float(value) for value in forward["singular_values_u"]]
    singular_values_d = [float(value) for value in forward["singular_values_d"]]
    m_c = singular_values_u[1]
    m_b = singular_values_d[2]
    if not (0.0 < m_c < m_b < mz_run):
        raise ValueError("expected m_c < m_b < m_z_run for the perturbative threshold chain")

    lambda5 = lambda_msbar_from_alpha(mz_run, alpha3_mz, n_f=5, loops=loops)
    alpha5_at_mb = alpha_from_lambda_msbar(m_b, lambda5, n_f=5, loops=loops)
    lambda4 = lambda_msbar_from_alpha(m_b, alpha5_at_mb, n_f=4, loops=loops)
    alpha4_at_mc = alpha_from_lambda_msbar(m_c, lambda4, n_f=4, loops=loops)
    lambda3 = lambda_msbar_from_alpha(m_c, alpha4_at_mc, n_f=3, loops=loops)

    return {
        "artifact": "oph_qcd_lambda_msbar3",
        "generated_utc": _timestamp(),
        "proof_status": "secondary_quantitative_descendant",
        "predictive_promotion_allowed": True,
        "source_artifacts": {
            "d10_observable_family": d10.get("artifact"),
            "forward_yukawas": forward.get("artifact"),
        },
        "law_id": "oph_qcd_lambda_msbar3_from_d10_alpha3",
        "borrowed_qft_piece": "perturbative_qcd_running_and_threshold_matching",
        "loops": loops,
        "source_inputs": {
            "mu_source_gev": mz_run,
            "alpha3_mz": alpha3_mz,
            "m_c_threshold_gev": m_c,
            "m_b_threshold_gev": m_b,
        },
        "threshold_chain": [
            {
                "n_f": 5,
                "mu_match_gev": mz_run,
                "alpha_s": alpha3_mz,
                "lambda_msbar_gev": lambda5,
                "role": "source_from_d10_alpha3",
            },
            {
                "n_f": 4,
                "mu_match_gev": m_b,
                "alpha_s": alpha5_at_mb,
                "lambda_msbar_gev": lambda4,
                "role": "bottom_threshold_continuity",
            },
            {
                "n_f": 3,
                "mu_match_gev": m_c,
                "alpha_s": alpha4_at_mc,
                "lambda_msbar_gev": lambda3,
                "role": "charm_threshold_continuity",
            },
        ],
        "outputs": {
            "Lambda_MSbar_5_gev": lambda5,
            "Lambda_MSbar_4_gev": lambda4,
            "Lambda_MSbar_3_gev": lambda3,
            "Lambda_MSbar_gev": lambda3,
        },
        "notes": [
            "This artifact is a secondary QCD descendant from the live OPH D10 alpha3 closure plus the current reference-free quark thresholds.",
            "It uses perturbative MS-bar running and simple threshold continuity at m_b and m_c.",
            "No hadron masses or PDG hadron inputs are used anywhere in this descendant.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the Lambda_MSbar^(3) descendant artifact.")
    parser.add_argument("--d10", default=str(DEFAULT_D10))
    parser.add_argument("--forward", default=str(DEFAULT_FORWARD))
    parser.add_argument("--loops", type=int, default=4)
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    d10 = json.loads(Path(args.d10).read_text(encoding="utf-8"))
    forward = json.loads(Path(args.forward).read_text(encoding="utf-8"))
    artifact = build_artifact(d10, forward, args.loops)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
