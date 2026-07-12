#!/usr/bin/env python3
"""Run the real lattice engine at diagnostic scale and export correlators.

This is genuine lattice QCD execution: quenched SU(3) ensemble generation
with the Wilson plaquette action, clover-Wilson valence propagators from CG
solves, pion and nucleon correlators from Wick contractions, and a small
two-flavor HMC validation chain for the dynamical branch. Nothing in the
output is anchored to target masses; the masses are whatever the lattice
dynamics produce at the chosen bare parameters.

Execution class: real_lattice_diagnostic_toy_scale. The export is
non-promoting by construction: the #425 production closure needs the frozen
seeded 2+1 family at HPC scale plus continuum, finite-volume, chiral, and
statistical budgets, none of which a laptop-scale run can supply. This
runner exists so the hadron lane has a physics-true engine and an executed
end-to-end path instead of a target-anchored surrogate.

Run:
    python3 code/particles/hadron/run_lattice_diagnostic_backend.py [--smoke]
writes code/particles/runs/hadron/lattice_diagnostic_backend_export.json.
"""

from __future__ import annotations

import argparse
import json
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from lattice_backend.core import average_plaquette, cold_start, sweep  # noqa: E402
from lattice_backend.dirac import WilsonClover, point_propagator  # noqa: E402
from lattice_backend.hmc import TwoFlavorHMC  # noqa: E402
from lattice_backend.spectroscopy import (  # noqa: E402
    effective_mass,
    nucleon_correlators,
    pion_correlator,
)

OUT_PATH = HERE.parent / "runs" / "hadron" / "lattice_diagnostic_backend_export.json"
CHECKPOINT = HERE.parent / "runs" / "hadron" / "lattice_diagnostic_backend_checkpoint.npz"

PRODUCTION_PARAMS = {
    "shape": (12, 6, 6, 6),
    "beta": 5.7,
    "kappas": [0.145, 0.150],
    "c_sw": 1.0,
    "n_therm": 200,
    "n_sep": 20,
    "n_cfg": 16,
    "n_or": 1,
    "cg_tol": 1e-7,
    "seed": 20260712,
    "plateau_window": [3, 5],
    "hmc": {"shape": (4, 4, 4, 4), "beta": 5.5, "kappa": 0.13,
            "n_steps": 48, "n_therm": 6, "n_meas": 12, "seed": 424242},
}

SMOKE_PARAMS = {
    "shape": (8, 4, 4, 4),
    "beta": 5.7,
    "kappas": [0.150],
    "c_sw": 1.0,
    "n_therm": 20,
    "n_sep": 5,
    "n_cfg": 2,
    "n_or": 1,
    "cg_tol": 1e-6,
    "seed": 99,
    "plateau_window": [2, 3],
    "hmc": {"shape": (4, 2, 2, 2), "beta": 5.5, "kappa": 0.12,
            "n_steps": 10, "n_therm": 1, "n_meas": 3, "seed": 4242},
}

PLAQUETTE_LITERATURE_B57 = 0.549


def _now_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _git_commit() -> str:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--short=12", "HEAD"],
            cwd=HERE, check=True, capture_output=True, text=True)
        return result.stdout.strip() or "unknown"
    except Exception:
        return "unknown"


def _jackknife_mean_err(samples: np.ndarray) -> tuple[float, float]:
    n = len(samples)
    if n < 2:
        return float(np.mean(samples)), float("nan")
    means = np.array([np.mean(np.delete(samples, i)) for i in range(n)])
    center = float(np.mean(means))
    err = float(np.sqrt((n - 1) / n * np.sum((means - center) ** 2)))
    return center, err


def generate_ensemble(params: dict, log) -> tuple[list[np.ndarray], list[float]]:
    rng = np.random.default_rng(params["seed"])
    u = cold_start(params["shape"])
    log(f"thermalizing {params['n_therm']} sweeps at beta={params['beta']}")
    for _ in range(params["n_therm"]):
        sweep(rng, u, params["beta"], n_or=params["n_or"])
    configs, plaquettes = [], []
    for i in range(params["n_cfg"]):
        for _ in range(params["n_sep"]):
            sweep(rng, u, params["beta"], n_or=params["n_or"])
        p = average_plaquette(u)
        configs.append(u.copy())
        plaquettes.append(p)
        log(f"cfg {i:03d}: plaquette = {p:.5f}")
    return configs, plaquettes


def measure_config(u: np.ndarray, params: dict) -> dict:
    out = {}
    shape = params["shape"]
    for kappa in params["kappas"]:
        op = WilsonClover(u, kappa=kappa, c_sw=params["c_sw"])
        prop, info = point_propagator(op, shape, tol=params["cg_tol"])
        c_pi = pion_correlator(prop)
        n_dir, n_exc = nucleon_correlators(prop)
        out[f"kappa_{kappa:.3f}"] = {
            "pi_iso": c_pi.tolist(),
            "N_iso_direct": n_dir.tolist(),
            "N_iso_exchange": n_exc.tolist(),
            "max_cg_relative_residual": float(max(info["cg_relative_residuals"])),
            "max_cg_iterations": int(max(info["cg_iterations"])),
        }
    return out


def dynamical_branch_validation(params: dict, log) -> dict:
    """Small two-flavor HMC chain: acceptance and the Creutz equality
    <exp(-dH)> = 1, the direct stochastic test of the dynamical branch."""
    hp = params["hmc"]
    rng = np.random.default_rng(hp["seed"])
    u = cold_start(hp["shape"])
    hmc = TwoFlavorHMC(beta=hp["beta"], kappa=hp["kappa"],
                       n_steps=hp["n_steps"], cg_tol=1e-10)
    for i in range(hp["n_therm"]):
        u, res = hmc.trajectory(rng, u)
        log(f"hmc therm {i}: dH={res.delta_h:+.4f} accept={res.accepted}")
    dhs, accepts = [], []
    for i in range(hp["n_meas"]):
        u, res = hmc.trajectory(rng, u)
        dhs.append(res.delta_h)
        accepts.append(res.accepted)
        log(f"hmc meas {i}: dH={res.delta_h:+.4f} accept={res.accepted}")
    exp_dh = np.exp(-np.array(dhs))
    center, err = _jackknife_mean_err(exp_dh)
    return {
        "branch": "two_flavor_wilson_pseudofermion_hmc",
        "parameters": {k: (list(v) if isinstance(v, tuple) else v) for k, v in hp.items()},
        "trajectories_measured": len(dhs),
        "acceptance_rate": float(np.mean(accepts)),
        "delta_h_values": [float(x) for x in dhs],
        "exp_minus_delta_h_mean": center,
        "exp_minus_delta_h_jackknife_error": err,
        "creutz_equality_within_2_sigma": bool(abs(center - 1.0) <= 2.0 * err) if err == err else None,
        "note": "the strange RHMC branch and the clover force are production "
                "steps outside this toy dynamical demo, per the engine docstring",
    }


def channel_analysis(per_cfg: list[dict], kappa: float, params: dict) -> dict:
    key = f"kappa_{kappa:.3f}"
    pi = np.array([c[key]["pi_iso"] for c in per_cfg])
    nd = np.array([c[key]["N_iso_direct"] for c in per_cfg])
    ne = np.array([c[key]["N_iso_exchange"] for c in per_cfg])
    n_iso = nd - ne
    lo, hi = params["plateau_window"]

    def plateau(corr_samples: np.ndarray) -> dict:
        n = len(corr_samples)
        estimates = []
        for i in range(n):
            mean_corr = np.mean(np.delete(corr_samples, i, axis=0), axis=0)
            am = effective_mass(mean_corr)
            estimates.append(np.nanmean(am[lo:hi + 1]))
        estimates = np.array(estimates)
        center = float(np.mean(estimates))
        err = float(np.sqrt((n - 1) / n * np.sum((estimates - center) ** 2)))
        mean_corr = np.mean(corr_samples, axis=0)
        return {
            "correlator_ensemble_mean": mean_corr.tolist(),
            "effective_mass": [None if np.isnan(x) else float(x)
                               for x in effective_mass(mean_corr)],
            "am_plateau": center,
            "am_plateau_jackknife_error": err,
            "plateau_window": [lo, hi],
        }

    return {
        "kappa": kappa,
        "pi_iso": plateau(pi),
        "N_iso": plateau(n_iso),
        "n_configs": len(per_cfg),
    }


def build_export(params: dict, smoke: bool) -> dict:
    t_start = time.time()
    log_lines: list[str] = []

    def log(msg: str) -> None:
        stamp = f"[{time.time() - t_start:8.1f}s] {msg}"
        log_lines.append(stamp)
        print(stamp, flush=True)

    configs, plaquettes = generate_ensemble(params, log)
    plaq_mean = float(np.mean(plaquettes))

    per_cfg = []
    for i, u in enumerate(configs):
        t0 = time.time()
        per_cfg.append(measure_config(u, params))
        log(f"cfg {i:03d}: propagators + contractions in {time.time() - t0:.0f}s")
        np.savez(CHECKPOINT, done=i + 1)

    analyses = [channel_analysis(per_cfg, kappa, params) for kappa in params["kappas"]]
    hmc_block = dynamical_branch_validation(params, log)

    ensemble_id = ("quenched_wilson_b{beta}_L{sx}T{st}_diag".format(
        beta=str(params["beta"]).replace(".", "p"),
        sx=params["shape"][1], st=params["shape"][0]))

    return {
        "artifact": "oph_lattice_diagnostic_backend_export",
        "format_version": 1,
        "generated_utc": _now_utc(),
        "execution_class": "real_lattice_diagnostic_toy_scale",
        "smoke_mode": smoke,
        "guards": {
            "real_lattice_execution": True,
            "target_anchored": False,
            "surrogate_hadron_artifact": False,
            "production_execution_class": False,
            "public_promotion_allowed": False,
            "satisfies_issue_425_closure": False,
        },
        "backend": {
            "family": "oph_local_lattice_engine",
            "name": "lattice_backend",
            "version": "1.0.0",
            "git_commit": _git_commit(),
            "run_id": f"lattice_diag_{params['seed']}",
            "build_id": "in_repo_numpy",
            "machine": f"{platform.node()} {platform.machine()} python{platform.python_version()}",
        },
        "physics_declaration": {
            "gauge_action": "isotropic Wilson plaquette, SU(3), periodic",
            "ensemble": "quenched (pure gauge) generation via Cabibbo-Marinari "
                        "heatbath with overrelaxation",
            "valence_fermions": "clover-improved Wilson, c_SW = 1.0, "
                                "antiperiodic time boundary",
            "dynamical_branch": "two-flavor unimproved Wilson pseudofermion HMC, "
                                "validated on a small chain; quenched ensembles "
                                "carry the spectroscopy in this run",
            "deviations_from_production_spec": [
                "quenched instead of seeded 2+1 dynamical ensembles",
                "no strange RHMC branch executed",
                "clover force absent from the toy dynamical branch",
                "toy volume and single coarse lattice spacing",
                "no continuum, finite-volume, chiral, or full statistical budgets",
            ],
        },
        "ensemble": {
            "ensemble_id": ensemble_id,
            "lattice_shape_TXYZ": list(params["shape"]),
            "beta": params["beta"],
            "n_therm_sweeps": params["n_therm"],
            "n_sep_sweeps": params["n_sep"],
            "n_configs": params["n_cfg"],
            "rng_seed": params["seed"],
            "plaquette_per_config": [float(p) for p in plaquettes],
            "plaquette_mean": plaq_mean,
            "plaquette_literature_beta57": PLAQUETTE_LITERATURE_B57,
            "plaquette_agrees_with_literature": bool(abs(plaq_mean - PLAQUETTE_LITERATURE_B57) < 0.02),
        },
        "channels": {
            "exported": ["pi_iso", "N_iso_direct", "N_iso_exchange"],
            "combination_rule": "N_iso = N_iso_direct - N_iso_exchange",
        },
        "correlators_per_config": per_cfg,
        "analysis": analyses,
        "dynamical_branch_validation": hmc_block,
        "runtime_seconds": round(time.time() - t_start, 1),
        "log": log_lines,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the diagnostic lattice backend.")
    parser.add_argument("--smoke", action="store_true", help="tiny fast run for tests")
    parser.add_argument("--output", default=str(OUT_PATH))
    args = parser.parse_args()
    params = SMOKE_PARAMS if args.smoke else PRODUCTION_PARAMS

    payload = build_export(params, smoke=args.smoke)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=1)
        f.write("\n")
    print(f"\nwrote {out}")
    for block in payload["analysis"]:
        pi_b, n_b = block["pi_iso"], block["N_iso"]
        print(f"kappa={block['kappa']}: am_pi = {pi_b['am_plateau']:.4f}"
              f" +- {pi_b['am_plateau_jackknife_error']:.4f}"
              f"   am_N = {n_b['am_plateau']:.4f}"
              f" +- {n_b['am_plateau_jackknife_error']:.4f}")
    hb = payload["dynamical_branch_validation"]
    print(f"HMC: acceptance {hb['acceptance_rate']:.2f}, "
          f"<exp(-dH)> = {hb['exp_minus_delta_h_mean']:.3f}"
          f" +- {hb['exp_minus_delta_h_jackknife_error']:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
