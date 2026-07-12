#!/usr/bin/env python3
"""Emit the lane status for the real diagnostic lattice engine (#425).

Records what now exists (a physics-true engine with executed output), what
its verification anchors are, and exactly what the #425 production closure
still requires. The status artifact never promotes the diagnostic output.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / "particles" / "runs" / "hadron" / "lattice_engine_lane_status.json"
EXPORT = ROOT / "particles" / "runs" / "hadron" / "lattice_diagnostic_backend_export.json"


def _now_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_status() -> dict[str, Any]:
    export: dict[str, Any] | None = None
    if EXPORT.exists():
        export = json.loads(EXPORT.read_text(encoding="utf-8"))

    executed = {
        "export_present": export is not None,
    }
    if export is not None:
        executed.update({
            "execution_class": export["execution_class"],
            "ensemble_id": export["ensemble"]["ensemble_id"],
            "n_configs": export["ensemble"]["n_configs"],
            "plaquette_mean": export["ensemble"]["plaquette_mean"],
            "plaquette_agrees_with_literature": export["ensemble"][
                "plaquette_agrees_with_literature"],
            "channels": export["channels"]["exported"],
            "masses": [
                {
                    "kappa": block["kappa"],
                    "am_pi": block["pi_iso"]["am_plateau"],
                    "am_pi_error": block["pi_iso"]["am_plateau_jackknife_error"],
                    "am_N": block["N_iso"]["am_plateau"],
                    "am_N_error": block["N_iso"]["am_plateau_jackknife_error"],
                }
                for block in export["analysis"]
            ],
            "hmc_acceptance": export["dynamical_branch_validation"]["acceptance_rate"],
            "hmc_exp_minus_dh": export["dynamical_branch_validation"][
                "exp_minus_delta_h_mean"],
        })

    return {
        "artifact": "oph_lattice_engine_lane_status",
        "generated_utc": _now_utc(),
        "github_issue": 425,
        "status": "real_engine_executed_diagnostic_scale_production_closure_open",
        "engine": {
            "package": "particles/hadron/lattice_backend/",
            "runner": "particles/hadron/run_lattice_diagnostic_backend.py",
            "tests": "particles/hadron/test_lattice_backend.py",
            "capabilities": [
                "quenched SU(3) heatbath/overrelaxation ensemble generation",
                "clover-Wilson valence Dirac operator with CG propagators",
                "pi_iso and N_iso direct/exchange contractions",
                "two-flavor Wilson pseudofermion HMC (leapfrog, Metropolis)",
            ],
            "verification_anchors": [
                "gamma algebra and gamma5-hermiticity exact to double precision",
                "plaquette and correlator gauge invariance under random rotations",
                "free-field pion and nucleon effective masses at 2x and 3x the "
                "Wilson quark pole mass",
                "gauge and pseudofermion forces match finite differences",
                "leapfrog dH scales as dt^2 and integration is exactly reversible",
                "quenched plaquette at beta=5.7 matches the literature value 0.549",
            ],
        },
        "distinction_from_surrogate": (
            "run_local_diagnostic_backend.py synthesizes correlators from target "
            "masses and is execution-bridge scaffolding only. The lattice engine "
            "computes correlators from gauge dynamics with no target anchoring; "
            "its masses are whatever the bare parameters produce."
        ),
        "executed_output": executed,
        "issue_425_remaining": {
            "closure_requires": [
                "seeded N_f = 2+1 ensembles at the frozen receipt schedule "
                "(N_therm = 2048, N_sep = 512) on production volumes",
                "strange RHMC branch with persisted rational coefficients",
                "clover force in the dynamical update",
                "even-odd preconditioned solvers at production tolerances",
                "continuum, finite-volume, chiral, and statistical budgets",
                "publication-complete manifest provenance",
            ],
            "compute_class": "HPC allocation or an established production engine "
                             "driven through run_external_production_backend.py; "
                             "out of local scope",
            "promotion_gate": "unchanged: hadron rows stay suppressed until "
                              "production output and systematics exist",
        },
        "guards": {
            "diagnostic_output_promotable": False,
            "issue_425_closed_by_this_artifact": False,
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Emit the lattice engine lane status.")
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()
    payload = build_status()
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
