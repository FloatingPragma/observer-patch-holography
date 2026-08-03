"""Collar-matrix realization probe for the THERMO-REALIZATION receipt.

The receipt requires a source-derived collar transition matrix that
either equals the conditional-resampling kernel or passes
stochasticity, stationarity, and protected-charge preservation, with
detailed balance where microscopic reversibility is claimed. This
probe consumes the empirical repair transition matrix that the
simulator's own transition-clock builder counts from stored
observer histories of the earned 64k universe run, on the declared
visible-packet quotient. The matrix is counted from the run's
repair dynamics; no entry is constructed from the target kernel
formula. The probe measures the receipt properties on that matrix and
records the result with content pins to the source run; it does not
close the receipt.

Measured here, with the protected datum declared as the record
family: row stochasticity, off-fibre leakage of the protected datum,
the equal-fibre-row comparison in both the pairwise-row and the
stationary-profile form, the detailed-balance error of the
reversibilized chain, and relative-entropy descent to the stationary
law of the reversibilized chain. The raw chain's reducibility, which
the simulator's own eligibility gate names as a blocker, is inherited
and recorded.

Run with --write to refresh the committed probe receipt.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from conditional_repair_certificate import (
    ThermoError,
    canonical_json_bytes,
    require,
    tagged_sha256,
)

HERE = Path(__file__).resolve().parent
ARTIFACT_DIR = HERE / "runtime" / "realization_probe_64k"
MATRIX_PATH = ARTIFACT_DIR / "finite_repair_transition_matrix.npz"
REPORT_PATH = ARTIFACT_DIR / "finite_repair_transition_matrix_report.json"
PROBE_PATH = HERE / "runtime" / "collar_matrix_realization_probe.json"

# Content pins established at production time. The matrix and report
# were produced by the simulator CLI from the earned 64k run snapshot;
# the source pins identify the exact inputs consumed.
PINS: dict[str, str] = {
    "matrix_npz_sha256": (
        "sha256:314114f430a0d9a93b6a573ec61bebec213193d5"
        "ee7ff9a8bbbd7d57390b16bf"
    ),
    "report_json_sha256": (
        "sha256:01f880b8b1d8bc77ac497efb8c6bdb1c88c489df"
        "15b7528179f9bb04ac988c54"
    ),
    "source_observer_views_sha256": (
        "sha256:3568f030412ac3db259f438ff6c36e5cee51214d"
        "5004984569814905a6abec44"
    ),
    "source_manifest_sha256": (
        "sha256:11ba89c21f612da1579f4aa24242a4e0933b85b8"
        "c4c91c5a421065e82f312d24"
    ),
    "producer_checkout_commit": "4aa2ce703b5cd172f687c7c9bc3d2d4aa04ed11e",
    "source_run": (
        "oph-workspace://oph-physics-sim/data/earned_runs/"
        "oph_universe_64k_3p1d_reearned"
    ),
    "producer_command": (
        "oph-fpe finite-repair-transition-clock --run-dir "
        "data/earned_runs/oph_universe_64k_3p1d_reearned --out <out> "
        "--packet-fields "
        "record_family,checkpoint_class,s3_sector_class,"
        "repair_load_bucket"
    ),
}

FIBRE_FIELD = "record_family"
ROW_SUM_TOL = 1e-12
DB_TOL = 1e-12
KL_MONOTONE_TOL = 1e-12
KL_STEPS = 16


def file_sha256(path: Path) -> str:
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


def load_inputs() -> tuple[dict[str, Any], dict[str, Any]]:
    """Load the pinned artifacts.

    Every numerical result downstream is computed in fixed-order pure
    Python float arithmetic, with logarithms through mpmath, so the
    receipt is byte-reproducible across platforms; numpy touches the
    data only to decode the npz container.
    """
    require(MATRIX_PATH.exists(), "matrix artifact missing")
    require(REPORT_PATH.exists(), "report artifact missing")
    require(
        file_sha256(MATRIX_PATH) == PINS["matrix_npz_sha256"],
        "matrix artifact drifted from its pin",
    )
    require(
        file_sha256(REPORT_PATH) == PINS["report_json_sha256"],
        "report artifact drifted from its pin",
    )
    with np.load(MATRIX_PATH, allow_pickle=True) as z:
        data = {
            "counts": np.asarray(z["counts"], dtype=float).tolist(),
            "raw": np.asarray(z["raw_empirical"], dtype=float).tolist(),
            "rev": np.asarray(
                z["reversible_empirical"], dtype=float
            ).tolist(),
            "labels": [json.loads(str(s)) for s in z["state_labels"]],
        }
    report = json.loads(REPORT_PATH.read_text())
    return data, report


def fibre_of(label: list[list[Any]]) -> Any:
    for field, value in label:
        if field == FIBRE_FIELD:
            return value
    raise ThermoError(f"label misses the fibre field {FIBRE_FIELD}")


POWER_ITERATIONS = 2048
STATIONARY_RESIDUAL_TOL = 1e-12


def push_vec(mu: list[float], matrix: list[list[float]]) -> list[float]:
    n = len(mu)
    return [
        sum(mu[x] * matrix[x][y] for x in range(n)) for y in range(n)
    ]


def stationary_of(matrix: list[list[float]]) -> list[float]:
    """Stationary law by fixed-count power iteration in fixed-order
    pure Python arithmetic."""
    n = len(matrix)
    mu = [1.0 / n] * n
    for _ in range(POWER_ITERATIONS):
        mu = push_vec(mu, matrix)
        total = sum(mu)
        require(total > 0, "power iteration lost mass")
        mu = [v / total for v in mu]
    pushed = push_vec(mu, matrix)
    residual = max(abs(pushed[i] - mu[i]) for i in range(n))
    require(
        residual <= STATIONARY_RESIDUAL_TOL,
        "power iteration failed to reach the stationary law",
    )
    return mu


def kl(p: list[float], q: list[float], mp) -> Any:
    """Relative entropy with mpmath logarithms, so the value does not
    depend on the platform's libm."""
    out = mp.mpf(0)
    for pi_, qi in zip(p, q):
        if pi_ > 0:
            require(
                qi > 0,
                "stationary law misses support of the iterate",
            )
            out += mp.mpf(pi_) * mp.log(mp.mpf(pi_) / mp.mpf(qi))
    return out


def build_probe() -> dict[str, Any]:
    import mpmath

    data, report = load_inputs()
    counts = data["counts"]
    raw = data["raw"]
    rev = data["rev"]
    labels = data["labels"]
    n = len(labels)
    require(
        len(raw) == n and all(len(row) == n for row in raw),
        "matrix shape drift",
    )
    fibres = [fibre_of(lab) for lab in labels]
    row_mass = [sum(counts[i]) for i in range(n)]
    total_mass = sum(row_mass)
    require(total_mass > 0, "empty transition count table")
    visit = [m / total_mass for m in row_mass]

    # stochasticity of counted rows
    counted = [i for i in range(n) if row_mass[i] > 0]
    row_sum_err = max(
        abs(sum(raw[i]) - 1.0) for i in counted
    )
    require(row_sum_err <= ROW_SUM_TOL, "raw rows are not stochastic")

    # protected-datum leakage: off-fibre mass of counted rows
    off_masses = []
    for i in counted:
        off = sum(
            raw[i][j] for j in range(n) if fibres[j] != fibres[i]
        )
        off_masses.append((off, visit[i]))
    off_max = max(m for m, _ in off_masses)
    off_visit_weighted = sum(m * w for m, w in off_masses) / sum(
        w for _, w in off_masses
    )

    # equal-fibre-row comparison, pairwise form: within a fibre, the
    # conditional-resampling kernel has identical rows after
    # restriction to the fibre
    pair_tv_max = 0.0
    pair_tv_count = 0
    profile_tv_max = 0.0
    for fib in sorted({str(f) for f in fibres}):
        members = [
            i for i in counted if str(fibres[i]) == fib
        ]
        if len(members) < 2:
            continue
        restricted = []
        for i in members:
            block = [
                raw[i][j] for j in range(n) if str(fibres[j]) == fib
            ]
            s = sum(block)
            if s > 0:
                restricted.append([v / s for v in block])
        for a in range(len(restricted)):
            for b in range(a + 1, len(restricted)):
                tv = 0.5 * sum(
                    abs(x - y)
                    for x, y in zip(restricted[a], restricted[b])
                )
                pair_tv_max = max(pair_tv_max, tv)
                pair_tv_count += 1
        # profile form: restricted rows against the visit-weight
        # profile of the fibre
        fibre_visit = [
            visit[i] for i in range(n) if str(fibres[i]) == fib
        ]
        fv_total = sum(fibre_visit)
        if fv_total > 0 and restricted:
            profile = [v / fv_total for v in fibre_visit]
            for row in restricted:
                tv = 0.5 * sum(
                    abs(x - y) for x, y in zip(row, profile)
                )
                profile_tv_max = max(profile_tv_max, tv)

    # reversibilized chain: detailed balance and relative-entropy
    # descent to its stationary law
    rev_entry = report["matrices"]["reversible_empirical"]
    rev_summary = rev_entry.get("summary", rev_entry)
    require(
        bool(rev_summary["irreducible"]) and bool(rev_summary["aperiodic"]),
        "reversibilized chain lost ergodicity",
    )
    pi = stationary_of(rev)
    db_err = max(
        abs(pi[i] * rev[i][j] - pi[j] * rev[j][i])
        for i in range(n)
        for j in range(n)
    )
    require(db_err <= DB_TOL, "reversibilized detailed balance fails")
    mp = mpmath.mp
    saved_dps = mp.dps
    try:
        mp.dps = 30
        mu = [1.0 / n] * n
        kl_seq = [kl(mu, pi, mp)]
        for _ in range(KL_STEPS):
            mu = push_vec(mu, rev)
            kl_seq.append(kl(mu, pi, mp))
        descents = [
            kl_seq[t] - kl_seq[t + 1] for t in range(len(kl_seq) - 1)
        ]
        require(
            min(descents) >= -mp.mpf(str(KL_MONOTONE_TOL)),
            "relative entropy to the stationary law fails to descend",
        )
        kl_initial = mpmath.nstr(kl_seq[0], 12)
        kl_final = mpmath.nstr(kl_seq[-1], 12)
        kl_min_descent = mpmath.nstr(min(descents), 12)
    finally:
        mp.dps = saved_dps

    raw_entry = report["matrices"]["raw_empirical"]
    raw_summary = raw_entry.get("summary", raw_entry)
    inherited_blockers = list(report.get("blockers", []))

    body: dict[str, Any] = {
        "schema": "oph.collar_matrix_realization_probe.v1",
        "status": (
            "MEASURED_PROBE__SOURCE_PRODUCED_MATRIX_ATTAINED__"
            "PROTECTED_DATUM_LEAKAGE_AND_FIBRE_PROFILES_MEASURED__"
            "REVERSIBILIZED_KL_DESCENT_VERIFIED__"
            "RAW_CHAIN_REDUCIBLE__RECEIPT_OPEN"
        ),
        "receipt_target": "THERMO-REALIZATION",
        "noncircularity": (
            "the transition matrix is counted by the simulator's "
            "transition-clock builder from stored observer repair "
            "histories of the earned run; no entry is constructed "
            "from the conditional-resampling formula, and the "
            "comparison profiles are built from the chain's own "
            "visit weights"
        ),
        "pins": PINS,
        "quotient": {
            "packet_fields": report["packet_fields"],
            "fibre_field": FIBRE_FIELD,
            "state_count": n,
            "counted_state_count": len(counted),
            "observer_count": report["observer_count"],
            "transition_count": report["transition_count"],
            "weight_field": report["weight_field"],
        },
        "measurements": {
            "row_sum_max_err": row_sum_err,
            "off_fibre_mass_max": off_max,
            "off_fibre_mass_visit_weighted": off_visit_weighted,
            "equal_fibre_row_pairwise_tv_max": pair_tv_max,
            "equal_fibre_row_pairwise_pairs": pair_tv_count,
            "fibre_profile_tv_max": profile_tv_max,
            "reversibilized_detailed_balance_max_err": db_err,
            "reversibilized_spectral_gap": (
                1.0 - float(rev_summary["lambda_2"])
            ),
            "kl_to_stationary_initial": kl_initial,
            "kl_to_stationary_final": kl_final,
            "kl_steps": KL_STEPS,
            "kl_min_stepwise_descent": kl_min_descent,
            "raw_chain_irreducible": bool(raw_summary["irreducible"]),
            "raw_chain_aperiodic": bool(raw_summary["aperiodic"]),
        },
        "inherited_blockers": inherited_blockers,
        "verdict": {
            "stochasticity": "pass",
            "reversibilized_detailed_balance": "pass",
            "reversibilized_kl_descent": "pass",
            "protected_datum_preservation": (
                "measured leakage; the off-fibre mass quantifies how "
                "far one repair step moves the record family on this "
                "quotient"
            ),
            "equal_fibre_row": (
                "measured deviation; the kernel is Metropolis repair, "
                "and coincidence with conditional resampling is a "
                "hypothesis to test, never an assumption"
            ),
            "receipt_state": (
                "open; the raw chain is reducible through freezeout "
                "absorption, so a certified irreducible primary chain "
                "or an explicit recurrent-class restriction stays "
                "required"
            ),
        },
    }
    body["receipt_sha256"] = tagged_sha256(canonical_json_bytes(body))
    return body


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args(argv)
    probe = build_probe()
    if args.write:
        PROBE_PATH.parent.mkdir(parents=True, exist_ok=True)
        PROBE_PATH.write_text(
            canonical_json_bytes(probe).decode() + "\n"
        )
        print(f"wrote {PROBE_PATH}")
    print(probe["status"])
    print(
        "off-fibre visit-weighted mass:",
        probe["measurements"]["off_fibre_mass_visit_weighted"],
    )
    print(
        "pairwise fibre-row TV max:",
        probe["measurements"]["equal_fibre_row_pairwise_tv_max"],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
