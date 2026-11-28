"""Determine whether a lottery of native probe flights can maintain the ordinary equilibrium."""

import argparse
from hashlib import sha256
import json
from os.path import relpath
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

from . import observation_density as obs

HERE, lab = obs.HERE, obs.lab
ARMS = ("dense", "sparse", "left", "right")
SEED = 28112113


def short_steps(h, start, end):
    # Each dyadic substep has ||abs(H dt)||_2 <= 1 by Hermiticity and the
    # row-sum bound. Evaluate its Taylor series through degree 26. Repeated
    # substeps avoid accumulating a phase discrepancy across the two arms.
    substeps = 2 ** max(0, int(np.ceil(np.log2(abs(h).sum(axis=1).max()))))
    A = h / substeps
    degree = 26
    row_terms = np.max(np.diff(A.tocsr().indptr))
    # A sparse product and real division incur gamma_(2d+2). In the sum
    # of powers their first-order amplification is at most e. Accumulating
    # 27 terms costs gamma_27. This allowance also bounds higher-order
    # products of those rounding factors for the finite dimensions here.
    step_error = 16 * np.e * (row_terms + 32) * np.finfo(float).eps + 1e-26
    # The four initial columns have Frobenius norm two. Composition of
    # the two flight halves and tracing the pointer cost less than 24
    # times the absolute column-matrix error in any unit head observable.
    flight_error = 2 * np.expm1(step_error * substeps * end)
    assert 24 * flight_error + 1e-11 < 3e-8
    state = start.astype(complex)
    out = [state]
    for _ in range(end):
        for _ in range(substeps):
            term = state.copy()
            total = state.copy()
            for k in range(1, degree + 1):
                term = (A @ term) / k
                total += term
            state = total
        out.append(state)
    return np.asarray(out)


def center(record):
    matrices, err = obs.reg.src.observed(record, lab.SHOTS, False)
    coord = np.einsum("kij,bji->bk", obs.pos.AXES[1:], matrices).real.ravel()
    return coord, np.linalg.norm(2 * err) + 1e-12


def separator(points):
    dim = points.shape[1]
    opt = linprog(
        np.r_[np.zeros(dim), -1],
        A_ub=np.column_stack((-points, np.ones(len(points)))),
        b_ub=np.zeros(len(points)),
        bounds=[(-1, 1)] * dim + [(0, None)],
        method="highs",
    )
    if not opt.success:
        raise RuntimeError(opt.message)
    w = opt.x[:-1]
    n = np.linalg.norm(w)
    if n <= 1e-8:
        return np.zeros(dim)
    return np.rint(w / n * 2**26) / 2**26


def interior(points, radius):
    # If the true columns can realize each signed coordinate vertex up to e,
    # their hull contains a ball of radius min(t_j)/sqrt(dim)-e.
    count, dim = points.shape
    trials = []
    for j in range(dim):
        for sign in (-1, 1):
            ray = sign * np.eye(dim)[j]
            equal = np.vstack(
                (np.column_stack((points.T, -ray)), np.r_[np.ones(count), 0])
            )
            opt = linprog(
                np.r_[np.zeros(count), -1],
                A_eq=equal,
                b_eq=np.r_[np.zeros(dim), 1],
                bounds=(0, None),
                method="highs",
            )
            if not opt.success:
                return {"status": "unresolved", "failed_axis": [j, sign]}
            weights = np.maximum(0, opt.x[:-1])
            total = weights.sum()
            if not np.isfinite(opt.x).all() or total <= 0:
                return {"status": "unresolved", "failed_axis": [j, sign]}
            weights /= total
            t = opt.x[-1]
            residual = (
                np.linalg.norm(np.einsum("i,ij->j", weights, points) - t * ray) + 1e-11
            )
            if not np.isfinite(residual):
                return {"status": "unresolved", "failed_axis": [j, sign]}
            trials.append(
                {
                    "axis": j,
                    "sign": sign,
                    "distance": t,
                    "residual": residual,
                    "weights": [[int(i), weights[i]] for i in np.flatnonzero(weights)],
                }
            )
    inner = (
        min(z["distance"] for z in trials) / np.sqrt(dim)
        - max(z["residual"] for z in trials)
        - radius
        - 1e-11
    )
    return {
        "status": "contains_reference" if inner > 0 else "unresolved",
        "true_inradius_lower": inner,
        "axis_witnesses": trials,
    }


def classify(points, radius):
    w = separator(points)
    margin = np.min(points @ w) - np.linalg.norm(w) * radius - 1e-11
    if margin > 0:
        return {
            "status": "excluded",
            "signed_Pauli_witness": w.tolist(),
            "minimum_true_drift": margin,
        }
    return interior(points, radius)


def run():
    rng = np.random.default_rng(SEED)
    _, sources, *_ = obs.c.trials()
    reference = sources["reference"]
    head = lab.pairs(reference, rng, False)
    x, ex = center(head)
    ports, _ = lab.six.survey.return_bins()
    tape = {}
    panel = []
    errors = []
    for port in ports:
        for sign in (-1, 1):
            _, _, local = obs.history(
                reference,
                (port,),
                rng,
                {"left": 1j, "right": -1j},
                signs=(sign,),
                windows=tuple(range(1, 13)),
                propagate=short_steps,
            )
            for tau, arms in local.items():
                key = f"{port}/{sign}/{tau}"
                counts = {k: lab.pairs(arms[k].sum(0), rng, False) for k in ARMS}
                tape[key] = counts
                obs_arm = [center(counts[k]) for k in ARMS]
                panel.append([z - x for z, _ in obs_arm])
                errors.extend(e + 3e-8 for _, e in obs_arm)
    panel = np.array(panel)
    certificate = Path(obs.c.dose.enclosure.__file__).with_name("certificate.json")
    h = json.loads(certificate.read_bytes())["state_HS_error_upper"] + 1e-12
    # One prepared copy, then one CPTP output. The true X_* lies h from
    # the original preparation, independently of the selected flight.
    stationary = (1 + np.sqrt(8)) * h
    radius = max(errors) + ex + stationary + 1e-11
    dense, sparse, left, right = panel.transpose(1, 0, 2)
    tested = {
        "dense": classify(dense, radius),
        "real_zero_to_one": classify(panel[:, :2].reshape(-1, 30), radius),
    }
    # Circumscribed regular polygon offers directions for a separating read.
    # Verify the proposal on the whole disk analytically, including angles
    # lying between polygon vertices.
    real, imag = sparse - dense, (left - right) / 2

    def polygon(count, scale):
        angle = 2 * np.pi * np.arange(count) / count
        return (
            dense[:, None]
            + scale
            * (
                np.cos(angle)[None, :, None] * real[:, None]
                + np.sin(angle)[None, :, None] * imag[:, None]
            )
        ).reshape(-1, 30)

    w = separator(polygon(16, 1 / np.cos(np.pi / 16)))
    window_min = dense @ w - np.hypot(real @ w, imag @ w)
    disk_radius = 4 * max(errors) + ex + stationary + 2e-11
    margin = min(window_min) - np.linalg.norm(w) * disk_radius - 2e-11
    tested["complex_unit_disk"] = (
        {
            "status": "excluded",
            "signed_Pauli_witness": w.tolist(),
            "minimum_true_drift": margin,
            "per_flight_lower_at_prepared_center": window_min.tolist(),
        }
        if margin > 0
        else interior(polygon(16, 1 - 1e-12), disk_radius)
    )
    return {
        "seed": SEED,
        "shots_per_head_setting": lab.SHOTS,
        "source_heads": head,
        "terminal_heads": tape,
        "settings": list(tape),
        "arm_order": list(ARMS),
        "geometry_reference_HS_radius": h,
        "long_flight_numerical_allowance": 3e-8,
        "true_drift_coordinate_error": radius,
        "whole_pointer_disk_error": disk_radius,
        "outcomes": tested,
        "inputs": {
            relpath(f, HERE): sha256(f.read_bytes()).hexdigest()
            for f in (
                HERE / "EQUILIBRIUM_SERVICE_PROTOCOL.md",
                HERE / "finite_pairing_dose.json",
                certificate,
            )
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    r = run()
    path = HERE / "equilibrium_service.json"
    if args.write:
        path.write_text(json.dumps(r, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), r)
    print(
        {
            k: [v["status"], v.get("minimum_true_drift"), v.get("true_inradius_lower")]
            for k, v in r["outcomes"].items()
        }
    )
