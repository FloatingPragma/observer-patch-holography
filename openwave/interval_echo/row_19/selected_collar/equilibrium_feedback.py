"""Calibrate native feedback using the actual click or miss in each serviced copy."""

import argparse
from hashlib import sha256
import json

import numpy as np
from scipy.optimize import linprog

from . import equilibrium_service as service

obs, lab, HERE = service.obs, service.lab, service.HERE
SEED = 28112115
N = lab.SHOTS
TIMES = tuple(range(13))


def head_effects():
    axes = lab.heads.AXES
    return np.array(
        [
            [
                np.kron((np.eye(2) + s * a) / 2, (np.eye(2) + t * b) / 2)
                for s in lab.heads.SIGNS
                for t in lab.heads.SIGNS
            ]
            for a in axes[1:]
            for b in axes[1:]
        ]
    )


def acquire():
    rng = np.random.default_rng(SEED)
    _, sources, *_ = obs.c.trials()
    x = sources["reference"]
    reference = np.array(lab.pairs(x, rng, False))
    ports, _ = lab.six.survey.return_bins()
    native = lab.six.survey.native
    board = native.iso.pairing.formation.federation(1)
    region = native.iso.pairing.formation.region(board)
    effects = head_effects()
    records = np.zeros((len(ports), 2, 12, 13, 3, 3, 2, 8), dtype=np.int64)
    for i, home in enumerate(ports):
        for q, sign in enumerate((-1, 1)):
            columns = []
            for b in (0, 1):
                h = native.generator(board, region, sign, bool(b)) * (-0.25j)
                initial = np.zeros((h.shape[0], 4))
                initial[4 * home : 4 * home + 4] = np.eye(4)
                columns.append(
                    service.short_steps(h, initial, 24).reshape(25, -1, 4, 4)
                )
            for t in range(1, 13):
                for end in TIMES:
                    arms = np.empty((2, 2, 4, 4), complex)
                    for b, path in enumerate(columns):
                        click = path[end] @ path[t, home]
                        miss = path[t + end] - click
                        for j, K in enumerate((click, miss)):
                            arms[j, b] = (K @ x[b] @ K.conj().swapaxes(-1, -2)).sum(
                                axis=0
                            )
                    prob = np.einsum("akij,zbji->azbk", effects, arms).real.reshape(
                        9, 16
                    )
                    for setting, law in enumerate(prob):
                        # Each cell is an actual joint midpoint and end reading.
                        assert min(law) > -1e-13 and abs(law.sum() - 1) < 1e-11
                        order = np.argsort(law, kind="stable")
                        p = np.maximum(0, law[order])
                        p[-1] = 1 - p[:-1].sum()
                        counts = rng.multinomial(2**52, p, size=N // 2**52).sum(axis=0)
                        a, c = divmod(setting, 3)
                        records[i, q, t - 1, end, a, c] = counts[
                            np.argsort(order)
                        ].reshape(2, 8)
    return reference, records


def drift_panel(reference, records):
    x, _ = service.center(reference)
    reports = records.reshape(144, 13, 3, 3, 2, 8)
    assert np.min(reports) >= 0 and np.all(reports.sum(axis=(-1, -2)) == N)
    branch = np.empty((144, 13, 2, 30))
    for history in range(144):
        for end in TIMES:
            for arm in (0, 1):
                C = reports[history, end, :, :, arm]
                tape = np.concatenate((C, (N - C.sum(axis=-1))[..., None]), axis=-1)
                matrices, _ = obs.reg.src.observed(tape, N, True)
                branch[history, end, arm] = np.einsum(
                    "kij,bji->bk", obs.pos.AXES[1:], matrices
                ).real.ravel()
    return (branch[:, :, 0, None] + branch[:, None, :, 1] - x).reshape(-1, 30)


def inner_certificate(points, radius):
    # Start with a deterministic subset to keep exact finite witnesses small.
    # Every accepted subset is contained in the full policy hull. If it does
    # not decide the question the separator and final attempt use all rules.
    for indices in (np.arange(0, len(points), 11), np.arange(len(points))):
        witness = service.interior(points[indices], radius)
        if witness["status"] == "contains_reference":
            for vertex in witness["axis_witnesses"]:
                vertex["weights"] = [[int(indices[i]), w] for i, w in vertex["weights"]]
            return witness
    return witness


def decide(points, radius):
    w = service.separator(points)
    margin = np.min(points @ w) - radius * np.linalg.norm(w) - 1e-11
    if margin > 0:
        return {
            "status": "excluded",
            "signed_Pauli_witness": w.tolist(),
            "minimum_drift": margin,
        }
    cert = inner_certificate(points, radius)
    if cert["status"] == "unresolved":
        # Record a nominal mixture with the smallest worst coordinate drift.
        # This preserves a diagnostic without treating central feasibility
        # as proof of an exactly stationary channel.
        n, d = points.shape
        A = np.column_stack((np.vstack((points.T, -points.T)), -np.ones(2 * d)))
        fit = linprog(
            np.r_[np.zeros(n), 1],
            A_ub=A,
            b_ub=np.zeros(2 * d),
            A_eq=[np.r_[np.ones(n), 0]],
            b_eq=[1],
            bounds=(0, None),
            method="highs",
        )
        if fit.success:
            weights = fit.x[:-1]
            cert["nominal_residual_HS"] = float(np.linalg.norm(weights @ points))
            cert["nominal_weights"] = [
                [int(i), float(weights[i])] for i in np.flatnonzero(weights)
            ]
    return cert


def analyze(reference, counts):
    src = HERE / "equilibrium_service.json"
    h = json.loads(src.read_bytes())["geometry_reference_HS_radius"]
    panel = drift_panel(reference, counts)
    # Each HS coordinate is a registered product-Pauli score divided by two,
    # with zero on any other branch/parity. Its single-shot range has length
    # one. The union bound includes every reference and branch coordinate.
    coordinate_error = np.sqrt(np.log(2 * 500000 / 0.01) / (2 * N)) + 2e-13
    base_error = 3 * np.sqrt(30) * coordinate_error + 2 * 3e-8 + 1e-11
    radius = base_error + (1 + np.sqrt(8)) * h
    groups = {"blind": panel, "even": panel[:, :15], "odd": panel[:, 15:]}
    return {
        "seed": SEED,
        "shots": N,
        "deterministic_policy_count": len(panel),
        "rule_index": "(((sender_index*2+sign_index)*12+midpoint-1)*13+click_duration)*13+miss_duration",
        "source_HS_localization": h,
        "per_read_coordinate_error": coordinate_error,
        "drift_HS_radius": radius,
        "uncertain_output_HS_at_prepared_reference": base_error,
        "outcomes": {k: decide(z, radius) for k, z in groups.items()},
        "inputs": {
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                src,
                HERE / "EQUILIBRIUM_FEEDBACK_PROTOCOL.md",
                HERE / "finite_pairing_dose.json",
            )
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    reference, counts = acquire()
    raw = HERE / "equilibrium_feedback_counts.npz"
    if args.write:
        np.savez_compressed(raw, reference=reference, terminal=counts)
    else:
        with np.load(raw, allow_pickle=False) as old:
            np.testing.assert_array_equal(old["reference"], reference)
            np.testing.assert_array_equal(old["terminal"], counts)
    result = analyze(reference, counts)
    result["inputs"][raw.name] = sha256(raw.read_bytes()).hexdigest()
    path = HERE / "equilibrium_feedback.json"
    if args.write:
        path.write_text(json.dumps(result, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), result)
    print(
        {
            k: [
                v["status"],
                v.get("minimum_drift"),
                v.get("true_inradius_lower"),
                v.get("nominal_residual_HS"),
            ]
            for k, v in result["outcomes"].items()
        }
    )
