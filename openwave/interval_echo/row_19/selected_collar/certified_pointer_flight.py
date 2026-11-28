"""Independently bound native flight error on the inconclusive pointer records."""

import argparse
from hashlib import sha256
import json
from math import factorial

import numpy as np

from . import equilibrium_pointer_price as price

HERE = price.HERE
UNIT = 2**56
DEGREE = 20


def integer_columns(H, home, last):
    # One native half-tick is H/4. Sixteen substeps keep ||abs(H)/64|| <= 1.
    assert H.dtype.kind in "iu" and (H - H.T).count_nonzero() == 0
    row = int(abs(H).sum(axis=1).max())
    assert row <= 64
    dim = H.shape[0]
    real = np.zeros((dim, 4), np.int64)
    real[4 * home : 4 * home + 4] = UNIT * np.eye(4, dtype=np.int64)
    imag = np.zeros_like(real)
    floats = price.service.short_steps(-0.25j * H, real.astype(float) / UNIT, last)
    # For term j the column error obeys e_j <= e_(j-1)/j + eps.
    # Hence e_j <= 2 eps. A previous-step error propagates unitarily.
    eps = np.sqrt(2 * dim) / (2 * UNIT)
    per_column = 2 * DEGREE * eps + 2 * np.e / factorial(DEGREE + 1)
    assert 16 * last * per_column < 1e-8
    by_tick = []
    for tick in range(1, last + 1):
        for _ in range(16):
            re, im = real, imag
            maximum = max(np.abs(re).max().item(), np.abs(im).max().item())
            assert row * (maximum + DEGREE) + 64 * DEGREE < 2**63
            assert 4 * (maximum + DEGREE) < 2**63
            total_r, total_i = re.copy(), im.copy()
            for k in range(1, DEGREE + 1):
                half = 32 * k
                re, im = (H @ im + half) // (64 * k), (-H @ re + half) // (64 * k)
                total_r += re
                total_i += im
            real, imag = total_r, total_i
        fixed = (real.astype(float) + 1j * imag.astype(float)) / UNIT
        diff = float(np.linalg.norm(floats[tick] - fixed)) + 2e-13
        # Four separate column errors cost at most twice their largest norm.
        by_tick.append(
            {
                "central_difference_HS_upper": diff,
                "flight_operator_error": diff + 2 * 16 * tick * per_column + 1e-14,
            }
        )
    return {
        "row_norm": row,
        "dimension": dim,
        "per_substep_column_error": per_column,
        "ticks": by_tick,
    }


def lower_drift(data, history, pad):
    w = np.array(data["signed_Pauli_coefficients"])
    vertices = price.disk_vertices()
    coeff = np.column_stack(
        [1 - vertices[:, 0], vertices[:, 0], vertices[:, 1] / 2, -vertices[:, 1] / 2]
    )
    n, scale = data["starts_per_setting"], data["binary_scale"]
    sampling = data["count_error_per_read"]
    error = sampling + pad * (np.linalg.norm(w) + 1e-12)
    record = data["history_results"][history]
    counts = np.array(record["terminal_plus_counts"])
    drift = coeff @ (2 * scale * ((counts - data["reference_plus_count"]) / n))
    measured = scale * (2 * (np.array(record["process_plus_counts"]) / n) - 1)
    coords = np.concatenate(
        [2 * measured[:, :, :1], measured[:, :, 2::2] - measured[:, :, 1::2]], axis=2
    )
    op = coords - np.column_stack((np.zeros(2), w))
    norm = np.linalg.norm(np.einsum("va,abk->vbk", coeff, op), axis=(1, 2))
    cost = np.abs(coeff).sum(axis=1)
    norm += cost * np.sqrt(32) * 2 * error + 1e-11
    bound = drift - cost * error - sampling - 2e-11 - data["reference_HS_radius"] * norm
    index = int(np.argmin(bound))
    return float(bound[index]), index


def run():
    base = HERE / "equilibrium_pointer_price.json"
    data = json.loads(base.read_bytes())
    undecided = [
        k
        for k, r in data["history_results"].items()
        if r["minimum_stationary_drift"] <= 0
    ]
    pairs = {}
    for k in undecided:
        h, q, t = map(int, k.split("/"))
        pairs[h, q] = max(pairs.get((h, q), 0), 2 * t)
    native = price.lab.six.survey.native
    board = native.iso.pairing.formation.federation(1)
    region = native.iso.pairing.formation.region(board)
    replay = {}
    for (home, sign), last in pairs.items():
        replay[f"{home}/{sign}"] = [
            integer_columns(native.generator(board, region, sign, bool(b)), home, last)
            for b in (0, 1)
        ]
    rows = {}
    for key, record in data["history_results"].items():
        original = record["minimum_stationary_drift"]
        row = {
            "prior_lower": original,
            "lower": original,
            "numerical_allowance": data["numerical_effect_allowance"],
        }
        if key in undecided:
            home, sign, t = map(int, key.split("/"))
            bank = replay[f"{home}/{sign}"]
            err = max(
                sector["ticks"][k - 1]["flight_operator_error"]
                for sector in bank
                for k in (t, 2 * t)
            )
            # C is a product of two contractions, with error (2+err)*err.
            # The physical Stinespring pair is [(free+(gamma-1)C), sqrt(1-|gamma|^2)C].
            # Its operator error is < 8 err, hence W's compression costs < 17 err.
            # A further 2e-11 covers the 4x4 products and sums over the native
            # terminal addresses (at most 960), followed by the Pauli projection.
            assert err < 1e-6 and board.n <= 960
            pad = min(row["numerical_allowance"], 24 * err + 2e-11)
            lower, index = lower_drift(data, key, pad)
            row.update(lower=lower, numerical_allowance=pad, worst_vertex=index)
        rows[key] = row
    gap = min(row["lower"] for row in rows.values())
    span = 2 * data["binary_scale"]
    return {
        "scale": UNIT,
        "degree": DEGREE,
        "substeps_per_half_tick": 16,
        "replayed_columns": replay,
        "per_history": rows,
        "argmin": min(rows, key=lambda k: rows[k]["lower"]),
        "uniform_drift_lower": gap,
        "fixed_point_trace_distance_lower": max(0.0, gap / span - 1e-15),
        "necessary_restoring_fraction": max(0.0, gap / (span + max(0.0, gap)) - 1e-15),
        "inputs": {
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (base, HERE / "CERTIFIED_POINTER_FLIGHT.md")
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    r = run()
    dst = HERE / "certified_pointer_flight.json"
    if args.write:
        dst.write_text(json.dumps(r, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(dst.read_bytes()), r)
    print(r["argmin"], r["uniform_drift_lower"], r["necessary_restoring_fraction"])
