"""Calibrate the drift of one signed Pauli observable on all native pointer settings."""

import argparse
from hashlib import sha256
import json

import numpy as np

from . import equilibrium_service as service

obs, lab, HERE = service.obs, service.lab, service.HERE
AXES = obs.pos.AXES
N = 2**60
SEED = 28112114
PHASES = (0, 1, 1j, -1j)
PAD = 3e-8
HISTORIES = 144


def disk_vertices():
    theta = np.arange(64) * (2 * np.pi / 64)
    # Outward dilation absorbs the rounding in cos and all vertex coordinates.
    return np.column_stack((np.cos(theta), np.sin(theta))) * (
        (1 + 1e-12) / np.cos(np.pi / 64)
    )


def choose_witness(base, vertices):
    x, _ = service.center(base["source_heads"])
    rows = np.array(
        [
            [
                service.center(base["terminal_heads"][key][arm])[0] - x
                for arm in service.ARMS
            ]
            for key in base["settings"]
        ]
    )
    dense, sparse, left, right = rows.transpose(1, 0, 2)
    values = (
        dense[:, None]
        + vertices[None, :, 0, None] * (sparse - dense)[:, None]
        + vertices[None, :, 1, None] * ((left - right) / 2)[:, None]
    )
    return service.separator(values.reshape(-1, 30)).reshape(2, 15)


def binary_tape(rng, expectation, scale):
    # Select a Pauli in the measured parity sector, proportional to its
    # coefficient. Fill unused settings with an independent fair sign.
    law = (1 + np.asarray(expectation) / scale) / 2
    assert np.min(law) >= -1e-13 and np.max(law) <= 1 + 1e-13
    shape = (N // 2**52, *law.shape)
    return rng.binomial(2**52, np.clip(law, 0, 1), size=shape).sum(axis=0)


def pulled_prices(columns, home, t, price):
    click = columns[t] @ columns[t, home]
    miss = columns[2 * t] - click

    def adjoint(k):
        return (k.conj().swapaxes(-1, -2) @ (price @ k)).sum(axis=0)

    hit = adjoint(click)
    return [
        adjoint(miss + gamma * click) + (1 - abs(gamma) ** 2) * hit for gamma in PHASES
    ]


def run():
    previous = HERE / "equilibrium_service.json"
    base = json.loads(previous.read_bytes())
    vertices = disk_vertices()
    w = choose_witness(base, vertices)
    W = np.einsum("bk,kij->bij", w, AXES[1:])
    scale = np.max(abs(w).sum(axis=1)) / 2 + 1e-12
    weight_norm = np.linalg.norm(w) + 1e-12
    rng = np.random.default_rng(SEED)
    # Simultaneous Hoeffding event for fewer than 40,000 binary means.
    count_error = scale * np.sqrt(2 * np.log(2 * 40000 / 0.01) / N) + 2e-12
    # The long-flight simulator's trace-error bound holds on every input.
    error = count_error + PAD * weight_norm
    _, sources, *_ = obs.c.trials()
    x = sources["reference"]
    base_hits = int(binary_tape(rng, obs.c.station.pair(W, x), scale))
    preparations = [np.eye(4) / 4]
    for P in AXES[1:]:
        preparations.extend([np.eye(4) / 4 - P / 2, np.eye(4) / 4 + P / 2])
    preparations = np.asarray(preparations)
    native = lab.six.survey.native
    board = native.iso.pairing.formation.federation(1)
    region = native.iso.pairing.formation.region(board)
    ports, _ = lab.six.survey.return_bins()
    out = {}
    for home in ports:
        for sign in (-1, 1):
            flights = []
            for b in (0, 1):
                h = native.generator(board, region, sign, bool(b)) * (-0.25j)
                start = np.zeros((h.shape[0], 4))
                start[4 * home : 4 * home + 4] = np.eye(4)
                flights.append(service.short_steps(h, start, 24).reshape(25, -1, 4, 4))
            for t in range(1, 13):
                G = np.array(
                    [pulled_prices(flights[b], home, t, W[b]) for b in (0, 1)]
                ).transpose(1, 0, 2, 3)
                signals = np.einsum("abij,bji->a", G, x).real
                counts = binary_tape(rng, signals, scale)
                # Equal-population product-Pauli lotteries directly measure
                # every coordinate of each conditional Heisenberg pullback.
                calibrators = np.einsum("abij,pji->abp", G, preparations).real
                cal_hits = binary_tape(rng, calibrators, scale)
                means = scale * (2 * (cal_hits / N) - 1)
                coords = np.concatenate(
                    [2 * means[:, :, :1], means[:, :, 2::2] - means[:, :, 1::2]], axis=2
                )
                incoming = np.column_stack((np.zeros(2), w))
                op = coords - incoming
                # All affine coefficients at a polygon vertex refer to the
                # four physical preparations. The input reference cancels once.
                coeff = np.column_stack(
                    (
                        1 - vertices[:, 0],
                        vertices[:, 0],
                        vertices[:, 1] / 2,
                        -vertices[:, 1] / 2,
                    )
                )
                drift = coeff @ (2 * scale * ((counts - base_hits) / N))
                drift_error = abs(coeff).sum(axis=1) * error + count_error + 1e-11
                opnorm = np.linalg.norm(
                    np.einsum("va,abk->vbk", coeff, op), axis=(1, 2)
                )
                opnorm += abs(coeff).sum(axis=1) * np.sqrt(32) * 2 * error + 1e-11
                lower = (
                    drift
                    - drift_error
                    - base["geometry_reference_HS_radius"] * opnorm
                    - 1e-11
                )
                k = int(np.argmin(lower))
                out[f"{home}/{sign}/{t}"] = {
                    "terminal_plus_counts": counts.tolist(),
                    "process_plus_counts": cal_hits.tolist(),
                    "worst_vertex": k,
                    "minimum_stationary_drift": float(lower[k]),
                    "prepared_drift_at_worst": float(drift[k]),
                    "count_and_flight_allowance_at_worst": float(drift_error[k]),
                    "drift_operator_norm_at_worst": float(opnorm[k]),
                }
    assert len(out) == HISTORIES
    worst = min(out, key=lambda k: out[k]["minimum_stationary_drift"])
    return {
        "seed": SEED,
        "starts_per_setting": N,
        "signed_Pauli_coefficients": w.tolist(),
        "binary_scale": float(scale),
        "count_error_per_read": float(count_error),
        "numerical_effect_allowance": PAD,
        "reference_plus_count": base_hits,
        "reference_HS_radius": base["geometry_reference_HS_radius"],
        "history_results": out,
        "least_history": worst,
        "minimum_drift_all_settings": out[worst]["minimum_stationary_drift"],
        "inputs": {
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (previous, HERE / "EQUILIBRIUM_POINTER_PRICE.md")
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    r = run()
    path = HERE / "equilibrium_pointer_price.json"
    if args.write:
        path.write_text(json.dumps(r, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), r)
    print(r["least_history"], r["minimum_drift_all_settings"])
