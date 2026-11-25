"""Read a modular variation from repeated preparations and two pump launches."""

import argparse
import json

import numpy as np
import sympy as sp
from scipy.linalg import qr

from . import hidden, readout
from .source import certificate as cert
from .source import measure as source


def certified_menu(board, route, axes):
    h, _, _ = source.system(board, route, -1)
    rows, errors, names = [], [], []
    for time, w, col_error in cert.propagate(h, source.launch(board, route)):
        symmetric, antisymmetric = cert.effects(w)
        mixed = [(2 * s + a) / 3 for s, a in zip(symmetric, antisymmetric, strict=True)]
        # A sharp endpoint is (s+a)/sqrt(2). Sum the integer columns and
        # carry the normalization using a rational 1/2 in the effects.
        re, im = np.split(w, 2, axis=1)
        sharp = np.column_stack(
            [
                p
                for part in (re, im)
                for p in (part[:, :4] + part[:, 4:], part[:, :4] - part[:, 4:])
            ]
        )
        left = [m / 2 for m in cert.effects(sharp)[0]]
        for name, ops, error in (
            ("mixed", mixed, 2 * col_error),
            ("left", left, 3 * col_error),
        ):
            for label, op in zip("XYZ", ops, strict=True):
                rows.append([sp.re(sp.trace(op * a)) for a in axes])
                errors.append(2 * error + error**2)
                names.append([time, name, label])
    return sp.Matrix(rows), errors, names


def verify_rank(matrix, errors):
    _, _, indices = qr(np.array(matrix.T, dtype=float), pivoting=True)
    selected = list(map(int, indices[:9]))
    minor = matrix[selected, :]
    mid = minor.applyfunc(lambda v: sp.Rational(round(v * 10**13), 10**13))
    row_error = sp.Matrix(
        [
            sum(abs(v) for v in (mid - minor).row(i)) + 9 * errors[j]
            for i, j in enumerate(selected)
        ]
    )
    bound = max(mid.inv().applyfunc(abs) * row_error)
    assert bound < sp.Rational(1, 10)
    return selected, cert.enclosed(bound, 0)[1]


def independent_populations(board, route, preparations):
    h, _, _ = source.system(board, route, -1)
    out = {charge: [] for charge in preparations}
    for _, wave in source.fast_flow(h, source.launch(board, route) / 2):
        blocks = source.compressed(wave)
        for rail in (np.diag([2 / 3, 1 / 3]), np.ones((2, 2)) / 2):
            for pauli in (readout.X, readout.Y, readout.Z):
                effect = np.einsum("ba,abij->ij", pauli, blocks)
                for charge, states in preparations.items():
                    out[charge].append(
                        [np.trace(effect @ np.kron(rail, rho)).real for rho in states]
                    )
    return out


def run():
    board = source.formation.federation(1)
    route = source.formation.region(board)
    axes = hidden.correlation_axes()
    mat, errors, names = certified_menu(board, route, axes)
    selected, rank_bound = verify_rank(mat, errors)
    (formed, _), err, _ = cert.formed_state(board, route, -1, cert.prior.weights()[7])
    target, log_error = cert.prior.entropy_coordinates(
        formed, err, sp.Rational(1, 10**8), axes
    )
    # Least-squares proposes a left decoder for the fixed target. The certificate
    # bounds every correlation coefficient of the residual in rational arithmetic.
    proposal = np.linalg.lstsq(
        np.array(mat.T, dtype=float), np.array(target, dtype=float), rcond=None
    )[0].ravel()
    coeff = sp.Matrix([sp.Rational(round(float(x) * 10**10), 10**10) for x in proposal])
    weights = np.array(coeff, dtype=float).ravel()
    residual = max(abs(v) for v in (mat.T * coeff - target))
    # Each Tr[(P x Q) delta] is bounded by two for a difference of densities.
    tolerance = cert.enclosed(
        18
        * (
            residual
            + log_error
            + sum(abs(c) * e for c, e in zip(coeff, errors, strict=True))
        ),
        0,
    )[1]
    assert tolerance < 0.001
    sources = {q: source.form(board, route, q) for q in (-1, 1)}
    reads = independent_populations(board, route, sources)
    k = -source.pairing.log_density(sources[-1][0])
    rng = np.random.default_rng(62281)
    shots = 2**36
    # Two held differences, independent shots per rail, setting and preparation.
    noise = 2 * np.linalg.norm(weights) * np.sqrt(np.log(400) / shots)
    records = []
    for q, means in reads.items():
        counts = rng.binomial(shots, (1 + np.array(means)) / 2)
        signed = 2 * counts / shots - 1
        inferred = float(weights @ (signed[:, 1] - signed[:, 0]))
        target_value = float(np.trace(k @ (sources[q][1] - sources[q][0])).real)
        records.append(
            dict(
                forming_charge=q,
                plus_counts=counts.tolist(),
                all_trials_per_arm=shots,
                decoded_modular_nats=inferred,
                direct_modular_nats=target_value,
                exact_population_prediction=float(
                    weights @ np.diff(means, axis=1).ravel()
                ),
            )
        )
    return dict(
        cells=80,
        pump_charge=-1,
        reference_forming_charge=-1,
        directions=[a + b for a in "XYZ" for b in "XYZ"],
        settings=names,
        spanning_setting_indices=selected,
        rank=9,
        inverse_error_upper=rank_bound,
        modular_decoder_weights=weights.tolist(),
        decoder_l2_norm=float(np.linalg.norm(weights)),
        uniform_decoder_error_upper=tolerance,
        seed=62281,
        simultaneous_count_error_99pct=noise,
        reads=records,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    data = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if p.parse_args().write:
        (hidden.HERE / "history.json").write_text(data)
    else:
        print(data, end="")
