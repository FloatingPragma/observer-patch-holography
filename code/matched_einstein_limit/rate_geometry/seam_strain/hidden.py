"""Certify a collar variation invisible to every read confined to the strain qubit."""

import argparse
import json

import numpy as np
import sympy as sp
from scipy.linalg import qr

from .source import certificate as cert
from .source import measure

HERE = measure.HERE.parent


def correlation_axes():
    x = sp.Matrix([[0, 1], [1, 0]])
    y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    z = sp.diag(1, -1)
    return [sp.kronecker_product(a, b) / 4 for a in (x, y, z) for b in (x, y, z)]


def certificate():
    # A single ordinary/collar reference and a single appointment.
    board = measure.formation.federation(1)
    route = measure.formation.region(board)
    (formed, _), state_error, _ = cert.formed_state(
        board, route, -1, cert.prior.weights()[7]
    )
    a, _, _ = measure.system(board, route, -1)
    for time, wave, column_error in cert.propagate(a, measure.launch(board, route)):
        if time == 4:
            break
    else:
        raise ValueError("missing four-tick pump")
    rail_ops = cert.effects(wave)
    dirs = correlation_axes()
    # Operators P⊗Q/4 have trace norm one and both partial traces vanish.
    effects = [(2 * s + t) / 3 for s, t in zip(*rail_ops, strict=True)]
    response = sp.Matrix([[sp.re(sp.trace(op * d)) for d in dirs] for op in effects])
    lower = sp.Rational(1, 10**8)
    log, log_error = cert.prior.entropy_coordinates(formed, state_error, lower, dirs)
    augmented = response.col_join(log.T)
    # The QR calculation proposes a minor. Only rational arithmetic accepts it.
    _, _, pivot = qr(np.array(augmented, dtype=float), pivoting=True)
    cols = list(map(int, pivot[:4]))
    minor = augmented[:, cols]
    midpoint = minor.applyfunc(lambda x: sp.Rational(round(x * 10**12), 10**12))
    op_error = 2 * column_error  # spectral norm of four columns per classical rail
    entry_errors = [2 * op_error + op_error**2] * 3 + [log_error]
    row_errors = sp.Matrix(
        [
            sum(abs(x) for x in (minor - midpoint).row(i)) + 4 * err
            for i, err in enumerate(entry_errors)
        ]
    )
    inverse_error = max(midpoint.inv().applyfunc(abs) * row_errors)
    assert inverse_error < sp.Rational(1, 10)
    # A concrete prospective direction is for inspection. Invertibility proves
    # existence of the exact kernel witness independently of this approximation.
    witness = np.zeros(9)
    witness[cols] = np.array(midpoint.inv()[:, 3], dtype=float).ravel()
    witness /= np.linalg.norm(witness)
    names = [i + j for i in "XYZ" for j in "XYZ"]
    return dict(
        cells=80,
        forming_charge=-1,
        pump_charge=-1,
        pump_time=4,
        basis="P x Q / 4",
        selected_columns=[names[i] for i in cols],
        rows=["qubit_X", "qubit_Y", "qubit_Z", "regional_modular_nats"],
        minor_midpoint=np.array(midpoint, dtype=float).tolist(),
        per_entry_bounds=[cert.enclosed(e, 0)[1] for e in entry_errors],
        inverse_error_upper=cert.enclosed(inverse_error, 0)[1],
        rank_qubit=3,
        rank_with_modular_read=4,
        example_unit_direction=witness.tolist(),
        example_modular_slope=float(np.array(log, dtype=float).ravel() @ witness),
        formed_state_min_eigenvalue_lower=float(lower),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    data = json.dumps(certificate(), indent=2, allow_nan=False) + "\n"
    if p.parse_args().write:
        (HERE / "hidden.json").write_text(data)
    else:
        print(data, end="")
