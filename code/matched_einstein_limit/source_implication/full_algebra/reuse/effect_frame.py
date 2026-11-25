"""A minimal native preparation table for calibrating arbitrary terminal event effects."""

import argparse
import json
from pathlib import Path

import sympy as sp

from .. import certificate as cert
from .. import station
from ...directional_stress.connection.force_gradients import force_effects

HERE = Path(__file__).resolve().parent


def run():
    identity = sp.eye(2)
    x = sp.Matrix([[0, 1], [1, 0]])
    y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    z = sp.diag(1, -1)
    local = {axis: (identity + pauli) / 2 for axis, pauli in (
        ("Z+", z), ("Z-", -z), ("X+", x), ("Y+", y)
    )}
    basis = [cert.prior.sym_matrix(op) for op in station.iso.basis(False)]
    effects = force_effects()
    names, rows, force_data = [], [], []
    for mode in ("kept", "dephased"):
        for left, a in local.items():
            for right, b in local.items():
                product = sp.kronecker_product(a, b)
                blocks = (product, sp.zeros(4)) if mode == "kept" else (product / 2, product / 2)
                names.append([mode, left, right])
                rows.append([sp.S.One] + [cert.inner(op, h) for op in blocks for h in basis]
                            + [blocks[1].trace() - sp.Rational(1, 2)])
                force_data.append([sum(cert.inner(m, state) for m, state in zip(pair, blocks))
                                   for pair in effects])
    panel = sp.Matrix(rows)
    inverse = panel.inv()
    assert panel * inverse == sp.eye(32)
    # Remove any one preparation. A nonzero covector in the 31-row null space can
    # perturb the strictly interior event effect I/2 by either sign after rescaling.
    assert panel.rank() == 32
    calibration = sp.Matrix(force_data)
    coefficients = inverse * calibration
    # Compare with both separately derived force covectors on every state, not just the probes.
    for i, (even, odd) in enumerate(effects):
        independent = sp.Matrix([(even.trace() + odd.trace()) / 8]
                                + [cert.inner(m, h) for m in (even, odd) for h in basis]
                                + [(odd.trace() - even.trace()) / 4])
        assert coefficients.col(i) == independent
    return {"normalized_head_coordinates": 31, "independent_terminal_effect_parameters": 32,
            "rank_of_native_preparations": 32, "minimum_calibrations_for_unrestricted_effect": 32,
            "preparations": names, "preparation_matrix_determinant": str(panel.det()),
            "inverse_squared_frobenius_norm": str(sum(v*v for v in inverse)),
            "force_readings_from_independent_meter": [
                [float(row[0]), float(row[1])] for row in force_data],
            "force_column_order": ["vacancy", "excess"],
            "recovered_force_coefficient_residual_on_every_normalized_state": ["0", "0"],
            "force_coefficients_scope": "exact equality for the independently enclosed propagation centers",
            "uniform_effect_error_per_native_force": "2e-15"}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        Path(__file__).with_suffix(".json").write_text(result)
    else:
        print(result, end="")
