"""Test the spatial representation of the interacting cubic two-load excitation."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from .cubic_resonance import coupling, green, punctured_cube


HERE = Path(__file__).resolve().parent


def address_turns(n):
    """Half-turns around the x and y axes on the punctured periodic cube."""
    x, y, z = np.indices((n, n, n)).reshape(3, -1)[:, 1:]
    shape = (n, n, n)
    turn_x = np.ravel_multi_index((x, -y, -z), shape, mode="wrap") - 1
    turn_y = np.ravel_multi_index((-x, y, -z), shape, mode="wrap") - 1
    return turn_x, turn_y


def cube(n, critical):
    kinetic, adjacent, distances = punctured_cube(n)
    turn_x, turn_y = address_turns(n)
    xy, yx = turn_x[turn_y], turn_y[turn_x]
    source = np.full(n**3 - 1, 1 / np.sqrt(n**3 - 1))
    rows = []
    for arm in ("critical", "exchange", "native"):
        kappa = coupling(arm, 1 / n, critical)
        hamiltonian = kinetic - sparse.diags(2 * kappa * adjacent)
        defects = [
            int((hamiltonian[turn][:, turn] - hamiltonian).count_nonzero())
            for turn in (turn_x, turn_y)
        ]
        state = expm_multiply((-1j * n**2 / 8) * hamiltonian, source)
        norm = float(np.vdot(state, state).real)
        if abs(norm - 1) > 1e-9:
            raise ArithmeticError("Two-load packet lost normalization")
        overlap = np.vdot(state[xy], state[yx]) / norm
        rows.append(dict(
            arm=arm, kappa_eff=kappa, turn_symmetry_mismatches=defects,
            order_overlap=[float(overlap.real), float(overlap.imag)],
            probability_separation_below_quarter=float(np.sum(abs(state[distances < 0.25])**2)),
        ))
    return dict(n=n, addresses=len(xy), order_mismatches=int(np.count_nonzero(xy != yx)), arms=rows)


def run():
    critical = 1 / (12 * (green(0.0) - 1 / 6))
    spin_x = -1j * np.array([[0, 1], [1, 0]])
    spin_y = -1j * np.array([[0, -1j], [1j, 0]])
    commutator = (spin_x @ spin_y).conj().T @ (spin_y @ spin_x)
    return dict(
        protocol_sha256=sha256((HERE / "electron_rotation_protocol.json").read_bytes()).hexdigest(),
        spinor_control_order_comparison_real=commutator.real.tolist(),
        spinor_control_order_comparison_imag=commutator.imag.tolist(),
        cubes=[cube(n, critical) for n in (8, 16, 32)],
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.parse_args().output.write_text(json.dumps(run(), indent=2, allow_nan=False) + "\n")
