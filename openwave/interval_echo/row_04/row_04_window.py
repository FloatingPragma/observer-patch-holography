"""Lower-bound the home-pair Born probability throughout the free interval."""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy import sparse

from row_04_check import screen


def continuous_window(board):
    edges, indices, exchange, adjacent = screen(board)
    internal = edges[(edges[:, 0] // 12 == 13) & (edges[:, 1] // 12 == 13)]
    gateways = np.unique(edges[edges[:, 0] // 12 != edges[:, 1] // 12])
    source = internal[np.isin(internal, gateways).any(axis=1)][0]
    patch = indices[internal[:, 0], internal[:, 1]]
    source_index = indices[tuple(source)]

    # Adding 2 times the identity only changes the overall phase.
    h = (exchange + sparse.diags(2 - 16 * adjacent)).tocsr()
    max_abs_row_sum = float(np.asarray(abs(h).sum(axis=1)).max())
    max_terms = int(np.diff(h.indptr).max())
    if max_abs_row_sum > 14 or max_terms > 13:
        raise ArithmeticError("The operator exceeds the numerical enclosure")
    column = h.getcol(source_index)
    variance = float(column.multiply(column).sum() - h[source_index, source_index] ** 2)
    dt, degree, steps = 1 / 8, 24, 296
    scaled_norm = max_abs_row_sum * dt
    taylor_tail = math.exp(scaled_norm) * scaled_norm ** (degree + 1) / math.factorial(degree + 1)
    # Budget for complex arithmetic at double precision in each Taylor step.
    u = np.finfo(float).eps / 2
    gamma = 128 * u / (1 - 128 * u)
    local_error = taylor_tail + math.exp(scaled_norm) * math.expm1((degree + 2) * math.log1p(gamma))
    vector_error = math.expm1(steps * math.log1p(local_error))
    probability_error = 1e-6
    if 2 * vector_error + vector_error**2 + gamma >= probability_error:
        raise ArithmeticError("Floating-point uncertainty exceeds the declared budget")

    state = np.zeros(h.shape[0], complex)
    state[source_index] = 1
    weights = [1.0]
    for _ in range(steps):
        term = state
        total = state.copy()
        for j in range(1, degree + 1):
            term = ((-1j * dt) * (h @ term)) / j
            total += term
        state = total
        local = state[patch]
        weights.append(float(np.vdot(local, local).real))

    grid_floor = min(weights) - probability_error
    # Distance in projective Hilbert space to this patch's material subspace
    # is acos(sqrt(p)). Delta-H bounds its change per unit exposure.
    angle = math.acos(math.sqrt(grid_floor)) + math.sqrt(variance) * dt / 2
    if angle >= math.pi / 2:
        raise ArithmeticError("The continuous bound is vacuous")
    guaranteed_percent = math.floor(100 * math.cos(angle) ** 2)
    return dict(source=source.tolist(), kappa=8, last_exposure=37, steps=steps,
                step=dt, variance=variance, polynomial_degree=degree, max_terms=max_terms,
                absolute_row_sum=max_abs_row_sum, grid_probability_error_budget=probability_error,
                grid_minimum=round(min(weights), 9), grid_minimum_exposure=dt * weights.index(min(weights)),
                sample_probabilities=[round(weights[round(t / dt)], 9) for t in (1, 7, 37)],
                all_exposures_lower_bound_percent=guaranteed_percent)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--board", type=Path, required=True)
    ap.add_argument("--record", action="store_true", help="create the supplemental receipt once")
    args = ap.parse_args()
    folder = Path(__file__).resolve().parent
    report = continuous_window(args.board)
    original = json.loads((folder / "row_04_receipt.json").read_text())
    for weight, reading in zip(report["sample_probabilities"], original["arms"][0]["readings"], strict=True):
        if abs(weight - sum(reading["probabilities"][:-1])) > report["grid_probability_error_budget"]:
            raise ArithmeticError("Independent grid evolution disagrees with the terminal read")
    receipt = folder / "row_04_window.json"
    if args.record:
        with receipt.open("x") as f:
            f.write(json.dumps(report, indent=2) + "\n")
    if report != json.loads(receipt.read_text()):
        raise SystemExit("The calculated bound differs from the supplemental receipt")
    print(json.dumps(report, indent=2))
