"""Directional transport between two independently enrolled balls in a common relay ruler."""

import json
from scipy import sparse
from scipy.optimize import linprog
from scipy.sparse.csgraph import shortest_path
import numpy as np

from ..common_front import measure as front
from ..null_optics import measure as optics


def uniform_transport(distances, i, j, radius):
    left = np.flatnonzero(distances[i] <= radius)
    right = np.flatnonzero(distances[j] <= radius)
    a, b = len(left), len(right)
    constraints = sparse.vstack(
        [
            sparse.kron(sparse.eye(a), np.ones((1, b))),
            sparse.kron(np.ones((1, a)), sparse.eye(b)),
        ],
        format="csr",
    )
    target = np.r_[np.full(a, 1 / a), np.full(b, 1 / b)]
    costs = distances[np.ix_(left, right)].ravel()
    fit = linprog(
        costs, A_eq=constraints, b_eq=target, bounds=(0, None), method="highs"
    )
    if not fit.success:
        raise ArithmeticError(fit.message)
    # Check a feasible primal plan and a dual lower bound.
    potentials = fit.eqlin.marginals
    violation = max(0.0, np.max(constraints.T @ potentials - costs))
    lower = target @ potentials - violation
    gap = fit.fun - lower
    assert max(abs(constraints @ fit.x - target)) < 1e-9 and -1e-9 < gap < 1e-8
    return dict(
        centers=[i, j],
        radius=radius,
        ball_counts=[a, b],
        separation=float(distances[i, j]),
        transport_cost=float(fit.fun),
        primal_dual_gap=float(gap),
        curvature=float(1 - fit.fun / distances[i, j]),
    )


def run():
    board = optics.pairing.formation.federation(1)
    route = optics.pairing.formation.region(board)
    ruler = front.Ruler(board, route, 3)
    state = optics.source()["state"]
    # The entire phase plane in the volume certificate is imaginary Hermitian.
    # A finite perturbation in either basis direction must leave every link unchanged.
    y = np.array([[0, -1j], [1j, 0]])
    z = np.diag([1, -1])
    dx = np.zeros_like(state)
    dx[1] = np.kron(y, z) / 2
    assert np.array_equal(ruler.coefficients(dx), np.zeros(len(ruler.menu)))
    weights = ruler.coefficients(state)
    assert np.array_equal(ruler.coefficients(state + 1e-5 * dx), weights)
    _, _, matrix = ruler.distances(state)
    distances = shortest_path(matrix, directed=False)
    return dict(
        order=3,
        graph_ports=board.n,
        direct_distance=float(distances[0, 39]),
        uniform_vertex_ball_reads=[
            uniform_transport(distances, 0, 39, r) for r in (0.75, 1.5)
        ],
    )


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, allow_nan=False))
