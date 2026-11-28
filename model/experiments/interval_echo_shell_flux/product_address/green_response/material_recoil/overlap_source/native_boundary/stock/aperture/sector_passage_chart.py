"""Protected blank-link service and calibrated passage times on the occupied board."""

import numpy as np
from scipy.linalg import solve
from scipy.sparse.csgraph import shortest_path

from . import sector_ball_source as source
from . import sector_boundary_clock as boundary


def apparatus():
    words = source.thermal.apparatus()[0]
    board, homes, seams, *_ = source.balls.energy.returned(words)
    graph, _, _ = source.balls.radar.chart(board, seams)
    rail = 38
    ports = np.delete(np.arange(board.n), rail)
    graph = graph[ports][:, ports].tocsr()
    graph.eliminate_zeros()
    adj = graph.toarray().astype(int)
    loc = np.searchsorted(ports, homes)
    # Physical neighbors establish clearance for independent launched material;
    # ordinary traffic uses every surviving bypass link on the survey chart.
    distance = shortest_path(board.adjacency(), unweighted=True)
    clearance = distance[
        np.ix_(ports, np.unique([rail, *homes, *np.ravel(seams)]))
    ].min(1)
    assert np.isfinite(shortest_path(graph, unweighted=True)).all()
    uv = np.transpose(np.nonzero(np.triu(adj, 1)))
    inventories, group = np.unique(words[:, :4], axis=0, return_inverse=True)
    # Include ordinary interventions from the rail even when an arm lies
    # outside the supported peer preparation. The observer's link clock is
    # defined by native gates on the larger physical carrier in both arms.
    cal = []
    for home in inventories[inventories[:, 1] == 0]:
        for t in (0, 1, 2):
            after = home.copy()
            after[1] = t
            cal.append(after)
    all_loads = np.unique(np.vstack((inventories, cal)), axis=0)
    times = boundary.rates(graph, all_loads, loc, *uv.T)
    periods = {}
    for home, ticks in zip(all_loads, times, strict=True):
        occ = np.zeros(len(ports), int)
        occ[loc] = home > 0
        assert np.array_equal(ticks, 1 + occ[uv[:, 0]] + occ[uv[:, 1]])
        assert max(ticks) <= 2
        periods[tuple(home > 0)] = ticks
    return words, board, ports, clearance, adj, uv, inventories[group], periods


def inverse(A):
    # Certificate for every entry of the grounded integer Laplacian inverse.
    # If B is the dyadic trial inverse and E=I-A*B, then
    # ||A^-1||_inf <= ||B||_inf/(1-||E||_inf).
    grid = 2**46
    n = len(A)
    approx = solve(A, np.eye(n), assume_a="pos")
    B = np.rint(approx * grid).astype(np.int64)
    assert np.abs(B).max() * np.abs(A).sum(1).max() < 2**61
    E = np.eye(n, dtype=np.int64) * grid - A @ B
    err = np.abs(E).sum(1).max() / grid
    assert err < 1
    cap = np.abs(B).sum(1).max() / grid / (1 - err)
    return B / grid, float(cap * (np.abs(E).max() / grid) + 1e-13)


def law(adj, uv, ticks):
    n = len(adj)
    # Two equiprobable phases for every chosen slot. An idle cycle retains
    # the current address. Each accepted edge also passes both material loads.
    step = np.zeros_like(adj)
    u, v = uv.T
    step[u, v] = step[v, u] = 2 // ticks
    step[np.diag_indices(n)] = 12 - step.sum(1)
    A = 12 * np.eye(n, dtype=int) - step
    B, eps = inverse(A[:-1, :-1])
    G = np.zeros((n, n))
    G[:-1, :-1] = 2 * B
    diag = np.diag(G)
    R = diag[:, None] + diag - 2 * G
    costs = np.full(n, 6)
    for ends in (u, v):
        np.add.at(costs, ends, ticks - 1)
    return step, G, R, 8 * eps + 1e-13, costs


def bearing(R, ab, err):
    a, b = ab
    span = R[a, b]
    prob = (1 + (R[:, a] - R[:, b]) / span) / 2
    assert min(prob) >= -1e-11 and max(prob) <= 1 + 1e-11 and span > err
    return prob, 1.5 * err / (span - err) + 1e-13


def serial(adj, reference, difference):
    # delta H - alpha*H0 is a pure difference of station offsets.
    # G is grounded; zero total on the right makes the field independent
    # of the ground up to an additive constant.
    L = np.diag(adj.sum(1)) - adj
    g, _ = inverse(L[:-1, :-1])
    G = np.zeros(adj.shape)
    G[:-1, :-1] = g
    jterm = np.diag(G) - G

    def times(cost):
        potential = G @ cost
        return potential[:, None] - potential + cost.sum() * jterm

    t0 = times(reference)
    dt = times(difference)
    rate = difference.sum() / reference.sum()
    field = G @ (difference - rate * reference)
    assert np.max(abs(dt - rate * t0 - (field[:, None] - field))) < 1e-7
    assert np.max(abs(dt + dt.T - rate * (t0 + t0.T))) < 1e-7
    return {
        "relative_commute_change": float(rate),
        "station_offsets_empty_link_times": field.tolist(),
        "maximum_offset_difference": float(np.ptp(field)),
    }
