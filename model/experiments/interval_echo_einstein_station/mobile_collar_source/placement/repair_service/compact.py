"""Enumerate the closed sharp Markov classes on four physically adjacent owners."""

from fractions import Fraction
from itertools import combinations, product
import numpy as np
from model.candidate_laws.interval_echo import observations as obs
from scipy import sparse
from . import local


def path(board):
    graph = {i: [] for i in range(board.cells)}
    for a, b in board.external:
        u, v = a // 12, b // 12
        graph[u].append(v)
        graph[v].append(u)
    seq = [0]
    while len(seq) < 4:
        seq.append(min(set(graph[seq[-1]]) - set(seq)))
    seams = []
    for u, v in zip(seq, seq[1:]):
        seams.append(
            next(e for e in board.external if {e[0] // 12, e[1] // 12} == {u, v})
        )
    groups = []
    for x, y, _, _ in obs.schedules(board)[0]:
        live = [
            j for j, (a, b) in enumerate(seams) if (a // 12, b // 12) in set(zip(x, y))
        ]
        groups.append(live)
    return seq, seams, groups


def solve(board, name, admission=None):
    owners, seams, groups = path(board)
    edges = [(owners.index(a // 12), owners.index(b // 12)) for a, b in seams]
    words = []
    for pair in combinations(range(4), 2):
        for mem in product(range(7), repeat=3):
            words.append(tuple(int(i in pair) for i in range(4)) + mem)
    lookup = {w: i for i, w in enumerate(words)}
    transition, write, exchange = local.operators()
    if admission is None:
        m = transition[name] / 4096
    else:
        m = (write / 64) @ (
            (1 - float(admission)) * np.eye(28) + float(admission) * exchange / 64
        )
    parts = []
    for i, (a, b) in enumerate(edges):
        row = []
        col = []
        prob = []
        for j, w in enumerate(words):
            code = 4 * w[4 + i] + 2 * w[a] + w[b]
            for end in np.flatnonzero(m[code]):
                new = list(w)
                new[a], new[b], new[4 + i] = (
                    int(local.LEFT[end]),
                    int(local.RIGHT[end]),
                    end // 4,
                )
                row.append(j)
                col.append(lookup[tuple(new)])
                prob.append(m[code, end])
        parts.append(
            sparse.csr_matrix((prob, (row, col)), shape=(len(words), len(words)))
        )
    T = sparse.eye(len(words), format="csr")
    for layer in groups:
        for i in layer:
            T = T @ parts[i]
    np.testing.assert_allclose(np.asarray(T.sum(axis=1)).ravel(), 1, atol=0)
    initial = lookup[(1, 1, 0, 0, 0, 0, 0)]
    nclass, labels = sparse.csgraph.connected_components(
        T, directed=True, connection="strong"
    )
    selected = labels[initial]
    if name != "fresh":
        np.testing.assert_allclose(np.asarray(T.sum(axis=0)).ravel(), 1, atol=0)
        recurrent = labels == selected
        assert T[recurrent, :][:, ~recurrent].nnz == 0
    else:
        sinks = set(range(nclass)) - set(
            labels[T.nonzero()[0]][labels[T.nonzero()[0]] != labels[T.nonzero()[1]]]
        )
        assert len(sinks) == 1
        recurrent = np.isin(labels, list(sinks))
    assert np.max(T.diagonal()[recurrent]) > 0
    # Closed retained classes are exactly uniform. Fresh replacement is
    # nonunital on the system whose archive is omitted; solve by iteration.
    eq = recurrent.astype(float)
    eq /= eq.sum()
    for _ in range(256):
        eq = T.T @ eq
    assert np.linalg.norm(T.T @ eq - eq, 1) < 1e-12
    p = np.zeros(len(words))
    p[initial] = 1
    for _ in range(256):
        p = T.T @ p
    w = np.array(words)
    stock = w[:, :4]
    mem = w[:, 4:]
    oriented = 1 if edges[1][0] == 1 else -1
    observables = {
        "near_minus_far": (stock[:, :2].sum(axis=1) - stock[:, 2:].sum(axis=1)) / 2,
        "signed_cut": oriented * local.SIGNED[mem[:, 1]],
        "cut_current": oriented * local.CURRENT[mem[:, 1]],
        "unit_heads": local.FLAG[mem].sum(axis=1),
        "absolute_cut_excess": local.FLAG[mem[:, 1]]
        - local.FLAG[mem[:, [0, 2]]].mean(axis=1),
        "load_loss_inner": 2 - stock[:, :2].sum(axis=1),
    }
    threshold = None
    if name == "retained_exchange":
        p0 = np.zeros(len(words))
        p0[initial] = 1
        step = T.T.tocsr()
        for time in range(1, 32769):
            p0 = step @ p0
            if p0 @ observables["near_minus_far"] < 0.25:
                threshold = time
                break
    thresholds = {}
    if name != "fresh":
        thresholds = {
            k: str(Fraction(int((2 * y[recurrent]).sum()), 2 * int(recurrent.sum())))
            for k, y in observables.items()
        }
    # Are all six possible K=2 all-blank patches in the original communicating class?
    common_blank = sum(
        labels[lookup[tuple(int(i in pair) for i in range(4)) + (0, 0, 0)]] == selected
        for pair in combinations(range(4), 2)
    )
    # Evaluate the *same* stationary start at each actual write, including
    # the within-cycle stage on the fresh and therefore nonunital control.
    now = eq.copy()
    birth = death = export = 0.0
    flags = local.FLAG[local.MEMORY]
    for stage in groups:
        for i in stage:
            a, b = edges[i]
            code = 4 * mem[:, i] + 2 * stock[:, a] + stock[:, b]
            newflag = m @ flags
            if name == "fresh":
                birth += now @ newflag[code]
                export += now @ flags[code]
            else:
                birth += now @ ((1 - flags) * newflag)[code]
                death += now @ (flags * (1 - newflag))[code]
            now = parts[i].T @ now
    assert abs(birth - death - export) < 1e-12
    return dict(
        owners=owners,
        physical_seams=seams,
        stages=groups,
        words=len(words),
        communicating_classes=nclass,
        recurrent_size=int(recurrent.sum()),
        blank_patches_in_initial_class=int(common_blank),
        first_cycle_contrast_below_quarter=threshold,
        exact_uniform_stationary=thresholds,
        stationary_cut_unit=float(eq @ local.FLAG[mem[:, 1]]),
        stationary_cut_replica_product=float((eq @ local.SIGNED[mem[:, 1]]) ** 2),
        stationary_birth_return_export=[float(birth), float(death), float(export)],
        stationary={k: float(eq @ y) for k, y in observables.items()},
        cycle256={k: float(p @ y) for k, y in observables.items()},
        tv_to_stationary_cycle256=float(abs(p - eq).sum() / 2),
    )
