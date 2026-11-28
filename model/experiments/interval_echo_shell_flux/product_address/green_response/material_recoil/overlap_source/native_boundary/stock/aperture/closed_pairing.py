"""Closed-inventory overlap readouts, local-stress matches and finite feasibility certificates."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from . import record_density as source

HERE = source.HERE
stock = source.stock
MENUS = ("energy", "single_records", "seams", "spectators")
SHOTS = 2**20


def setup():
    full = source.classes()[1][2]
    mask = ((full[:, 4:6] != 0) & (stock.MAG[full[:, 4:6]] != 1)).all(axis=1)
    w = full[mask]
    _, _, _, edges, *_ = stock.carrier()
    _, _, H, _ = stock.kernels()
    codes = [9 * w[:, 4 + i] + 3 * w[:, a] + w[:, b] for i, (a, b) in enumerate(edges)]
    bits = (stock.MAG[w[:, 4:]] == 2).astype(int)
    cells = 2 * bits[:, 0] + bits[:, 1]
    fuel = (w[:, :4] == 2).sum(axis=1) + bits[:, 2]
    np.testing.assert_array_equal(fuel, 2 - bits[:, :2].sum(axis=1))
    np.testing.assert_array_equal(bits[:, 0] * bits[:, 1], fuel == 0)
    # Equalize sharp weights on each connected off-diagonal component of
    # all three simultaneous free Hamiltonians. Check closure on every column.
    lookup = {tuple(x): i for i, x in enumerate(w)}
    adj = np.eye(len(w), dtype=bool)
    for i, (a, b) in enumerate(edges):
        for idx, word in enumerate(w):
            for end in np.flatnonzero(H[:, codes[i][idx]]):
                dst = word.copy()
                dst[4 + i], material = divmod(end, 9)
                dst[a], dst[b] = divmod(material, 3)
                adj[idx, lookup[tuple(dst)]] = True
    assert np.array_equal(adj, adj.T)
    groups, ids = connected_components(csr_matrix(adj), directed=False)
    return full, w, cells, fuel, np.array(codes), H, ids, groups


def flags(values):
    return [values == x for x in np.unique(values)]


def constraints(w, cells, codes, h, menu):
    x, y = cells // 2, cells % 2
    # The last head cell remains free; its coefficient is the objective.
    rows = [np.ones(len(w), int), x, y, *h.diagonal()[codes]]
    if menu != "energy":
        for v in w.T:
            rows.extend(flags(v))
    if menu in ("seams", "spectators"):
        for v in codes:
            rows.extend(flags(v))
    if menu == "spectators":
        _, idx = np.unique(w[:, [0, 1, 2, 3, 6]], axis=0, return_inverse=True)
        rows.extend(flags(idx))
    return np.asarray(rows, dtype=np.int64)


def rational(values):
    return [Q(float(v)).limit_denominator(10**7) for v in values]


def group_matrix(rows, ids, groups):
    return np.array(
        [np.bincount(ids, weights=v, minlength=groups) for v in rows], dtype=np.int64
    )


def certify(rows, ids, groups, double, goal):
    A = group_matrix(rows, ids, groups)
    b = A.sum(axis=1)
    c = group_matrix(double[None, :], ids, groups)[0]
    opt = linprog(-c, A_eq=A, b_eq=b, bounds=(0, None), method="highs")
    assert opt.success
    # The nonnegative primal and upper supporting plane are both checked
    # exactly. The optimizer only proposes a witness.
    z = np.array(rational(opt.x), dtype=object)
    dual = np.array(rational(-opt.eqlin.marginals), dtype=object)
    assert all(z >= 0) and np.all(A @ z == b)
    assert np.all(dual @ A >= c) and dual @ b == c @ z
    native = int(double.sum())
    highest = c @ z
    witness = {
        "largest_double_count": str(highest),
        "fixed_native_count": native,
        "dual_nonzero": [[i, str(v)] for i, v in enumerate(dual) if v],
    }
    if highest >= goal:
        # Find the least severe rejection filter at the physical target.
        target_rows = np.vstack([A, c])
        problem = linprog(
            np.r_[np.zeros(groups), 1],
            A_eq=np.c_[target_rows, np.zeros(len(target_rows))],
            b_eq=np.r_[b, goal],
            A_ub=np.c_[np.eye(groups), -np.ones(groups)],
            b_ub=np.zeros(groups),
            bounds=(0, None),
            method="highs",
        )
        assert problem.success
        z = np.array(rational(problem.x[:-1]), dtype=object)
        assert all(z >= 0) and np.all(A @ z == b) and c @ z == goal
        # Same maximum weight for every member of each free component.
        witness.update(weights=[str(v) for v in z], largest_weight=str(max(z)))
    return witness


def closed_counts(w, weights, rng):
    draw = rng.multinomial(SHOTS, weights / weights.sum())
    big = stock.MAG[w[:, 4:]] == 2
    fuel = (w[:, :4] == 2).sum(1) + big[:, 2]
    table = np.bincount(2 * big[:, 0] + big[:, 1], weights=draw, minlength=4).astype(
        int
    )
    observed = np.bincount(fuel, weights=draw, minlength=3).astype(int)
    # Audit preparation IDs directly via the truth table on every accepted word.
    assert np.array_equal(observed, [table[3], table[1] + table[2], table[0]])
    middle_nonunit = int(draw @ (w[:, 1] != 1))
    assert middle_nonunit == table[0]
    return dict(
        pair_counts=table.tolist(),
        spectator_counts_fuel_0_1_2=observed.tolist(),
        middle_nonunit=middle_nonunit,
    )


def stationary_target():
    full, w, cells, fuel, codes, h, ids, groups = setup()
    double = (cells == 3).astype(int)
    n = len(w)
    t = np.bincount(cells, minlength=4).reshape(2, 2)
    goal = Q(int(t[1].sum() * t[:, 1].sum()), n)
    assert goal.denominator == 1
    panels = {}
    for m in MENUS:
        rows = constraints(w, cells, codes, h, m)
        panels[m] = certify(rows, ids, groups, double, goal.numerator)
    return full, w, cells, fuel, ids, panels


def run():
    full, w, cells, fuel, ids, panels = stationary_target()
    rng = np.random.default_rng(5947032)
    n = len(w)
    controls = {}
    table = np.bincount(cells, minlength=4).reshape(2, 2)
    product = np.outer(table.sum(1), table.sum(0))
    simple = [Q(int(product.ravel()[k]), n * int(table.ravel()[k])) for k in cells]

    # At fixed one-head marginals the spectator count bins change by D,-2D,D.
    def tv(weights, labels):
        return sum(
            abs(sum(weights[j] - 1 for j in np.flatnonzero(labels == a)))
            for a in np.unique(labels)
        ) / (2 * n)

    _, spect_id = np.unique(w[:, [0, 1, 2, 3, 6]], axis=0, return_inverse=True)
    source_shift = np.array([1, -1, -1, 1]) / 14

    def deposit(name, weights):
        weights = np.array(weights, dtype=object)
        actual = np.array([sum(weights[cells == j]) / n for j in range(4)], float)
        original = table.ravel() / n
        result = closed_counts(w, np.asarray(weights, float), rng)
        result["true_cell_shift"] = (actual - original).tolist()
        result["spectator_total_variation"] = str(tv(weights, spect_id))
        result["spectator_fuel_total_variation"] = str(tv(weights, fuel))
        controls[name] = result

    deposit("native", [Q(1)] * n)
    deposit("cell_reweighted", simple)
    assert controls["cell_reweighted"]["spectator_total_variation"] == "1/7"
    for m, v in panels.items():
        if "weights" in v:
            deposit(m, [Q(v["weights"][i]) for i in ids])
            np.testing.assert_allclose(
                controls[m]["true_cell_shift"], source_shift, atol=1e-16
            )
    # Independent second-head reports can be created after the records were
    # made. Keep the actual first copy's spectator readings as a blind check.
    i, j = rng.integers(n, size=(2, SHOTS))
    observed = (cells[i] // 2) * (cells[j] % 2)
    predicted = fuel[i] == 0
    controls["cross_copy"] = {
        "native_spectator_and_edited_pair": np.bincount(
            2 * predicted + observed, minlength=4
        )
        .reshape(2, 2)
        .tolist()
    }
    # The initial number of doubles may instead be an observed preparation tag.
    test = {}
    for k, words in enumerate(source.classes()[1]):
        mag = stock.MAG[words[:, 4:]]
        pred = k - (words[:, :4] == 2).sum(1) - (mag[:, 2] == 2)
        assert np.array_equal(pred, (mag[:, :2] == 2).sum(1))
        assert np.all(words[:, :4].sum(1) == 2 * k)
        assert np.array_equal(
            pred == 2, (words[:, :4] == 1).all(axis=1) & (mag[:, 2] != 2)
        )
        accept = ((words[:, 4:6] != 0) & (mag[:, :2] != 1)).all(axis=1)
        both_zero = (mag[:, :2] == 0).all(axis=1)
        middle = words[:, 1] != 1
        test[str(k)] = dict(
            middle_nonunit_vs_both_zero=np.bincount(
                (2 * middle + both_zero)[accept], minlength=4
            )
            .reshape(2, 2)
            .tolist(),
            words=len(words),
            accepted=int(sum(accept)),
            double_two=int(sum((pred == 2) & accept)),
        )
    # Removing the native echo leaves all three original blank memories blank.
    # Check the exchange control's complete transition columns for that subspace.
    _, exchange, _, _ = stock.kernels()
    assert not exchange[9:, :9].any()
    return dict(
        exchange_only_accepted_from_blank=0,
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                HERE / "CLOSED_PAIRING.md",
                HERE / "pair_contact.json",
                HERE.parent / "parked.json",
            )
        },
        total_words=len(full),
        accepted=n,
        free_components=len(set(ids)),
        shots_per_accepted_bank=SHOTS,
        menus=panels,
        readouts=controls,
        inventory_check=test,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    opt = cli.parse_args()
    result = run()
    path = HERE / "closed_pairing.json"
    if opt.write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        assert json.loads(path.read_bytes()) == result
    print(
        {
            k: (v["largest_double_count"], v.get("largest_weight"))
            for k, v in result["menus"].items()
        }
    )
    print(
        {k: v.get("spectator_total_variation") for k, v in result["readouts"].items()}
    )
