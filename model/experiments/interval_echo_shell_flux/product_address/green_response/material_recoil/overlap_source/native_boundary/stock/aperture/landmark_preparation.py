"""A modular repair at fixed ordinary energies and fixed selective-contact records."""

from fractions import Fraction as Q

import numpy as np
from sympy import Matrix, diag

from . import collective_source as source
from . import fluctuation_return
from . import landmark_entropy
from . import collective_retention as retention


def persistence(words, prob):
    full, original, _, service, _ = retention.setup()
    ids = {tuple(w): i for i, w in enumerate(full)}
    repaired = np.full(len(full), Q(0), object)
    repaired[[ids[tuple(w)] for w in words]] = prob
    for _, exchange in service:
        row, col = exchange.nonzero()
        assert np.array_equal(repaired[row], repaired[col])
    grid = retention.GRID
    delta = np.array([round(x * grid) for x in repaired - original], np.int64)
    assert int(max(abs(delta))) * 128 + 64 < 2**63
    blank = full[:, -1] == 0
    result = {}
    for i in range(4097):
        if i in (0, 1, 4, 16, 64, 256, 1024, 4096):
            error = Q((1 + 3 * i) * len(full), 2 * grid)
            tv = (Q(sum(abs(int(v)) for v in delta), grid) + error) / 2
            val = Q(sum(map(int, delta[blank])), grid)
            result[str(i)] = dict(
                tv_upper=str(tv),
                third_blank_change=[float(val - error), float(val + error)],
            )
        if i < 4096:
            for mat, _ in service:
                delta = (mat @ delta + 64) // 128
    return result


def prepare():
    energy = source.tidal.lab.energy
    _, w, cells, _, codes, h, ids, n = energy.pair.setup()
    _, probs, _ = source.records()
    original = probs["original"]
    target = probs["source"]
    partial = np.isin(w[:, -1], source.PARTIAL)
    blank = w[:, -1] == 0
    # The registered ordinary-unit comparison at 13,20 writes the signed
    # middle reading to the fresh head. Every native write ticket gives it.
    _, write, contrast, _ = fluctuation_return.native()
    signed = -contrast[write[3 + w[:, 1]] // 9]
    np.testing.assert_array_equal(
        signed, np.broadcast_to((w[:, 1] - 1)[:, None], signed.shape)
    )
    free = ~np.isin(np.arange(n), ids[partial])
    assert np.all(free[ids] == ~partial)
    sizes = np.bincount(ids)
    rows = np.vstack(
        (
            [cells == i for i in range(4)],
            w[:, :4].T,
            h.diagonal()[codes],
            energy.returned(w)[-1],
        )
    )
    target_read = rows @ original
    target_read[:4] = (rows @ target)[:4]
    gram_weight = diag(*[Q(1, int(g)) for g in sizes[free]])

    def solve(rows, b):
        # For squared change per *word*, a group of n_g equiprobable
        # records costs n_g times its squared per-word change.
        C = Matrix(energy.pair.group_matrix(rows, ids, n)[:, free])
        rhs = Matrix((b - rows @ original).tolist())
        independent = list(C.T.rref()[1])
        T = C[independent, :]
        lagrange = (T * gram_weight * T.T).inv() * rhs[independent, :]
        diff = gram_weight * T.T * lagrange
        assert C * diff == rhs and gram_weight.inv() * diff == T.T * lagrange
        groups = np.full(n, Q(1, len(w)), object)
        groups[free] += np.array([Q(int(v.p), int(v.q)) for v in diff], object)
        z = groups[ids]
        assert min(z) > 0 and sum(z) == 1
        assert np.all(z[partial] == original[partial]) and np.all(rows @ z == b)
        return z, groups

    base_a = sum(original[blank])
    records = {}
    for name, desired in (
        ("minimum_change", None),
        ("zero", Q(0)),
        ("reversed", -Q(1, 112)),
    ):
        if desired is None:
            z, groups = solve(rows, target_read)
        else:
            z, groups = solve(
                np.vstack((rows, blank)), np.r_[target_read, base_a + desired]
            )
        middle = w[:, 1] - 1
        assert z @ middle == 0 and z @ (middle**2) == z @ (cells == 0)
        assert z @ (middle**2) - original @ (middle**2) == Q(1, 14)
        rate = sum(z[blank]) - base_a
        heads = [sum(z[w[:, -1] == j]) for j in range(11)]
        # Rejection from an equiprobable returned word uses the maximum
        # group weight as a common denominator and retains both outcomes.
        accept = 1 / (len(w) * max(groups))
        records[name] = dict(
            blank_change=str(rate),
            min_word_weight=str(min(groups)),
            uniform_proposal_acceptance=str(accept),
            group_weights=list(map(str, groups)),
            head_fractions=list(map(str, heads)),
            ordinary_reads=list(map(str, rows[4:] @ z)),
            modular_cells=list(map(str, rows[:4] @ z)),
            middle_write_variance=str(z @ (middle**2)),
        )
        if desired is not None:
            assert rate == desired
    ideal, primary = landmark_entropy.certify(True)
    unrestricted, _ = landmark_entropy.certify(False)
    records["entropy_stationary"] = dict(
        blank_change=str((primary - original) @ blank),
        min_word_weight=str(min(primary)),
        head_fractions=[str(sum(primary[w[:, -1] == j])) for j in range(11)],
        ordinary_reads=list(map(str, rows[4:] @ primary)),
        modular_cells=list(map(str, rows[:4] @ primary)),
        middle_write_variance=str(primary @ (middle**2)),
        **ideal,
    )
    # Acquire the actual word proposal and acceptance, keeping the count
    # of declined preparations. The original arm returns every proposal.
    rng = np.random.default_rng(572114)
    live = {}
    for name, z in (("original", original), ("entropy_stationary", primary)):
        proposals = 2**22
        candidates = rng.multinomial(proposals, np.full(len(w), 1 / len(w)))
        accepted = rng.binomial(candidates, np.array(z / max(z), float))
        N = int(sum(accepted))
        live[name] = dict(
            proposals=proposals,
            accepted=N,
            declined=proposals - N,
            blank=int(accepted @ blank),
            middle_squared=int(accepted @ (signed[:, 0] ** 2)),
        )
    stats = {}
    for what in ("blank", "middle_squared"):
        vals = [a[what] / a["accepted"] for a in live.values()]
        err = sum(
            np.sqrt(np.log(8 / 0.003) / (2 * a["accepted"])) for a in live.values()
        )
        stats[what] = [vals[1] - vals[0] - err, vals[1] - vals[0] + err]
    assert (
        0
        < stats["blank"][0]
        <= Q(records["entropy_stationary"]["blank_change"])
        <= stats["blank"][1]
    )
    assert 0 < stats["middle_squared"][0] <= Q(1, 14) <= stats["middle_squared"][1]
    return dict(
        preparations=records,
        unfrozen_selective_entropy_control=unrestricted,
        reopened_write_sweeps=persistence(w, primary),
        source_counts=live,
        measured_changes=stats,
        family_error=0.003,
    ), primary
