"""Pool charge-blind return histories and test the resulting four-station distance matrix."""

from itertools import combinations
import json
from math import comb, factorial
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.linalg import helmert

from .measure import formation, generator, parity_blocks
from .certificate import integer_flight


def pooled_ranges(panels, arm, horizon):
    """Mix complete return laws before applying the capped-clock estimator."""
    step = panels[0]["step"]
    stop = round(horizon / step)
    clock = np.array(panels[0]["clock"])
    elapsed = np.unwrap(np.arctan2(2*clock[:, 1]-1, 1-2*clock[:, 0]))[:stop+1]
    result = np.zeros((4,4))
    for choices in zip(*(p["pairs"] for p in panels),strict=True):
        p,q = choices[0]["p"],choices[0]["q"]
        assert all((c["p"],c["q"]) == (p,q) for c in choices)
        returns = np.mean([c["arms"][arm]["first_return"][:stop+1] for c in choices],axis=0)
        relays = np.mean([c["arms"][arm]["first_relay"][:stop+1] for c in choices],axis=0)
        miss_relay = 1-relays.sum()
        wait_return = relays.sum()-returns.sum()
        assert min(miss_relay,wait_return) > -1e-12
        # Both censored outcomes contribute the same time at the cap.
        result[p,q] = (returns @ elapsed + (miss_relay + wait_return)*elapsed[-1]) / 2
    return result


def distance_check(ranges):
    """Fit endpoint gains by unweighted least squares in log distance."""
    log_range = np.zeros_like(ranges,dtype=float)
    off = ~np.eye(4,dtype=bool)
    log_range[off] = np.log(ranges[off])
    skew = log_range - log_range.T
    cycles = np.array([skew[i,j]+skew[j,k]+skew[k,i] for i,j,k in combinations(range(4),3)])
    log_gain = skew.mean(axis=1)
    corrected = log_range - log_gain[:,None]
    symmetric = np.exp((corrected + corrected.T)/2)
    np.fill_diagonal(symmetric,0.)
    directional_residual = (corrected-corrected.T) / 2
    margins = []
    for i,j,k in combinations(range(4),3):
        sides = (symmetric[i,j],symmetric[j,k],symmetric[k,i])
        margins.append(sum(sides) - 2*max(sides))
    centered = helmert(4)
    eigs = np.linalg.eigvalsh(-.5*centered @ (symmetric**2) @ centered.T)
    return dict(ranges=ranges.tolist(), cycle_ratios=np.exp(cycles).tolist(),
                gains=np.exp(log_gain).tolist(), log_residual_max=float(abs(directional_residual).max()),
                necessary_log_error_floor=float(abs(cycles).max()/6), projected=symmetric.tolist(),
                triangle_margins=margins, centered_eigenvalues=eigs.tolist())


def leading_reciprocity():
    """Full-coin equality of pooled first coefficients along every unique shortest path."""
    board = formation.federation(1)
    route = formation.region(board)
    for parity in (1,-1):
        identity = np.eye(4,dtype=np.int64)
        legs = {}
        for q in (1,-1):
            sequence = []
            for i in range(3):
                if i == 1:
                    sequence.append(-identity)
                    continue
                blocks = 2*parity_blocks(q,8,parity if i == 0 else 1)
                dst = 2 if route[i] < route[i+1] else 1
                entry = blocks[dst,3-dst]
                np.testing.assert_array_equal(entry,entry.astype(np.int64))
                eye = np.eye(2,dtype=np.int64)
                sequence.append(np.kron(entry.astype(np.int64),eye) if i == 0
                                else np.kron(eye,entry.astype(np.int64)))
            legs[q] = sequence
        for i,j in combinations(range(4),2):
            effects = []
            for q in (1,-1):
                path = identity
                for hop in legs[q][i:j]:
                    path = hop @ path
                forward, reverse = path.T @ path, path @ path.T
                effects.append((forward @ forward, reverse @ reverse))
            np.testing.assert_array_equal(effects[0][0],effects[1][1])
            np.testing.assert_array_equal(effects[1][0],effects[0][1])
    return 12


def generator_difference():
    board = formation.federation(1)
    route = formation.region(board)
    coin_flip = sparse.kron(sparse.eye(board.n), np.eye(4)[::-1],format="csr")
    target = 32 * np.repeat(6-board.degree,4)
    for odd in (False,True):
        d = generator(board,route,1,8,odd) - coin_flip @ generator(board,route,-1,8,odd) @ coin_flip
        d = d - sparse.diags(target-target.mean())
        assert not len(d.data) or max(abs(d.data)) < 2e-12
    # In the even-head, exchange-only arm the signs agree even without changing the coins.
    assert (generator(board,route,1,0,False)-generator(board,route,-1,0,False)).nnz == 0
    return int(np.ptp(target))


def correction_at_sixth_order():
    """Charge-pooled forward minus reverse two-leg effect at ports 38,39."""
    board = formation.federation(1)
    route = formation.region(board)
    p, q = route[2:]
    coefficients = np.zeros((7,4,4),dtype=object)
    for charge in (-1,1):
        h = integer_flight(board,route,charge)
        assert int(abs(h).sum(axis=1).max())**6 < np.iinfo(np.int64).max
        v = np.zeros((h.shape[0],4),dtype=np.int64)
        v[4*p:4*p+4] = np.eye(4)
        blocks = []
        for order in range(7):
            blocks.append(v[4*q:4*q+4].astype(object))
            if order != 6:
                v = h @ v
        for orientation, a in ((1,blocks),(-1,[m.T for m in blocks])):
            # Derivatives in -it/2. m describes one detection at each of t and 2t.
            m = [sum(comb(n,j)*a[j].T@a[n-j] for j in range(n+1)) for n in range(7)]
            for n in range(7):
                coefficients[n] += orientation * sum((-1)**j*comb(n,j)*m[j].T@m[n-j]
                                                       for j in range(n+1))
    assert not np.any(coefficients[:6])
    # Recover the coefficient of t^6 with (-i)^6, the factorial, the factor
    # 2 from the doubled Hamiltonian at every order and equal species weights.
    target = np.diag([-1,1,-1,1])
    np.testing.assert_array_equal(-3*coefficients[6],2048*2*2**6*factorial(6)*target)
    return "−(2048/3) Z_right"


def run():
    archived = json.loads(Path(__file__).with_name("results.json").read_text())
    groups = {}
    for run in archived["runs"]:
        groups.setdefault((run["contrast"],run["left_odd"]),{})[run["probe"]] = run
    rows = []
    for (k,odd), group in groups.items():
        choices = [(str(s),[r]) for s,r in group.items()]
        # The missing negative exchange-only control has the identical generator on the same state.
        choices.append(("pooled",[group[1],group.get(-1,group[1])]))
        for label, sources in choices:
            for panels in zip(*(s["panels"] for s in sources), strict=True):
                for arm in ("formed","independent"):
                    for horizon in (2,4,6):
                        ranges = pooled_ranges(panels,arm,horizon)
                        rows.append(dict(contrast=k,left_odd=odd,species=label,arm=arm,
                                         step=panels[0]["step"],horizon=horizon,
                                         **distance_check(ranges)))
    return dict(stations=archived["ports"],charge_flip_potential_span=generator_difference(),
                reversible_leading_pairs_in_both_settings=leading_reciprocity(),
                right_seam_sixth_order=correction_at_sixth_order(),panels=rows)


if __name__ == "__main__":
    print(json.dumps(run(),indent=2,allow_nan=False))
