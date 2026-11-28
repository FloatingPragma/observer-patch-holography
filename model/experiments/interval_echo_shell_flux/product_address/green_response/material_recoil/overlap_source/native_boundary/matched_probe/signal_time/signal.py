"""Clock sign returns on two source-sensitive neutral memories of the native board."""

import argparse
from functools import cache
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.linalg import eigh
from scipy.optimize import brentq
from scipy.sparse.linalg import expm_multiply

from .. import measure as m

HERE = Path(__file__).resolve().parent
NAMES = ("empty", "mass+", "paired", "erased", "feedback", "orbit", "dephased")
HANDS = ("near", "far")
ARMS = {
    "mass": ("empty", "mass+"),
    "orbital": ("orbit", "dephased"),
    "erasure": ("paired", "erased"),
    "phase_feedback": ("paired", "feedback"),
}


@cache
def local_four_h(stiffness=8):
    law = m.source.LAW
    h = m.source.formation.emission.four_h(stiffness)
    pairs = ((0, 0), (1, 0), (0, 1), (1, 1))
    blocks = []
    for sign in (0, 1):
        emb = np.zeros((99, 8), dtype=int)
        for i, pair in enumerate(pairs):
            for c in range(2):
                for value, coef in ((-1, 1), (1, (-1) ** sign)):
                    emb[
                        (1 + law.LABELS.index((c, value))) * 9 + law.SRC.index(pair),
                        2 * i + c,
                    ] = coef
        op = emb.T @ h @ emb // 2
        np.testing.assert_array_equal(h @ emb, emb @ op)
        blocks.append(op)
    # Changing reversal parity adds one positive projector with norm one.
    delta = blocks[1] - blocks[0]
    np.testing.assert_array_equal(delta @ delta, 4 * delta)
    return blocks


@cache
def model(total, parities):
    board, source, _, _ = m.fp.make_board()
    marked = (source, (0, 1), (2, 3))
    local = local_four_h()
    sites = list(range(board.n)) if total else [-1]
    dim = len(sites) * 8
    h = np.zeros((dim, dim), dtype=int)
    for i, p in enumerate(sites):
        for cc in range(8):
            src = 8 * i + cc
            if p >= 0:
                h[src, src] += 32 * (6 - int(board.degree[p]))
            for e in board.seams:
                if e not in marked:
                    if p in e:
                        h[src, src] += 34
                        h[8 * (e[0] + e[1] - p) + cc, src] -= 2
                    continue
                k = marked.index(e)
                mask = 1 << k
                c = (cc >> k) & 1
                word = (1 + e.index(p)) if p in e else 0
                col = local[(parities >> k) & 1][:, 2 * word + c]
                for target in np.flatnonzero(col):
                    tw, tc = divmod(int(target), 2)
                    q = e[tw - 1] if tw else p
                    h[8 * sites.index(q) + (cc & ~mask) + mask * tc, src] += col[target]
    np.testing.assert_array_equal(h, h.T)
    indexes = [
        8 * sites.index(p) + c + 2 * nc + 4 * fc
        for c in (0, 1)
        for p in (source if total else (-1,))
        for nc, fc in product(range(2), repeat=2)
    ]
    return h, indexes


def get_arms():
    s = m.stocks()
    states = {
        "empty": s["balanced_mass"][0],
        "paired": s["target"][0],
        "erased": s["target"][1],
        "feedback": s["phase_feedback"][1],
        "orbit": s["exchange_1_1"][0],
        "dephased": s["exchange_1_1"][1],
    }
    full = [np.zeros((2, 2)), np.eye(4) / 4]

    def energy(state):
        val = 0.0
        for n, rho in enumerate(state):
            for par in range(8):
                h, ix = model(n, par)
                val += np.sum(np.kron(rho, np.eye(4) / 4) * h[np.ix_(ix, ix)].T) / 32
        return val

    raw = energy(full) - energy(states["empty"])
    pi = (17 / 8) / raw
    assert 0 < pi < 1
    states["mass+"] = [(1 - pi) * states["empty"][0], pi * full[1]]
    energies = {name: energy(states[name]) for name in NAMES}
    gains = {key: energies[y] - energies[x] for key, (x, y) in ARMS.items()}
    np.testing.assert_allclose(list(gains.values()), 17 / 8, atol=1e-12)
    return states, dict(
        full_ordinary_cost=raw, admission=pi, absolute_energy=energies, gains=gains
    )


def bank(states):
    all_freq = [[], []]
    all_coef = [[], []]
    checked = {}
    for n in (0, 1):
        mats = [model(n, s) for s in range(8)]
        ix = mats[0][1]
        hs = [h / 4 for h, _ in mats]
        parts = [eigh(h) for h in hs]
        for k, hand in enumerate(HANDS):
            mask = 2 << k
            for s in range(8):
                if s & mask:
                    continue
                e, v = parts[s]
                f, w = parts[s ^ mask]
                overlap = v.T @ w
                coeff = np.column_stack(
                    [
                        (
                            overlap
                            * (
                                v[ix].T
                                @ np.kron(states[name][n], np.eye(4) / 4)
                                @ w[ix]
                            )
                        ).ravel()
                        / 4
                        for name in NAMES
                    ]
                )
                all_freq[k].append(abs(e[:, None] - f[None, :]).ravel())
                all_coef[k].append(coeff)
        # Independent direct propagation at two times with the same local sign read.
        for t in (1.0, 4.0):
            initial = np.eye(len(hs[0]))[:, ix]
            waves = [
                expm_multiply(
                    -1j
                    * t
                    * sparse.csr_matrix(h - np.trace(h) / len(h) * np.eye(len(h))),
                    initial,
                )
                * np.exp(-1j * t * np.trace(h) / len(h))
                for h in hs
            ]
            for k, hand in enumerate(HANDS):
                mask = 2 << k
                gram = (
                    sum(
                        (waves[s].conj().T @ waves[s ^ mask]).real
                        for s in range(8)
                        if not (s & mask)
                    )
                    / 4
                )
                vals = np.array(
                    [
                        np.sum(gram * np.kron(states[name][n], np.eye(4) / 4))
                        for name in NAMES
                    ]
                )
                checked[t, k] = checked.get((t, k), 0) + vals
    banks = []
    for k in (0, 1):
        true = np.concatenate(all_freq[k])
        weights = np.vstack(all_coef[k])
        freq, which = np.unique(np.round(true, 11), return_inverse=True)
        grouped = np.stack(
            [np.bincount(which, weights=c, minlength=len(freq)) for c in weights.T],
            axis=1,
        )
        keep = np.max(abs(grouped), axis=1) > 2e-15
        error = 10 * np.sum(
            abs(true - freq[which])[:, None] * abs(weights), axis=0
        ) + np.sum(abs(grouped[~keep]), axis=0)
        banks.append((freq[keep], grouped[keep], float(max(error))))
        assert max(error) < 1e-8
    discrepancy = max(
        np.max(abs(np.cos(banks[k][0] * t) @ banks[k][1] - z))
        for (t, k), z in checked.items()
    )
    assert discrepancy < 1e-9
    return banks, float(discrepancy)


def signal(b, t):
    return np.cos(np.array(t)[..., None] * b[0]) @ b[1]


def crossings(banks):
    ts = np.arange(321) / 32
    data = {name: {} for name in NAMES}
    scan = []
    for k, b in enumerate(banks):
        z = signal(b, ts)
        np.testing.assert_allclose(z[:, 0], np.cos(ts), atol=5e-11)
        # |dZ/dt| <= 1 exactly: two parity generators differ by the
        # norm-one local projector, including foreign contacts. If a cell
        # can reach zero, bound its derivative variation by the spectral
        # absolute second-derivative sum to check monotonicity there.
        dz = -np.sin(ts[:, None] * b[0]) @ (b[1] * b[0][:, None])
        second = np.sum(abs(b[1]) * b[0][:, None] ** 2, axis=0)
        for j, i in np.ndindex(len(ts) - 1, len(NAMES)):
            if max(abs(z[j, i]), abs(z[j + 1, i])) <= 1 / 32 + 1e-8:
                assert abs(dz[j, i]) > second[i] / 32 + 1e-6
        scan.append(z.tolist())
        for i, name in enumerate(NAMES):
            ids = np.flatnonzero(z[:-1, i] * z[1:, i] < 0)
            data[name][HANDS[k]] = [
                float(
                    brentq(
                        lambda t: np.cos(t * b[0]) @ b[1][:, i],
                        ts[j],
                        ts[j + 1],
                        xtol=1e-13,
                    )
                )
                for j in ids
            ]
    return data, dict(times=ts.tolist(), signed_return_by_hand=scan)


def samples(banks, events, recorded=None):
    # A common local grid anchored to the coarse scan rather than the calculated root.
    count = 2**30
    rng = np.random.default_rng(661823)
    grid = {}
    families = 0
    for k, hand in enumerate(HANDS):
        for n in NAMES:
            for i, root in enumerate(events[n][hand][:3]):
                start = np.floor(32 * root) / 32
                t = start + np.arange(-32, 289) / 8192
                grid[k, n, i] = t
                families += len(t)
    eps = np.sqrt(np.log(2 * families / 0.01) / (2 * count))
    rows = {}
    for (k, n, i), t in grid.items():
        b = banks[k]
        idx = NAMES.index(n)
        # Avoid evaluating every other stock in the repeated dense windows.
        probs = (1 + np.cos(t[:, None] * b[0]) @ b[1][:, idx]) / 2
        probs = np.clip(probs, 0, 1)
        label = f"{n}/{HANDS[k]}/{i + 1}"
        if recorded is None:
            raw = rng.binomial(count, probs)
        else:
            raw = np.asarray(recorded["records"][label]["same_sign"])
        assert raw.shape == t.shape and np.issubdtype(raw.dtype, np.integer)
        assert np.all((0 <= raw) & (raw <= count))
        np.testing.assert_array_less(abs(raw / count - probs), eps)
        center = raw / count - 0.5
        orientation = 1 if i % 2 == 0 else -1
        before = np.flatnonzero(orientation * center > eps)
        after = np.flatnonzero(orientation * center < -eps)
        bounds = (
            [None, None]
            if not len(before) or not len(after)
            else [float(t[before[-1]]), float(t[after[0]])]
        )
        if bounds[0] is not None:
            assert bounds[0] <= bounds[1]
        rows[label] = dict(
            first_time=float(t[0]),
            step=1 / 8192,
            same_sign=raw.tolist(),
            bracket=bounds,
        )
    return dict(
        seed=661823,
        shots=count,
        appointments=families,
        probability_half_width=float(eps),
        records=rows,
    )


def run(recorded=None):
    states, energies = get_arms()
    banks, err = bank(states)
    events, scan = crossings(banks)
    shot = samples(banks, events, recorded)
    delays = {}
    for k, (lo, hi) in ARMS.items():
        rows = []
        for i in range(3):
            terms = [(hi, "near", 1), (hi, "far", -1), (lo, "near", -1), (lo, "far", 1)]
            if any(len(events[tag][hand]) <= i for tag, hand, _ in terms):
                rows.append(dict(delta=None, interval=None))
                continue
            delta = sum(s * events[tag][hand][i] for tag, hand, s in terms)
            br = [
                (s, shot["records"][f"{n}/{h}/{i + 1}"]["bracket"]) for n, h, s in terms
            ]
            interval = (
                None
                if any(None in z for s, z in br)
                else [sum(s * z[int((s < 0) ^ j)] for s, z in br) for j in (0, 1)]
            )
            # Constant relative rescaling of the two clock axes requires this
            # near/far duration ratio of ratios to agree at successive events.
            ratio = np.prod([events[tag][hand][i] ** sign for tag, hand, sign in terms])
            ratio_interval = (
                None
                if any(None in z for sign, z in br)
                else [
                    float(np.prod([z[int((sign < 0) ^ j)] ** sign for sign, z in br]))
                    for j in (0, 1)
                ]
            )
            rows.append(
                dict(
                    delta=delta,
                    interval=interval,
                    near_far_duration_ratio=float(ratio),
                    ratio_interval=ratio_interval,
                )
            )
        delays[k] = rows
    # The original bridge still conserves first-owner load count minus its
    # coin. Timing reads at either entirely internal clock preserve it.
    one_q = np.array([1, 0, 0, -1])
    mask = one_q[:, None] == one_q[None, :]
    delta = (states["feedback"][1] - states["paired"][1]) - (
        states["dephased"][1] - states["orbit"][1]
    )
    assert np.max(abs(mask * delta)) < 2e-15
    for n in (0, 1):
        q = np.tile([-(i & 1) for i in range(8)], 24 if n else 1)
        if n:
            q[: 12 * 8] += 1
        for parity in range(8):
            h, _ = model(n, parity)
            i, j = np.nonzero(h)
            np.testing.assert_array_equal(q[i], q[j])
    return dict(
        full_clock_waveform_linear_transfer="orbital dephasing minus its own baseline equals feedback minus paired by bridge conservation",
        protocol_sha256=sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        names=NAMES,
        energy=energies,
        register_crossings=events,
        delay_contrasts=delays,
        direct_check_max_error=err,
        spectrum_consolidation_error=[b[2] for b in banks],
        resolved_lines=[len(b[0]) for b in banks],
        scan=scan,
        pointer_counts=shot,
    )


def replay(archive):
    result = run(recorded=archive["pointer_counts"])
    old, new = dict(archive), json.loads(json.dumps(result))
    for packet in (old, new):
        # Near-degenerate frequencies can fall into different numerical bins.
        # The consolidation and independent propagation errors bound the readout.
        lines = packet.pop("resolved_lines")
        assert len(lines) == len(HANDS) and all(type(n) is int and n > 0 for n in lines)
        assert max(packet["spectrum_consolidation_error"]) < 1e-8
        assert packet["direct_check_max_error"] < 1e-9
    m.fp.check(old, new)
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    path = HERE / "results.json"
    if args.write:
        result = run()
        path.write_text(
            json.dumps(result, separators=(",", ":"), allow_nan=False) + "\n"
        )
    else:
        result = replay(json.loads(path.read_text()))
    print(
        result["energy"],
        result["delay_contrasts"],
        result["direct_check_max_error"],
        result["resolved_lines"],
        flush=True,
    )
