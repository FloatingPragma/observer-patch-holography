"""Counterpropagating sharp-hub winding channels under free native evolution."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply
from scipy.stats import beta

from ..live_trimer.measure import LAW, blocks
from .apparatus import fixture, native_graph
from .measure import BASE, ONE, X, Y, Z

HERE = Path(__file__).resolve().parent
TIME = 32
SHOTS = 2**39
AXES = (X, Y, Z)
OUTCOMES = (
    "plus",
    "minus",
    "symmetric_pointer",
    "head_or_return_refusal",
    "no_remote_vacancy",
    "source_refusal",
)


def mark_and_lap():
    """Check both matched and blocked directions at both local apparatuses."""
    f = fixture()
    memo = []
    for terminal in ("hub", "tip"):
        for k, leg in enumerate(f["legs"]):
            h = f[terminal]
            a, c, g = leg[:3] if terminal == "hub" else leg[:-4:-1]
            assert (a, g) == f["marked"][k + 3 * (terminal == "tip")]
            motions = [(h, a), (a, g), (g, c), (c, a), (a, h)]
            ports = (h, a, c, g)
            gates = []
            for left, right in motions:
                marked = (left, right) == (a, g)
                assert (
                    marked or (left, right) in f["blank"] or (right, left) in f["blank"]
                )
                term = blocks()[1, 1] if marked else blocks()["blank"]
                local = np.zeros((8, 8))
                for site in (left, right):
                    for bit in (0, 1):
                        # Route-relative c and ordered pair reverse together in
                        # the physically stored native convention.
                        u, v = sorted((left, right))
                        pair = tuple(int(p != site) for p in (u, v))
                        coin = (bit ^ int(left > right)) if marked else 0
                        native = 9 * coin + LAW.SRC.index(pair)
                        host = 9 * coin + LAW.SRC.index((1, 1))
                        j = 2 * ports.index(site) + bit
                        local[j, j] -= term[host, host] / 8
                        for end in np.flatnonzero(term[:, native]):
                            found = (u, v)[LAW.SRC[end % 9].index(0)]
                            newbit = (end // 9) ^ int(left > right) if marked else bit
                            local[2 * ports.index(found) + newbit, j] += (
                                term[end, native] / 8
                            )
                gates.append(expm(-1j * (np.pi / 8 if marked else np.pi) * local))
            lap = np.eye(8, dtype=complex)
            for gate in gates:
                lap = gate @ lap
            np.testing.assert_allclose(lap[:, 0], -np.eye(8)[:, 1], atol=1e-12)
            np.testing.assert_allclose(lap.conj().T[:, 1], -np.eye(8)[:, 0], atol=1e-12)
            np.testing.assert_allclose(lap[0, 0], 0, atol=1e-12)
            memo.append(
                dict(
                    terminal=terminal,
                    arm=k,
                    matched_source_amplitude=-1,
                    matched_reverse_amplitude=-1,
                )
            )
    return memo


def station_swap(h):
    """This is a chart isomorphism, independent of the evolution and receiving filter."""
    f = fixture()
    swap = {f["hub"]: f["tip"], f["tip"]: f["hub"]}
    for leg in f["legs"]:
        swap.update(zip(leg, reversed(leg), strict=True))
    pp = {p: i for i, p in enumerate(f["ports"])}
    perm = np.array(
        [
            64 * pp[swap[x]] + ((mask & 7) << 3) + (mask >> 3)
            for x in f["ports"]
            for mask in range(64)
        ]
    )
    assert np.array_equal(perm[perm], np.arange(len(perm)))
    assert (h[perm, :][:, perm] != h).nnz == 0
    return perm


def channel(Y, terminal):
    """Input two columns, output only the opposite owner's local E head code."""
    remote = 1 - terminal
    slab = Y[64 * remote : 64 * (1 + remote)].reshape(8, 8, 2)
    # Native memory convention (B,A). Convert to (locally held,discarded,input).
    receiver = slab if remote else slab.transpose(1, 0, 2)
    one = receiver[ONE]
    K = np.einsum("ro,rli->loi", BASE, one)
    return (
        K,
        one.reshape(-1, 2),
        slab.reshape(-1, 2),
        Y[64 * terminal : 64 * (1 + terminal)],
    )


def vec(axis, sign):
    if axis == 2:
        return np.array([1, 0] if sign == 1 else [0, 1], complex)
    return np.array([1, sign * (1 if axis == 0 else 1j)], complex) / np.sqrt(2)


def chpt(J):
    return J.reshape(2, 2, 2, 2).transpose(2, 1, 0, 3).reshape(4, 4)


def acquire(K, one, hub, returned, seed):
    rng = np.random.default_rng(seed)
    # Native carrier, source, final hub and code capture. Rows label discarded
    # source memory; retaining that memory could change the quantum channel.
    J = sum(np.outer(k.T.ravel(), k.T.ravel().conjugate()) for k in K) / 2
    success = float(J.trace().real)
    P = hub.conj().T @ hub
    R = returned.conj().T @ returned
    actual = []
    density = []
    for axis in range(3):
        for sign in (1, -1):
            u = vec(axis, sign)
            omega = np.einsum("lij,j->li", K, u)
            rho = omega.T @ omega.conj()
            qE = float(rho.trace().real)
            q1 = float(np.linalg.norm(one @ u) ** 2)
            qhub = float(np.linalg.norm(hub @ u) ** 2)
            ret = float(np.linalg.norm(returned @ u) ** 2)
            measures = []
            for read, m in enumerate(AXES):
                plus = float(np.trace(rho @ (np.eye(2) + m) / 2).real)
                p = np.array(
                    [
                        plus / 9,
                        (qE - plus) / 9,
                        (q1 - qE) / 9,
                        (qhub - q1 / 3) / 3,
                        (1 - qhub) / 3,
                        2 / 3,
                    ]
                )
                np.testing.assert_allclose(sum(p), 1, atol=1e-12)
                assert min(p) > -1e-13
                counts = rng.multinomial(
                    SHOTS, np.maximum(p, 0) / np.maximum(p, 0).sum()
                )
                measures.append(counts.tolist())
            arr = np.array(measures)
            trace = np.mean(arr[:, 0] + arr[:, 1]) / SHOTS
            means = (arr[:, 0] - arr[:, 1]) / SHOTS
            density.append(
                (trace * np.eye(2) + sum(z * m for z, m in zip(means, AXES))) / 2
            )
            # Separate destructive terminal check without receiver head decoder;
            # the preparation erasure refusal is retained on every launch.
            gross = rng.multinomial(
                SHOTS, [ret / 3, qhub / 3, (1 - ret - qhub) / 3, 2 / 3]
            )
            actual.append(
                dict(
                    input_axis=axis,
                    input_sign=sign,
                    settings=measures,
                    bare_ports_home_remote_else_refused=gross.tolist(),
                )
            )
    cross = (density[0] - density[1] + 1j * (density[2] - density[3])) / 2
    reconstructed = np.block([[density[4], cross], [cross.conj().T, density[5]]]) / 2
    exact_min = float(np.linalg.eigvalsh(chpt(J))[0] / 9)
    read_min = float(np.linalg.eigvalsh(chpt(reconstructed))[0])
    # Marginal Clopper-Pearson intervals resolve the rare accepted bins.
    # 2 outcomes x 3 bases x 6 inputs x 2 launch directions, both tails.
    alpha = 0.01 / (2 * 72)
    bin_radius = 0.0
    for item in actual:
        for counts in item["settings"]:
            for k in counts[:2]:
                lo = beta.ppf(alpha, k, SHOTS - k + 1) if k else 0.0
                hi = beta.isf(alpha, k + 1, SHOTS - k) if k < SHOTS else 1.0
                bin_radius = max(bin_radius, k / SHOTS - lo, hi - k / SHOTS)
    # Trace and Pauli contrasts use sums/differences of two marginals.
    # For one reconstructed state the Frobenius radius is 2 sqrt(2)*bin_radius.
    # The two off-diagonal Choi blocks cost at most twice this per block.
    radius = np.sqrt(10) * np.sqrt(2) * bin_radius
    assert read_min + radius < 0
    assert abs(read_min - exact_min) < radius
    return dict(
        seed=seed,
        shots_per_setting=SHOTS,
        count_columns=OUTCOMES,
        bare_local_return=float(R.trace().real / 2),
        bare_remote_arrival=float(P.trace().real / 2),
        receiver_one_hot_given_launch=float(np.trace(one.conj().T @ one).real / 2),
        receiver_standard_given_launch=success,
        all_launches_code_acceptance=success / 9,
        all_launches_PT_minimum=exact_min,
        observed_PT=read_min,
        PT_radius99=radius,
        observed_upper99=read_min + radius,
        observations=actual,
    ), J


def run():
    f = fixture()
    h = native_graph()
    perm = station_swap(h)
    dim = h.shape[0]
    initial = np.zeros((dim, 4))
    initial[np.array(ONE), :2] = BASE
    initial[64 + 8 * np.array(ONE), 2:] = BASE
    np.testing.assert_array_equal(initial[perm, :2], initial[:, 2:])
    evol = -1j * TIME * h
    y = expm_multiply(evol, initial, traceA=evol.diagonal().sum())
    np.testing.assert_allclose(y.conj().T @ y, np.eye(4), atol=4e-11)
    np.testing.assert_allclose(y[perm, :2], y[:, 2:], atol=4e-11)
    directions = []
    js = []
    for idx, seed in enumerate((857711, 590239)):
        params = channel(y[:, 2 * idx : 2 * idx + 2], idx)
        row, J = acquire(*params, seed)
        directions.append(row)
        js.append(J)
    np.testing.assert_allclose(js[0], js[1], atol=4e-12)
    return dict(
        time=TIME,
        dimension=dim,
        terminal_ports=[int(f["hub"]), int(f["tip"])],
        native_lap_controls=mark_and_lap(),
        forward=directions[0],
        reciprocal=directions[1],
        largest_Choi_discrepancy=float(np.max(abs(js[0] - js[1]))),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    text = json.dumps(run(), indent=2) + "\n"
    if parser.parse_args().write:
        (HERE / "reciprocal.json").write_text(text)
    else:
        print(text, end="")
