"""Positive formation filters and receiver-only channel on a common-port charge."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.stats import beta as beta_distribution

from .certificate import isolate

from .apparatus import fixture, native_graph, permutations_of_chart

BASE = np.array([[0, 2], [1, -1], [-1, -1]], float) / np.array([np.sqrt(2), np.sqrt(6)])
ONE = [1, 2, 4]
CY = (1, 2, 0)
SWAP = (0, 2, 1)
SHOTS = 2**40


def reduce(h, perms):
    """Project one standard-representation row in each at-most-six-word orbit."""
    n = h.shape[0]
    visited = set()
    rows = []
    cols = []
    vals = []
    s = perms[SWAP]
    c = perms[CY]
    ci = c[c]
    k = 0
    for i in range(n):
        if i in visited:
            continue
        orbit = sorted({int(p[i]) for p in perms.values()})
        visited.update(orbit)

        def mat(perm):
            return np.eye(len(orbit))[:, [orbit.index(perm[t]) for t in orbit]]

        ident = np.eye(len(orbit))
        p = ((ident - mat(s)) / 2) @ ((2 * ident - mat(c) - mat(ci)) / 3)
        assert np.max(abs(p - p.T)) < 1e-13
        ev, v = np.linalg.eigh(p)
        for w in v[:, ev > 0.5].T:
            for j, x in zip(orbit, w):
                if abs(x) > 1e-13:
                    rows.append(j)
                    cols.append(k)
                    vals.append(x)
            k += 1
    q = coo_matrix((vals, (rows, cols)), shape=(n, k)).tocsr()
    small = (q.T @ h @ q).toarray()
    return q, small


def build():
    h = native_graph()
    ps = permutations_of_chart()
    assert connected_components(h)[0] == 1
    for p in ps.values():
        assert (h[p, :][:, p] != h).nnz == 0
    q, small = reduce(h, ps)
    ev, v = eigh(small, driver="evd")
    # The second column of the standard irrep is fixed by the source frame.
    c = ps[CY]
    rows = np.array([64 + b + (r << 3) for b in range(8) for r in ONE])
    # Coherent charge arrival followed by a receiver-only local winding read.
    b0 = q[rows] @ v
    b1 = (q[c[rows]] - q[c[c[rows]]]) @ v / np.sqrt(3)
    beta = BASE[:, 0] @ (q[ONE] @ v)
    # (C^-1 - C)/sqrt(3) maps the first input vector to the second.
    check = BASE[:, 1] @ ((q[c[ONE]] - q[c[c[ONE]]]) @ v) / np.sqrt(3)
    if np.linalg.norm(check + beta) < 1e-10:
        b1 = -b1
        check = -check
    np.testing.assert_allclose(check, beta, atol=4e-13)
    g = np.stack((b0, b1), axis=1).reshape(8, 3, 2, len(ev))
    amplitudes = np.einsum("ro,brik,k->boik", BASE, g, beta)
    return ev, amplitudes, beta, q, v, h, ps


def pt(K):
    # b is the discarded source memory; the source half of the original
    # entangled pair is retained, and the receiving owner sees only local E.
    P = np.einsum("boi,bpj->iojp", K, K.conjugate()) / 2
    J = P.reshape(4, 4)
    val = np.linalg.eigvalsh(P.transpose(2, 1, 0, 3).reshape(4, 4))
    return float(J.trace().real), float(val[0]), J


def scan(ev, A, beta):
    breaks = np.r_[0, np.where(np.diff(ev) > 1e-8)[0] + 1, len(ev)]
    groups = []
    for start, end in zip(breaks, breaks[1:]):
        if ev[start] <= 0 or np.linalg.norm(beta[start:end]) < 1e-7:
            continue
        success, minimum, _ = pt(A[:, :, :, start:end].sum(axis=-1))
        groups.append((minimum, success, int(start), int(end)))
    return sorted(groups)


X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]])
Z = np.diag([1, -1])


def bank(K, clock_success, rng):
    # Input path erasure costs 1/3 and receiver demultiplexing costs 1/3.
    # On each axis retain both detector signs, every missed acceptance, the
    # global filter refusal and the source refusal on original launches.
    success, minimum, J = pt(K)
    rows = []
    shots = SHOTS
    for axis, mat in zip(("X", "Y", "Z"), (X, Y, Z)):
        corr = np.trace(J @ np.kron(mat, mat)).real
        same = (success + corr) / 18
        different = (success - corr) / 18
        p = np.maximum(
            [
                same,
                different,
                clock_success / 3 - success / 9,
                (1 - clock_success) / 3,
                2 / 3,
            ],
            0,
        )
        counts = rng.multinomial(shots, p / sum(p))
        rows.append(dict(axis=axis, counts=counts.tolist()))
    # Simultaneous exact binomial intervals for 6 marginals × 4 banks. The
    # linear Bell witness uses two real correlations and the negative Y one.
    alpha = 0.01 / (2 * 6 * 4)
    bounds = []
    for row in rows:
        bins = []
        for k in row["counts"][:2]:
            low = beta_distribution.ppf(alpha, k, shots - k + 1) if k else 0.0
            high = beta_distribution.isf(alpha, k + 1, shots - k) if k < shots else 1.0
            bins.append((float(low), float(high)))
        bounds.append(bins)
    # W=s/2 - <Φ|J|Φ> = (s-Cxx+Cyy-Czz)/4. Estimate s from the X bank.
    a, b, c = bounds
    low = (2 * a[1][0] + b[0][0] - b[1][1] - c[0][1] + c[1][0]) / 4
    high = (2 * a[1][1] + b[0][1] - b[1][0] - c[0][0] + c[1][1]) / 4
    return dict(
        starts_per_axis=shots,
        ideal_total_yield=success / 9,
        original_PT_minimum=minimum / 9,
        conditional_PT_minimum=minimum / success if success else 0,
        clock_admission_or_unfiltered_unit=clock_success,
        counts=rows,
        original_bell_witness_99=[low, high],
    )


def run():
    ev, A, beta, q, v, h, groups = build()
    scan_rows = scan(ev, A, beta)
    minimum, s, i, j = scan_rows[0]
    assert j == i + 1
    proof = isolate(h, groups, q, v, ev, i)
    # Use the clock setting rounded independently of the transmitted input.
    center = round(float(ev[i]), 8)
    M = 2**25
    step = 1 / 32
    angles = (ev - center) * step / (2 * np.pi)
    f = np.exp(-1j * np.pi * (M - 1) * angles) * np.sinc(M * angles) / np.sinc(angles)
    admit = float(np.vdot(beta * f, beta * f).real)
    K = A @ f
    rng = np.random.default_rng(308591)
    narrow = bank(K, admit, rng)
    # Worst-case polar/eigenvector error on the ideal-pole Bell witness.
    error = proof["bipartite_amplitude_error"]
    witness = float(
        s / 2
        - np.trace(A[:, :, :, i], axis1=1, axis2=2).conj()
        @ np.trace(A[:, :, :, i], axis1=1, axis2=2)
        / 4
    )
    proof["ideal_witness_upper_bound"] = witness + (2 * np.sqrt(s) + error) * error
    assert proof["ideal_witness_upper_bound"] < 0
    nullK = np.zeros((16, 2, 2), complex)
    nullK[:8, :, 0] = K[:, :, 0]
    nullK[8:, :, 1] = K[:, :, 1]
    null = bank(nullK, admit, rng)
    free = {}
    for time in (16, 32):
        free[str(time)] = bank(A @ np.exp(-1j * time * (ev - ev[i])), 1.0, rng)
    return dict(
        material=fixture(),
        physical_dimension=h.shape[0],
        scan_positive_reached=len(scan_rows),
        positive_scope="positive poles with source amplitude norm > 1e-7 in measured standard row",
        selected_standard_index=i,
        selected_energy=float(ev[i]),
        ideal_pole_admission=float(beta[i] ** 2),
        ideal_joint_before_apertures=s,
        ideal_PT_before_apertures=minimum,
        isolated_energy_certificate=proof,
        clock=dict(
            center=center,
            tick=step,
            number_delays=M,
            measured_line_weight=float(abs(beta[i] * f[i]) ** 2 / admit),
        ),
        filtered=narrow,
        dephased_source=null,
        autonomous=free,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    path = parser.parse_args().output
    data = run()
    with path.open("x") as fp:
        json.dump(data, fp, indent=2)
        fp.write("\n")
