"""Source spectroscopy on separate copies predicts an undisturbed two-body release."""

import argparse
from hashlib import sha256
import json
from math import sqrt

import numpy as np
from scipy.linalg import eigh
from scipy.sparse.linalg import expm_multiply

from . import mediated_release as lab

HERE = lab.HERE
READS = 2**48
STEP = 1 / 96


def density_factors(hb):
    energies, modes = eigh(hb.toarray())
    weights = np.exp(-lab.cold.apparatus.BETA * energies)
    thermal = modes * np.sqrt(weights / weights.sum())
    p = (thermal**2).sum(axis=1)
    return dict(equilibrium=thermal, aligned=np.sqrt(p)[:, None]), energies, modes


def fly(h, i, j, shells, factor):
    n = len(shells)
    acc = np.zeros(shells.max() + 1)
    laws = np.zeros((len(lab.RELEASE_TIMES), len(acc)))
    # Purification columns are summed without observing the ordinary address.
    for k in range(0, factor.shape[1], 12):
        psi = factor[j, k : k + 12] / np.sqrt(n - 1)
        hpsi = h @ psi
        local = 2 * (hpsi**2 - psi * (h @ hpsi)).sum(axis=1)
        acc += np.bincount(shells[i], weights=local)
        last = 0
        for index, time in enumerate(lab.RELEASE_TIMES):
            if time != last:
                psi = expm_multiply((-1j * (time - last)) * h, psi)
            last = time
            laws[index] += np.bincount(shells[i], weights=(abs(psi) ** 2).sum(axis=1))
    np.testing.assert_allclose(laws.sum(axis=1), 1, atol=2e-12, rtol=0)
    return acc.cumsum(), laws


def cuts(board, shells):
    result = []
    for r in range(1, max(shells)):
        edges = [
            (u, v) if shells[u] <= r else (v, u)
            for u, v in board.seams
            if (shells[u] <= r) != (shells[v] <= r)
        ]
        result.append(np.array(edges).T)
    return result


def quadrature(factor, energies, modes, crossings, rng):
    coeff = modes.T @ factor
    bank = [[] for _ in crossings]
    for t in (STEP, STEP / 2):
        psi = modes @ (np.exp(-1j * t * energies)[:, None] * coeff)
        rho = psi @ psi.conj().T
        p = rho.diagonal().real
        for row, (u, v) in zip(bank, crossings, strict=True):
            mass = np.mean(p[u] + p[v])
            sign = 2 * np.mean(rho[u, v].imag)
            prob = [(mass + sign) / 2, (mass - sign) / 2, 1 - mass]
            row.append(rng.multinomial(READS, prob).tolist())
    return bank


def current_interval(rows, n_edges):
    assert all(min(x) >= 0 and sum(x) == READS for x in rows)
    y = np.array([(pos - neg) * n_edges / (2 * READS) for pos, neg, _ in rows])
    deriv = (8 * y[1] - y[0]) / (3 * STEP)
    # The cut observable is a sum of signed bond quadratures of norm <=3.
    # Spectrum of the single-load Hamiltonian lies in [-8,6]. Subtracting
    # its midpoint bounds the fifth commutator by 14^5 ||observable||.
    width = (
        3 * sqrt(2 * lab.RISK / READS) * n_edges / (2 * STEP)
        + 3 * 14**5 * STEP**4 / 288
    )
    return np.array([deriv - width, deriv + width])


def read(data, medium):
    board = lab.cold.federation(0)
    n = board.n
    marks, _ = lab.cold.read_map(board, medium["spectral_tape"])
    shells = lab.cold.distances(board, marks[0])
    cut = cuts(board, shells)
    g = np.array(medium["decoded"]["ordinary_neighbor_gain"])
    out = {}
    for kind in ("equilibrium", "aligned"):
        result = []
        for index, (u, v) in enumerate(cut):
            div = np.bincount(u, minlength=n) - np.bincount(v, minlength=n)
            # Every boundary port here is exterior, hence the same live
            # neighbor list applies in both media.
            weight = (board.adjacency() @ div) + div
            ci = {}
            for name in ("blank", "written"):
                count = medium["source_counts"][name]
                value = weight @ count / sum(count)
                error = np.ptp(weight) * np.sqrt(lab.RISK / (2 * sum(count)))
                ci[name] = value + np.array([-1, 1]) * error
            grad = lab.signed_product(g - 0.5, lab.delta(ci["written"], ci["blank"]))
            current = {}
            for name in ("blank", "written"):
                current[name] = current_interval(
                    data["arms"][name + "_" + kind]["preparation_current"][index],
                    len(u),
                )
            contrast = lab.delta(current["written"], current["blank"])
            response = lab.delta(grad, contrast) / (n - 1)
            result.append(
                dict(
                    radius=index + 1,
                    from_counts=response.tolist(),
                    stationary_prediction=(grad / (n - 1)).tolist(),
                    measured_bath_current_derivative_contrast=contrast.tolist(),
                )
            )
        tapes = []
        for name in ("blank", "written"):
            tape = data["arms"][name + "_" + kind]["counts"]
            assert all(
                sum(row) == data["release_trials"] and min(row) >= 0 for row in tape
            )
            masses = np.cumsum(tape, axis=1) / data["release_trials"]
            tapes.append(masses - masses[0])
        observed = tapes[1] - tapes[0]
        error = 4 * np.sqrt(lab.RISK / (2 * data["release_trials"]))
        for row in result:
            r = row["radius"]
            row["observed_enclosed_changes"] = [
                (observed[t, r] + np.array([-1, 1]) * error).tolist()
                for t in range(1, len(tape))
            ]
        out[kind] = result
    return out


def shell_field(medium):
    board = lab.cold.federation(0)
    marks, _ = lab.cold.read_map(board, medium["spectral_tape"])
    shells = lab.cold.distances(board, marks[0])
    g = np.mean(medium["decoded"]["ordinary_neighbor_gain"]) - 0.5
    field = np.zeros(board.n)
    for name, mask in (("blank", []), ("written", marks)):
        count = medium["source_counts"][name]
        p = np.array(count) / sum(count)
        adj = board.adjacency().tolil()
        for u, v in mask:
            adj[u, v] = adj[v, u] = 0
        field += (1 if name == "written" else -1) * g * (adj.tocsr() @ p + p)
    return np.bincount(shells, weights=field).tolist()


def run():
    origin = HERE / "mediated_release.json"
    previous = json.loads(origin.read_bytes())
    board = lab.cold.federation(0)
    mark, _ = lab.cold.read_map(board, previous["spectral_tape"])
    shells = lab.cold.distances(board, mark[0])
    boundary = cuts(board, shells)
    rng = np.random.default_rng(738915)
    arms = {}
    for name, mask in (("blank", []), ("written", mark)):
        hb, _, _ = lab.cold.motion(board, mask, 1, lab.K)
        h, i, j = lab.ordinary.two_loads(board, lab.K, mask)
        factors, e, modes = density_factors(hb)
        for kind, factor in factors.items():
            rho = factor @ factor.T
            p = rho.diagonal()
            density = board.adjacency() @ p + p
            # At an exterior edge no written potential touches either end.
            # Re[H_b,rho] is minus the time derivative of the imaginary fringe.
            comm = hb @ rho - (hb @ rho).T
            expected = [
                (4 * lab.K * sum(density[u] - density[v]) + sum(comm[u, v]))
                / (board.n - 1)
                for u, v in boundary
            ]
            a, laws = fly(h, i, j, shells, factor)
            np.testing.assert_allclose(a[1:-1], expected, rtol=0, atol=1e-11)
            if kind == "equilibrium":
                np.testing.assert_allclose(comm, 0, atol=2e-14)
            arms[name + "_" + kind] = dict(
                acceleration=a.tolist(),
                counts=[rng.multinomial(lab.SHOTS, r / r.sum()).tolist() for r in laws],
                preparation_current=quadrature(factor, e, modes, boundary, rng),
            )
    data = dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (origin, HERE / "UNDISTURBED_MEDIUM.md")
        },
        seed=738915,
        arms=arms,
        times=list(lab.RELEASE_TIMES),
        release_trials=lab.SHOTS,
        quadrature_trials=READS,
        quadrature_step=STEP,
    )
    data["stationary_shell_field_from_counts"] = shell_field(previous)
    data["decoded"] = read(data, previous)
    return data


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    g = cli.add_mutually_exclusive_group()
    g.add_argument("--write", action="store_true")
    g.add_argument("--counts-only", action="store_true")
    args = cli.parse_args()
    path = HERE / "undisturbed_medium.json"
    if args.counts_only:
        data = json.loads(path.read_bytes())
        lab.cold.apparatus.check(
            data["decoded"],
            read(data, json.loads((HERE / "mediated_release.json").read_bytes())),
        )
    else:
        data = run()
        if args.write:
            path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
        else:
            lab.cold.apparatus.check(json.loads(path.read_bytes()), data)
    for k, r in data["decoded"].items():
        print(k, r[0])
