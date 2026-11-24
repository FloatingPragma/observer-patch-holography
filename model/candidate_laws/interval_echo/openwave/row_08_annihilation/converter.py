"""Time-independent pair-conversion comparison on two independently grounded neutral guides."""

import argparse
import json
import numpy as np
from scipy import sparse
from scipy.linalg import eigh
from scipy.sparse.linalg import expm_multiply

from . import experiment as ex

TIMES = (0.5, 2, 8, 32, 64, 128)


def system(k, n, gain):
    h, hs, _, _, ground = ex.basis(k, n)
    d = len(ground)
    neutral = sparse.kronsum(sparse.csr_matrix(hs[1]), sparse.csr_matrix(hs[1]))
    h0 = sparse.block_diag([h, neutral, neutral], format="csr")
    conn = np.zeros((14, 14))
    for i, p in enumerate(ex.READY):
        conn[p, 12 + i] = conn[12 + i, p] = gain
    return h0 + sparse.kron(conn, sparse.eye(d * d), format="csr"), hs[1], ground


def measurement(psi, h, e0, ground, odd, rng):
    d = len(ground)
    states = psi.reshape(14, d, d)
    neutral = states[12:]
    rest = states[:12]
    remain = float(np.vdot(rest, rest).real)
    excitation = np.eye(d) - np.outer(ground, ground)
    energy = ex.read_joint(neutral, excitation, excitation)
    detectors = ex.read_joint(neutral, odd, odd)
    matter = h - e0 * np.eye(d)
    v = eigh(h)[1]
    eigen_distribution = np.sum(abs(v.T @ neutral @ v) ** 2, axis=0)
    return dict(
        energy_outcomes=ex.sample(np.r_[energy, remain], rng),
        parity_outcomes=ex.sample(np.r_[detectors, remain], rng),
        expectation_left_right=[
            float(np.vdot(neutral, matter @ neutral).real),
            float(np.vdot(neutral, neutral @ matter.T).real),
        ],
        energy_shell_weights=eigen_distribution.round(12).tolist(),
    )


def study(k, n, gain, rng):
    h, wing, g = system(k, n, gain)
    d = len(g)
    s = d**2
    starts = np.zeros((14 * s, 2))
    for j, pair in enumerate(((0, 3), (3, 0))):
        idx = ex.PAIRS.index(pair)
        starts[idx * s : (idx + 1) * s, j] = np.kron(g, g)
    h = h - sparse.eye(h.shape[0]) * h.diagonal().mean()
    energy = np.sum(starts * (h @ starts), axis=0)
    times = []
    states = starts.astype(complex)
    prev = 0
    ground = float(g @ wing @ g)
    for t in TIMES:
        states = expm_multiply((-1j * (t - prev)) * h, states)
        prev = t
        tot = np.sum(states.conj() * (h @ states), axis=0).real
        norm = np.sum(abs(states) ** 2, axis=0)
        assert max(abs(norm - 1)) < 1e-8 and max(abs(tot - energy)) < 1e-7
        times.append(
            dict(
                time=t,
                norm_residual=max(abs(norm - 1)),
                conserved_energy_residual=max(abs(tot - energy)),
                launches=[
                    measurement(states[:, j], wing, ground, g, ex.effects(n), rng)
                    for j in range(2)
                ],
            )
        )
    return dict(
        kappa=k,
        dimers=n,
        coupling=gain,
        wing_energies=(np.linalg.eigvalsh(wing) - ground).round(12).tolist(),
        records=times,
    )


def primitive():
    u = ex.LAW.echo()
    cols = [ex.LAW.SRC.index(t) for t in ((2, 0), (0, 2))]
    w = u[:, cols]
    left = np.eye(ex.LAW.DIM)[:, cols]
    f = w @ left.T
    g = (f + f.T) * abs(ex.term(8)[cols[0], cols[1]])
    dark = np.eye(ex.LAW.DIM)[:, ex.LAW.SRC.index((1, 1))]
    assert np.max(abs(g @ dark)) == 0 and np.max(abs(g @ left - 0.5 * w)) < 1e-15
    return float(np.linalg.norm(g @ dark))


def no_wings(k, gain, rng):
    board, ports, _ = ex.geometry()
    native = ex.collision(k, board, ports)
    full = np.zeros((14, 14))
    full[:12, :12] = native
    for i, index in enumerate(ex.READY):
        full[index, 12 + i] = full[12 + i, index] = gain
    energy, v = eigh(full)
    source = ex.PAIRS.index((0, 3))
    values = []
    for t in TIMES:
        psi = v @ (np.exp(-1j * t * energy) * v[source])
        p = float(np.sum(abs(psi[12:]) ** 2))
        values.append(
            dict(
                time=t,
                extinct_quiet_extinct_bright_survived=ex.sample([p, 0, 1 - p], rng),
            )
        )
    return dict(kappa=k, coupling=gain, measurements=values)


def run():
    rng = np.random.default_rng(81125)
    cases = [(8, 2, 0.5), (8, 3, 0.5), (8, 3, 0), (0, 3, 0.5)]
    return dict(
        seed=81125,
        shots=ex.SHOTS,
        blank_vacuum_coupling=primitive(),
        studies=[study(*case, rng) for case in cases],
        no_wings=[no_wings(8, gain, rng) for gain in (0, 0.5)],
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    data = json.loads(json.dumps(run()))
    path = ex.HERE / "converter.json"
    if args.write:
        path.write_text(json.dumps(data, indent=1, allow_nan=False) + "\n")
    else:
        from .verify import compare

        compare(data, json.loads(path.read_text()), "converter")
        print(
            "Verified unmonitored coincidence conversion and both neutral-wing spectra."
        )


if __name__ == "__main__":
    main()
