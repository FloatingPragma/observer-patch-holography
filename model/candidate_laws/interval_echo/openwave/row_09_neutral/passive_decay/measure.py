"""Ordinary formation energies and complete joint records of a native passive decay."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from .material import columns, commission, hardware

HERE = Path(__file__).resolve().parent
SHOTS = 2**22
WAITS = [0, 1, 2, 4, 8, 16, 32, 64]
ACQUISITION = HERE / "acquisition.npz"


class AcquisitionLaw:
    """Detector laws in draw order, as used at acquisition.

    Linear algebra libraries change the last bits of a propagated probability
    and the binomial sampler can then round its mode differently. A replay
    therefore certifies each recomputed law against the archived binary64 law
    and draws from the archived values with the original seed.
    """

    def __init__(self, archive=None, tolerance=1e-11):
        self.archived = None
        if archive is not None:
            with np.load(archive, allow_pickle=False) as data:
                self.archived = np.split(data["laws"], np.cumsum(data["sizes"])[:-1])
        self.tolerance = tolerance
        self.laws = []

    def __call__(self, law):
        law = np.asarray(law, dtype=float)
        if self.archived is not None:
            if len(self.laws) == len(self.archived):
                raise ArithmeticError("replay draws more often than the acquisition")
            saved = self.archived[len(self.laws)]
            if saved.size != law.size or np.max(abs(saved - law.ravel())) > self.tolerance:
                raise ArithmeticError("recomputed detector law differs from the acquisition law")
            law = saved.reshape(law.shape)
        self.laws.append(law)
        return law if law.ndim else float(law)

    def complete(self):
        if self.archived is not None and len(self.laws) != len(self.archived):
            raise ArithmeticError("replay leaves archived acquisition laws unused")

    def save(self, archive):
        np.savez_compressed(
            archive,
            laws=np.concatenate([law.ravel() for law in self.laws]),
            sizes=np.array([law.size for law in self.laws]),
        )


def hamiltonians(n):
    states = np.arange(2**n)
    occ = (states[:, None] >> np.arange(n - 1, -1, -1)) & 1
    a, b, c = occ[:, :3].T

    def form(size, potential):
        labels = np.arange(2**size)
        matrix = np.diag(potential.astype(float))
        for j in range(size):
            matrix[labels, labels ^ (1 << j)] -= 0.5
        return matrix

    source = form(2, np.array([0, 33, 32, 33]))
    reservoir = form(n - 2, np.array([0, 0, 0, -1]) if n == 4 else np.zeros(2))
    split = np.kron(source, np.eye(len(reservoir))) + np.kron(np.eye(4), reservoir)
    joining = -np.diag(a * c)
    return source, reservoir, split, joining


def physical(board, dimers, heads):
    n = len(dimers)
    states = np.arange(2**n)
    rails = np.ones((len(states), board.n), dtype=int)
    for i, (a, b) in enumerate(dimers):
        rails[:, a] = (states >> (n - 1 - i)) & 1
        rails[:, b] = 1 - rails[:, a]
    # Matching-degree references contribute a scalar, retained in the check.
    potentials = (8 * rails**2) @ (6 - board.degree)
    for (a, b), energy in heads.items():
        potentials += energy * (1 - rails[:, a] * rails[:, b])
    return potentials


def number(x):
    return np.round(x, 11).tolist()


def energies(h, stock):
    e, v = eigh(h)
    weight = abs(v.T @ stock) ** 2
    mean = weight @ e
    return dict(
        mean=float(mean),
        spread=float(np.sqrt(weight @ (e - mean) ** 2)),
        weights=number(weight),
    )


def run(law=None):
    """Acquire every arm. With an archived law, replay its draws."""
    # Flights of at most 64 units on sixteen levels below 36 units accumulate
    # eigenvalue roundoff below 1e-11 in any detector probability, far below
    # one count in 2**22.
    law = AcquisitionLaw() if law is None else law
    board, dimers, heads = hardware()
    n = len(dimers)
    s, r, split, j = hamiltonians(n)
    full = split + j
    p = physical(board, dimers, heads)
    np.testing.assert_allclose(p - p[0], np.diag(full) - full[0, 0])
    es, vs = eigh(s)
    er, vr = eigh(r)
    basis = np.kron(vs, vr)
    nr = len(er)
    delta = es[-1] - es[0]
    rng = np.random.default_rng(26092609)
    arms = []
    # Far X eigenvalue -1 (the isolated, excited end level).
    far_minus = np.kron(np.eye(2 ** (n - 1)), np.array([[0.5, -0.5], [-0.5, 0.5]]))
    for name, parent, initial_line in [
        ("parent", 3, 0),
        ("ground", 0, 0),
        ("lower", 2, 0),
        ("packet", 0, 1),
    ]:
        for attached in (True, False):
            h = full if attached else split
            e, v = eigh(h)
            stock = basis[:, parent * nr + initial_line]
            amp0 = v.T @ stock
            samples = []
            for wait in WAITS:
                psi = v @ (np.exp(-1j * e * wait) * amp0)
                outcomes = abs(basis.T @ psi) ** 2
                outcomes /= outcomes.sum()
                counts = rng.multinomial(SHOTS, law(outcomes)).reshape(4, nr)
                far = np.vdot(psi, far_minus @ psi).real
                samples.append(
                    dict(
                        wait=wait,
                        joint=counts.tolist(),
                        far_minus=int(rng.binomial(SHOTS, law(far))),
                        interaction=round(
                            float(np.vdot(psi, (h - split) @ psi).real), 11
                        ),
                        total_energy=round(float(np.vdot(psi, h @ psi).real), 11),
                    )
                )
            arms.append(
                dict(
                    source=name,
                    joined=attached,
                    initial_energy=energies(h, stock),
                    settings=samples,
                )
            )
    law.complete()
    rel = er - er[0]
    lines = []
    for k in range(3):
        for line in range(1, nr):
            outlet = basis[:, k * nr + line]
            lines.append(
                dict(
                    lower_source_index=k,
                    receiver_line=line,
                    source_release=float(es[3] - es[k]),
                    receiver_formation=float(rel[line]),
                    mismatch=float(es[3] - es[k] - rel[line]),
                    joined_outlet=energies(full, outlet),
                )
            )
    ground_e, ground_v = eigh(full)
    initial = basis[:, 3 * nr]
    outlet = basis[:, 2 * nr + 1]
    low = ground_v[:, :nr]
    tail = dict(
        joined_cold_band_inlet_weight=float(np.linalg.norm(low.T @ initial) ** 2),
        joined_cold_band_outlet_weight=float(np.linalg.norm(low.T @ outlet) ** 2),
        all_time_cold_band_amplitude_bound=float(
            np.sum(abs((low.T @ initial) * (low.T @ outlet)))
        ),
    )

    return dict(
        original_starts=SHOTS,
        waits=WAITS,
        dimers=dimers,
        material_heads=[[a, b, c] for (a, b), c in sorted(heads.items())],
        native_max_column_residual=round(columns(), 13),
        cells=board.cells,
        stock=commission(board, dimers, heads),
        source_levels=number(es - es[0]),
        reservoir_levels=number(rel),
        parent_formation=float(delta),
        light_cut=float(0.1 * delta),
        joined_levels=number(ground_e - ground_e[0]),
        cold_band=tail,
        lines=lines,
        arms=arms,
    )


def compare_records(saved, observed):
    if isinstance(saved, dict):
        assert saved.keys() == observed.keys()
        for key in saved:
            compare_records(saved[key], observed[key])
    elif isinstance(saved, list):
        for a, b in zip(saved, observed, strict=True):
            compare_records(a, b)
    elif isinstance(saved, float):
        np.testing.assert_allclose(saved, observed, atol=5e-10, rtol=1e-11)
    else:
        assert saved == observed


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.output:
        law = AcquisitionLaw()
        out = run(law)
        args.output.write_text(json.dumps(out, indent=1) + "\n")
        law.save(args.output.with_name(ACQUISITION.name))
    else:
        out = run(AcquisitionLaw(ACQUISITION))
        old = json.loads((HERE / "records.json").read_text())
        compare_records(old, out)
        print(
            "Verified every material address, probability count, energy and appointment."
        )
