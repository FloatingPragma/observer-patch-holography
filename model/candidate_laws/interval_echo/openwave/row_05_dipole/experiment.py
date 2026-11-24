"""Stationary signed currents and flux precession of a native compact vacancy."""

import argparse
import json
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.linalg import null_space
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.row_29_maxwell.experiment import (
    geometry,
)
from model.executable_model import cell, federation, rotations
from model.experiments.interval_echo_encounter.conjugation import defect_generator

HERE = Path(__file__).parent
SHOTS = 2**22
TIMES = (0, 4, 8, 16, 32, 64)
FIELDS = (-1 / 32, -1 / 64, 0, 1 / 64, 1 / 32)
CHARGE = -1


def sample(prob, rng):
    p = np.round(np.clip(prob, 0, 1), 12)
    if abs(sum(p) - 1) > 1e-9:
        raise ValueError("Incomplete read")
    return rng.multinomial(SHOTS, p / p.sum()).tolist()


@lru_cache(maxsize=None)
def source(gates):
    board = cell()
    xyz = geometry()
    axis = xyz[list(gates)].sum(axis=0)
    if np.linalg.norm(axis) < 1e-10:
        axis = xyz[gates[0]].copy()
    axis /= np.linalg.norm(axis)
    basis = null_space(
        np.vstack((np.eye(12) + board.adjacency().toarray(), np.eye(12)[list(gates)]))
    )
    moment = np.zeros((12, 12), complex)
    areas = []
    for a, b in board.seams:
        area = np.cross(xyz[a], xyz[b]) @ axis / 2
        areas.append(area)
        moment[b, a] = 0.5j * CHARGE * area
        moment[a, b] = moment[b, a].conjugate()
    values, modes = np.linalg.eigh(basis.T @ moment @ basis)
    plus = basis @ modes[:, -1]
    pivot = np.flatnonzero(abs(plus) > max(abs(plus)) - 1e-10)[0]
    plus *= np.exp(-1j * np.angle(plus[pivot]))
    return plus, np.array(areas), values, moment


def currents(state, areas, rng):
    records = []
    mean = 0.0
    for (a, b), area in zip(cell().seams, areas, strict=True):
        # Eigenvalues +/- 1/2 of the probability-current operator.
        p = abs(state[a] - 1j * state[b]) ** 2 / 2
        m = abs(state[a] + 1j * state[b]) ** 2 / 2
        counts = sample([p, m, 1 - p - m], rng)
        records.append(counts)
        mean += CHARGE * area * (counts[0] - counts[1]) / (2 * SHOTS)
    return dict(counts=records, moment=mean)


def gates_of(board, o):
    return tuple(map(int, np.flatnonzero(board.degree[12 * o : 12 * (o + 1)] == 6)))


def embedding(plus, n):
    modes = np.zeros((n, 2), complex)
    modes[:12, 0], modes[:12, 1] = plus, plus.conj()
    return modes


def turn(plus, gates, rng):
    for perm in rotations():
        p = np.asarray(perm)
        if (
            set(p[list(gates)]) != set(gates)
            or np.all(p == np.arange(12))
            or not np.all(p[p[p]] == np.arange(12))
        ):
            continue
        overlap = np.vdot(plus[p], plus)  # <psi|U|psi> for |a> -> |p[a]>
        out = []
        for z in (overlap.real, overlap.imag):
            out.append(sample([(1 + z) / 2, (1 - z) / 2], rng))
        return dict(
            mapping=p.tolist(), ramsey_xy=out, amplitude=[overlap.real, overlap.imag]
        )
    return None


def frames(amps, rng):
    p, m = amps
    norm = abs(p) ** 2 + abs(m) ** 2
    score = 2 * p.conj() * m
    data = []
    for z in (score.real, score.imag):
        data.append(sample([(norm + z) / 2, (norm - z) / 2, 1 - norm], rng))
    return data


def precession(board, gates, k, isolated, rng):
    plus, areas, _, _ = source(gates)
    h = defect_generator(board, CHARGE, k, 0).astype(complex)
    energy = 3 + 4 * k
    h -= sparse.eye(board.n) * energy
    code = embedding(plus, board.n)
    if np.linalg.norm(h @ code) > 3e-11:
        raise ValueError("Preparation leaves the stationary aperture")
    x = code.sum(axis=1) / np.sqrt(2)
    rows = []
    a, b = np.array(cell().seams).T
    for B in FIELDS:
        delta = -0.5 * np.expm1(1j * CHARGE * B * areas)
        drive = sparse.csr_matrix(
            (np.r_[delta, delta.conj()], (np.r_[b, a], np.r_[a, b])), shape=h.shape
        )
        # Each exposure in this series shares the same generator and preparation.
        wave = expm_multiply(-1j * (h + drive), x, start=0, stop=64, num=17)
        read = []
        for t in TIMES:
            vec = wave[t // 4]
            amps = code.conj().T @ vec
            counts = frames(amps, rng)
            read.append(
                dict(
                    time=t,
                    counts_xy=counts,
                    probability_outside_owner=1 - np.sum(abs(vec[:12]) ** 2),
                    probability_outside_doublet=1 - np.sum(abs(amps) ** 2),
                )
            )
        rows.append(dict(field=B, read=read))
    return dict(owners=board.cells, kappa=k, isolated=isolated, rows=rows)


def census():
    rows = []
    for level in range(3):
        board = federation(level)
        h = defect_generator(board, CHARGE, candidate.KAPPA, 0)
        layouts = {}
        for o in range(board.cells):
            gates = gates_of(board, o)
            layouts.setdefault(gates, []).append(o)
        for g, owners in layouts.items():
            p, _, ev, _ = source(g)
            for o in owners:
                v = np.zeros(board.n, complex)
                v[o * 12 : (o + 1) * 12] = p
                if np.linalg.norm(h @ v - (3 + 4 * candidate.KAPPA) * v) > 2e-11:
                    raise ValueError(
                        "Compact source is not stationary on its federation"
                    )
            rows.append(
                dict(
                    owners=board.cells,
                    addresses=owners,
                    gates=g,
                    dimension=len(ev),
                    moments=ev.tolist(),
                )
            )
    return rows


def run():
    rng = np.random.default_rng(510241)
    main_gates = gates_of(federation(0), 0)
    plus, areas, ev, _ = source(main_gates)
    read = {
        key: currents(state, areas, rng)
        for key, state in [
            ("plus", plus),
            ("minus", plus.conj()),
            ("real", np.sqrt(2) * plus.real),
        ]
    }
    banks = []
    for k in (0, candidate.KAPPA):
        for board in (cell(), federation(0), federation(2)):
            banks.append(precession(board, main_gates, k, board.cells == 1, rng))
    return dict(
        shots=SHOTS,
        seed=510241,
        census=census(),
        owner_zero_gates=main_gates,
        axial_areas=areas.tolist(),
        selected_moments=ev.tolist(),
        current=read,
        rotation=turn(plus, np.array(main_gates), rng),
        precession=banks,
    )


def portable(data):
    if isinstance(data, dict):
        return {k: portable(v) for k, v in data.items()}
    if isinstance(data, (list, tuple)):
        return [portable(v) for v in data]
    if isinstance(data, (float, np.floating)):
        return round(float(data), 10)
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    data = portable(run())
    path = HERE / "results.json"
    if parser.parse_args().write:
        path.write_text(json.dumps(data, indent=1) + "\n")
    elif data != json.loads(path.read_text()):
        raise ValueError("Dipole records differ")
    print(
        "Verified complete-board compact spectra, current and rotation counts and applied-flux precession."
    )


if __name__ == "__main__":
    main()
