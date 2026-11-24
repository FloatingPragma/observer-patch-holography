"""Test the full exterior probe channel and the compact band's recoil."""

import argparse
from collections import defaultdict
import json

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply
from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.executable_model import federation
from model.experiments.interval_echo_encounter.conjugation import (
    defect_generator,
    free_seam,
)
from .experiment import HERE, SHOTS, gates_of, portable, sample, source


def projected_column(board, vacancy, probe, h, exterior):
    # Subtract the all-one host. Both defects must lie on distinct ports.
    shifts = defaultdict(complex)
    shifts[vacancy, probe] = candidate.KAPPA * (
        (6 - board.degree[probe]) * 3 - (6 - board.degree[vacancy])
    )
    edges = {
        tuple(sorted((p, q))) for p in (vacancy, probe) for q in board.neighbours[p]
    }
    loads = {vacancy: 0, probe: 2}
    for a, b in edges:
        i, j = loads.get(a, 1), loads.get(b, 1)
        src = law.SRC.index((i, j))
        shifts[vacancy, probe] += h[src, src]
        dst = law.SRC.index((j, i))
        swap = {a: b, b: a}
        v, p = swap.get(vacancy, vacancy), swap.get(probe, probe)
        if p >= exterior and dst != src:
            shifts[v, p] += h[dst, src]
    return shifts


def verify_intertwiner(board, vec):
    h = free_seam(candidate.KAPPA)
    energy = 3 + 4 * candidate.KAPPA
    hp = defect_generator(board, 1, candidate.KAPPA, 0).tocsc()
    # Pick one representative for each way the exterior unit can meet
    # the source and its gateway, with every exterior port on 20 cells.
    near = sorted(b if a < 12 else a for a, b in board.external if (a < 12) != (b < 12))
    selected = range(12, board.n) if board.cells == 20 else [*near, board.n - 1]
    largest = 0
    for y in selected:
        diff = defaultdict(complex)
        for x, amp in enumerate(vec):
            for i, val in projected_column(board, x, y, h, 12).items():
                diff[i] += amp * val
            diff[x, y] -= energy * amp
            for pos in range(hp.indptr[y], hp.indptr[y + 1]):
                j = hp.indices[pos]
                if j >= 12:
                    diff[x, j] -= hp.data[pos] * amp
        largest = max(largest, max(map(abs, diff.values())))
    if largest > 2e-11:
        raise ArithmeticError("Exterior unit resolves the compact doublet")
    return len(selected), largest


def aperture(board, rng):
    v, _, _, _ = source(gates_of(board, 0))
    count, res = verify_intertwiner(board, v)
    x, y = next(e for e in board.external if e[0] < 12)
    r = next(i for i in board.neighbours[y] if i // 12 == y // 12)
    hp = defect_generator(board, 1, candidate.KAPPA, 0)[12:, 12:]
    wave = np.zeros(board.n - 12, complex)
    wave[y - 12] = 1
    state = expm_multiply(-1j * hp, wave, start=0, stop=16, num=17)
    recorded = []
    for t in (0, 1, 4, 8, 16):
        a, b = state[t, [y - 12, r - 12]]
        p = (abs(a - 1j * b) ** 2 / 2, abs(a + 1j * b) ** 2 / 2)
        probs = [*p, 1 - sum(p)]
        recorded.append(
            dict(
                time=t,
                exact=probs,
                directions_and_absent=[sample(probs, rng) for _ in range(3)],
            )
        )
    # Controlled branch of the external probe on the gateway neighbour,
    # versus another site. Protected parked source and individually held rails.
    c = free_seam(candidate.KAPPA)

    def cost(left, right):
        i = law.SRC.index((left, right))
        return c[i, i]

    delta = cost(0, 2) + cost(1, 1) - cost(0, 1) - cost(1, 2)
    pulses = []
    for d in (0, delta):
        angle = d / 16
        pulses.append(
            [
                sample([(1 + f(angle)) / 2, (1 - f(angle)) / 2], rng)
                for f in (np.cos, np.sin)
            ]
        )
    return dict(
        owners=board.cells,
        aperture=[x, y, r],
        native_columns_checked=count,
        max_residual=res,
        receiver=recorded,
        control_connected_energy=delta,
        sham_and_contact_control_xy=pulses,
    )


def recoil(board, rng):
    x, y = next(e for e in board.external if e[0] < 12)
    other = y // 12
    u = source(gates_of(board, 0))[0]
    v = source(gates_of(board, other))[0]
    start = np.zeros(board.n, complex)
    start[:12] = u / np.sqrt(2)
    start[12 * other : 12 * other + 12] = v / np.sqrt(2)
    h = defect_generator(board, -1, candidate.KAPPA, 0)
    cols = []
    for spin in (-1, 1):
        vec = start.copy()
        vec[12 * other : 12 * other + 12] *= np.exp(1j * spin * np.pi / 3)
        cols.append(vec)
    sharp = np.zeros(board.n, complex)
    sharp[x] = 1
    inputs = np.column_stack([*cols, sharp])
    wave = expm_multiply(
        -1j * (h - (3 + 4 * candidate.KAPPA) * sparse.eye(board.n)),
        inputs,
        start=0,
        stop=16,
        num=17,
    )
    rows = []
    for t in (0, 1, 4, 16):
        cases = []
        for i, vec in enumerate(wave[t].T):
            retain = abs(np.vdot(inputs[:, i], vec)) ** 2
            prob = abs(vec) ** 2
            first = prob[:12].sum()
            second = prob[12 * other : 12 * other + 12].sum()
            cases.append(
                dict(
                    return_counts=sample([retain, 1 - retain], rng),
                    owner_counts=sample([first, second, 1 - first - second], rng),
                )
            )
        rows.append(dict(time=t, sources=cases))
    return dict(owners=board.cells, other=other, rows=rows)


def run():
    rng = np.random.default_rng(510242)
    return dict(
        shots=SHOTS,
        seed=510242,
        exterior=[aperture(federation(level), rng) for level in (0, 2)],
        kick=[recoil(federation(level), rng) for level in (0, 2)],
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    data = portable(run())
    path = HERE / "exterior.json"
    if parser.parse_args().write:
        path.write_text(json.dumps(data, indent=1) + "\n")
    elif json.loads(path.read_text()) != data:
        raise ValueError("Exterior and aperture records differ")
    print(
        "Verified exterior vacancy/probe separation, active controls and flat-band kicks."
    )


if __name__ == "__main__":
    main()
