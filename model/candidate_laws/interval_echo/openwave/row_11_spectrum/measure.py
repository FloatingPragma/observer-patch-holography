"""Complete composition survey on the registered cell and its native ten-port belt."""

import argparse
from itertools import combinations
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.linalg import eigh

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.quark_winding.belt_walls import belt
from model.executable_model import cell, rotations
from model.candidate_laws.interval_echo.openwave.electron_field.replay import compare

HERE = Path(__file__).parent
SHOTS = 262144
WAITS = 2.0 ** np.arange(-12, 17)
TARGET = (3477 - 1) / (207 - 1)


def word_table(n):
    powers = 3 ** np.arange(n)
    codes = np.arange(3**n)
    return powers, codes, ((codes[:, None] // powers) % 3).astype(np.int8)


def sector(board, p, book):
    powers, codes, all_words = book
    keep = ((all_words == 0).sum(axis=1) == p + 1) & ((all_words == 2).sum(axis=1) == p)
    words, ids = all_words[keep], codes[keep]
    pos = np.arange(len(words))
    cost = ((words.astype(int) ** 2 - 1) @ (6 - board.degree)).astype(int)
    rows, cols = [], []
    for a, b in board.seams:
        delta = words[:, a] - words[:, b]
        cost += delta**2
        changed = delta != 0
        rows.append(pos[changed])
        dest = ids[changed] + delta[changed] * (powers[b] - powers[a])
        cols.append(np.searchsorted(ids, dest))
    rows, cols = np.concatenate(rows), np.concatenate(cols)
    swaps = sparse.csr_matrix(
        (np.full(len(rows), -0.5), (rows, cols)), shape=(len(words), len(words))
    )
    degree = np.bincount(rows, minlength=len(words))
    assert (swaps != swaps.T).nnz == 0
    return words, ids, swaps + sparse.diags(degree / 2), cost, degree


def symmetry_basis(words, ids, powers, group, momentum):
    transformations = np.array(
        [np.searchsorted(ids, words @ powers[list(perm)]) for perm in group]
    )
    rep = np.min(transformations, axis=0)
    rows, cols, vals = [], [], []
    count = 0
    for origin in np.unique(rep):
        orbit, first = np.unique(transformations[:, origin], return_index=True)
        if momentum and (momentum * len(orbit)) % len(group):
            continue
        phase = (
            np.ones(len(orbit))
            if momentum == 0
            else np.exp(2j * np.pi * momentum * first / len(group))
        )
        rows.extend(orbit)
        cols.extend([count] * len(orbit))
        vals.extend(phase / np.sqrt(len(orbit)))
        count += 1
    return sparse.csr_matrix((vals, (rows, cols)), shape=(len(words), count))


def phase_meter(energy, rng):
    ph = round(float(energy), 9) * WAITS
    counts = rng.binomial(
        SHOTS, np.round((1 + np.array([np.cos(ph), np.sin(ph)])) / 2, 8)
    )
    re, im = 2 * counts / SHOTS - 1
    wrapped = np.angle(re + 1j * im)
    value = wrapped[0] / WAITS[0]
    for t, angle in zip(WAITS[1:], wrapped[1:], strict=True):
        value = (angle + 2 * np.pi * np.round((t * value - angle) / (2 * np.pi))) / t
    return dict(real_imag_plus=counts.tolist(), energy=value)


def invariant(x):
    a = np.array(x)
    return float(np.sum(a) / np.sum(np.sqrt(a)) ** 2)


def triples(energies):
    rows = []
    for indic in combinations(np.argsort(energies), 3):
        a, b, c = (energies[i] for i in indic)
        rows.append(
            dict(
                p=[int(i) for i in indic],
                ratios=[float(b / a), float(c / a)],
                koide=invariant((a, b, c)),
                offset_free=(c - a) / (b - a),
            )
        )
    return rows


def scan(board, name, group, rng):
    book = word_table(board.n)
    levels = {8: [], 4: []}
    for p in range((board.n + 1) // 2):
        words, ids, kinetic, cost, degree = sector(board, p, book)
        sharp = int(np.argmin(cost))
        for mom in (0, 1, 2) if name == "belt" else (0,):
            u = symmetry_basis(words, ids, book[0], group, mom)
            t = (u.conj().T @ kinetic @ u).toarray()
            v = (u.conj().T @ sparse.diags(cost, dtype=float) @ u).diagonal().real
            for k in levels:
                h = t.copy()
                h[np.diag_indices_from(h)] += k * v
                e, vec = eigh(h, subset_by_index=(0, 0), check_finite=False)
                state = np.asarray(u @ vec[:, 0]).ravel()
                ev = float(e[0])
                residual = np.linalg.norm(kinetic @ state + (k * cost - ev) * state)
                if residual > 1e-8:
                    raise ArithmeticError(
                        "Sector vector is not an eigenstate of the free generator"
                    )
                if mom == 0:
                    prob = abs(state[sharp]) ** 2
                    accepted = int(rng.binomial(SHOTS, np.round(prob, 12)))
                    lv = dict(
                        p=p,
                        size=len(words),
                        sharp_source=words[sharp].tolist(),
                        spectral_floor=k * int(min(cost)),
                        spectral_ceiling=int(max(k * cost + degree)),
                        min_energy=ev,
                        sharp_overlap=prob,
                        filter_accepted=accepted,
                        formation_read=phase_meter(ev, rng),
                        recoils=[],
                    )
                    levels[k].append(lv)
                else:
                    lv = levels[k][-1]
                    gap = ev - lv["min_energy"]
                    lv["recoils"].append(
                        dict(
                            momentum=mom, energy=ev, gap=gap, read=phase_meter(gap, rng)
                        )
                    )
    result = {}
    for k, rows in levels.items():
        energies = [r["formation_read"]["energy"] for r in rows]
        for r in rows:
            if r["recoils"]:
                a, b = (rec["read"]["energy"] for rec in r["recoils"])
                r["finite_inertia"] = (1 - np.cos(2 * np.pi / board.n)) / a
                r["held_recoil_prediction"] = (
                    a
                    * (1 - np.cos(4 * np.pi / board.n))
                    / (1 - np.cos(2 * np.pi / board.n))
                )
                r["held_recoil_residual"] = b - r["held_recoil_prediction"]
        result[str(k)] = dict(compositions=rows, formation_triples=triples(energies))
    return result


def run():
    rng = np.random.default_rng(110424)
    table = dict(
        shots=SHOTS, waits=WAITS.tolist(), native_strength=candidate.KAPPA, seed=110424
    )
    table["cell"] = scan(cell(), "cell", rotations(), rng)
    table["belt_axis"], table["physical_belt"] = next(iter(belt.belts().items()))
    group = [np.roll(np.arange(10), -i) for i in range(10)]
    table["belt"] = scan(belt.ring(10), "belt", group, rng)
    return table


def portable(data):
    if isinstance(data, dict):
        return {key: portable(value) for key, value in data.items()}
    if isinstance(data, list):
        return [portable(value) for value in data]
    return round(data, 8) if isinstance(data, float) else data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    path = HERE / "results.json"
    fresh = run()
    data = json.loads(json.dumps(fresh, allow_nan=False))
    if parser.parse_args().write:
        path.write_text(json.dumps(data, indent=1) + "\n")
    else:
        # Compare every reported field after rounding only numerical diagnostics.
        compare(portable(data), portable(json.loads(path.read_text())), "spectrum")
    print(
        "Reproduced all six whole-cell compositions and five belt dispersions at both strengths."
    )


if __name__ == "__main__":
    main()
