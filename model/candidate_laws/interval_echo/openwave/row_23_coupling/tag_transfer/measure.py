"""Infer one local tag change and forecast the other tag's remote return records."""

import argparse
import json
from itertools import combinations
from math import factorial
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.sparse.csgraph import shortest_path

from model.candidate_laws.interval_echo import candidate
from model.executable_model import carrier
from model.candidate_laws.interval_echo.openwave.row_23_coupling.experiment import (
    inverse_error,
    norm,
    N,
)

HERE = Path(__file__).resolve().parent
STEP = 1 / 8
NODES, W = np.polynomial.legendre.leggauss(6)
TIMES = (np.arange(48 / STEP)[:, None] * STEP + (NODES + 1) * STEP / 2).ravel()
WEIGHTS = np.tile(W * STEP / 2, int(48 / STEP))
PAIRS = list(combinations(range(4), 2))
FREQUENCIES = [-31 + 0.5j, -10 + 0.5j, -1 + 0.5j, 1 + 0.5j]


def local_blocks(kappa):
    law = candidate.native_law
    ex, co = law.local_term_parts()
    h = ex + kappa * co
    basis = np.eye(law.DIM)
    host = law.SRC.index((0, 0))
    memory = (law.echo() @ basis[:, host]).reshape(-1, law.N_SRC)[:, host]
    odd = np.array([[0.5, -0.5], [-0.5, 0.5]])
    for q in (1, 2):
        rows = [law.SRC.index(w) for w in ((0, q), (q, 0))]
        source = basis[:, rows]
        stored = np.column_stack(
            [np.kron(memory, np.eye(law.N_SRC)[:, r]) for r in rows]
        )
        np.testing.assert_allclose(
            h @ source, source @ (kappa * q * q * np.eye(2) + odd), atol=2e-13
        )
        np.testing.assert_allclose(h @ stored, 0, atol=2e-13)


def setup(level):
    board = carrier.federation(level)
    seam = carrier.EMITTER_SEAM[level]
    a, b = seam
    r = min(set(board.neighbours[a]) & set(board.neighbours[b]))
    dist = shortest_path(
        board.adjacency(), indices=[a, b], directed=False, unweighted=True
    ).min(0)
    shell = (dist == 2).astype(float)
    profiles = np.vstack([np.eye(board.n)[[a, b, r]], shell / np.linalg.norm(shell)])
    np.testing.assert_allclose(profiles @ profiles.T, np.eye(4), atol=1e-14)
    sources = np.stack(
        [*profiles, *[(profiles[i] + profiles[j]) / np.sqrt(2) for i, j in PAIRS]]
    )
    return board, seam, sources


def acquire(board, seam, sources, tag, kappa, written, rng):
    h = board.laplacian().toarray() / 2
    if written:
        h[np.ix_(seam, seam)] -= kappa * tag**2 * np.eye(2) + np.array(
            [[0.5, -0.5], [-0.5, 0.5]]
        )
    e, v = eigh(h)
    # The same-tag uniform reference removes 6*kappa*tag² from both arms.
    amp = np.exp(-1j * TIMES[:, None] * e) @ ((sources @ v) ** 2).T
    probs = np.clip(np.round((1 + np.stack([amp.real, amp.imag])) / 2, 11), 0, 1)
    return rng.binomial(N, probs)


def decode(raw, z):
    sig = 2 * raw / N - 1
    w = -1j * WEIGHTS * np.exp(1j * z * TIMES)
    vals = w @ (sig[0] + 1j * sig[1])
    mat = np.diag(vals[:4])
    for i, (a, b) in enumerate(PAIRS):
        mat[a, b] = mat[b, a] = vals[4 + i] - (vals[a] + vals[b]) / 2
    return mat


def precision(z):
    eta = z.imag
    rate = max(abs(z + 33), abs(z - 6))
    coeff = factorial(6) ** 4 / (13 * factorial(12) ** 3)
    quad = coeff * STEP**13 * rate**12 / (-np.expm1(-eta * STEP))
    tail = np.exp(-48 * eta) / eta
    weights = WEIGHTS * np.exp(-eta * TIMES)
    noise = np.sqrt(2 * np.log(2 * 20000 / 0.01) * (weights @ weights) / N)
    return np.sqrt(2) * (noise + quad + tail) + 2e-11 / eta


def calibrated(raw):
    z = -16 + 1j
    small, big = (decode(raw[f"20_written_8_{q}"], z)[:2, :2] for q in (1, 2))
    shift = np.linalg.inv(small) - np.linalg.inv(big)
    value = shift.trace().real / 2
    err = sum(inverse_error(mat, 3 * precision(z)) for mat in (small, big))
    return dict(
        value=float(value),
        error=float(err),
        scalar_residual=float(norm(shift - value * np.eye(2))),
    )


def forecast(g, delta, uncertainty, eps):
    v = delta * np.eye(2)
    b, x = g[:2, :2], g[2:, :2]
    resolvent = np.linalg.inv(np.eye(2) - v @ b)
    vertex = resolvent @ v
    change = uncertainty * (norm(b) + 3 * eps) + abs(delta) * 3 * eps
    if norm(resolvent) * change >= 1:
        raise ValueError("Unresolved remote forecast")
    derr = (
        norm(resolvent)
        / (1 - norm(resolvent) * change)
        * (uncertainty + change * norm(vertex))
    )
    xerr = 2 * np.sqrt(2) * eps
    out = g[2:, 2:] + x @ vertex @ x.T
    radii = np.array([[eps, 2 * eps], [2 * eps, eps]])
    for i in range(2):
        for j in range(2):
            ni, nj = np.linalg.norm(x[i]), np.linalg.norm(x[j])
            radii[i, j] += (ni + xerr) * (nj + xerr) * derr + xerr * (
                ni + nj + xerr
            ) * norm(vertex)
    return out, radii


def summarize(raw):
    fit = calibrated(raw)
    rows = []
    for n in (20, 80):
        for stock, k in (("blank", 8), ("written", 8), ("written", 4)):
            tag = f"{n}_{stock}_{k}"
            delta = 0 if stock == "blank" else fit["value"] * k / 8
            unc = 0 if stock == "blank" else fit["error"] * k / 8
            for z in FREQUENCIES:
                eps = precision(z)
                a, b = [decode(raw[f"{tag}_{q}"], z) for q in (1, 2)]
                idx = (np.array([0, 1, 0]), np.array([0, 1, 1]))
                target = b[2:, 2:]
                bounds = np.array([[eps, 2 * eps], [2 * eps, eps]])
                pred, r = forecast(a, delta, unc, eps)
                mismatch = abs(pred - target)[idx]
                copy = abs(a[2:, 2:] - target)[idx]
                tests = dict(
                    residual=mismatch.tolist(),
                    joint_radius=(r + bounds)[idx].tolist(),
                    copy_residual=copy.tolist(),
                    copy_excluded=(copy > (2 * bounds)[idx]).tolist(),
                )
                if k == 4:
                    foil, rf = forecast(a, fit["value"], fit["error"], eps)
                    tests["full_contrast_excluded"] = (
                        abs(foil - target) > (rf + bounds)
                    )[idx].tolist()
                rows.append(
                    dict(cells=n, stock=stock, contrast=k, frequency=z.real, **tests)
                )
    return dict(calibration=fit, columns=["neighbor", "shell_two", "cross"], held=rows)


def run():
    rng = np.random.default_rng(2351108)
    banks = {}
    for level in (0, 1):
        board, seam, sources = setup(level)
        for stock, k in (("blank", 8), ("written", 8), ("written", 4)):
            local_blocks(k)
            for q in (1, 2):
                banks[f"{board.cells}_{stock}_{k}_{q}"] = acquire(
                    board, seam, sources, q, k, stock == "written", rng
                )
    return summarize(banks), banks


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    path = parser.parse_args().output
    result, raw = run()
    if path:
        path.mkdir(exist_ok=False)
        np.savez_compressed(path / "counts.npz", **raw)
        (path / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    else:
        from model.candidate_laws.interval_echo.openwave.electron_field.replay import (
            compare,
        )

        compare(result, json.loads((HERE / "results.json").read_text()), "tag transfer")
        with np.load(HERE / "counts.npz") as saved:
            if set(saved.files) != set(raw):
                raise ValueError("Acquisition bank names differ")
            for k in raw:
                np.testing.assert_array_equal(saved[k], raw[k])
    print(result["calibration"])
    for r in result["held"]:
        print(
            r["cells"],
            r["stock"],
            r["contrast"],
            r["frequency"],
            max(np.array(r["residual"]) / r["joint_radius"]),
            r["copy_excluded"],
        )
