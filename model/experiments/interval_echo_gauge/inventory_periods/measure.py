"""Interrogate the determinant phases of a protected register using two material arms."""

import argparse
import json
from itertools import combinations, product
from math import comb, log, sqrt, pi
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo.candidate import native_law
from model.executable_model import federation

HERE = Path(__file__).resolve().parent
SHOTS = 2**28


def character(n0, n1, n2):
    m = n0 + n1 + n2
    d = comb(m, n0) * comb(m - n0, n1)
    a, ra = divmod(d * (n1 + 4 * n2), m)
    b, rb = divmod(d * (17 * n0 * n1 + 65 * n0 * n2 + 17 * n1 * n2), m * (m - 1))
    if ra or rb:
        raise ArithmeticError("Noninteger local trace")
    return d, a, b


def local_blocks(words, physical, edge, k):
    X, C = native_law.local_term_parts()
    h = (X + k * C)[:9, :9]
    index = {w: j for j, w in enumerate(words)}
    mat = np.zeros((len(words), len(words)))
    for j, w in enumerate(words):
        a, b = edge
        s = native_law.SRC.index((w[a], w[b]))
        for z, (u, v) in enumerate(native_law.SRC):
            if abs(h[z, s]) < 1e-12:
                continue
            target = list(w)
            target[a], target[b] = u, v
            mat[index[tuple(target)], j] += h[z, s]
    private = np.diag([k * w[physical] ** 2 for w in words])
    return mat, private


def geometry():
    c = federation(0)
    private = [p for p in range(12) if c.degree[p] == 5]
    for group in combinations(private, 4):
        edges = [
            (i, j)
            for i, j in combinations(range(4), 2)
            if group[j] in c.neighbours[group[i]]
        ]
        # The phase experiment uses an incident powered edge. Require a connected protected region.
        seen = {0}
        for _ in group:
            for i, j in edges:
                if i in seen or j in seen:
                    seen.update((i, j))
        if len(seen) == 4:
            return group, next(e for e in edges if 0 in e)
    raise ValueError("No valid native aperture")


def acquire(u, rng):
    banks = []
    for i, j in product(range(len(u)), repeat=2):
        a = u[i, j]
        for ref in (1, 1j):
            probs = np.array(
                [
                    abs(a + ref) ** 2 / 4,
                    abs(a - ref) ** 2 / 4,
                    max(0, 1 - abs(a) ** 2) / 2,
                ]
            )
            probs = np.round(probs, 13)
            banks.append(rng.multinomial(SHOTS, probs / probs.sum()).tolist())
    return banks


def decode(banks, n):
    bank = np.array(banks).reshape(n, n, 2, 3)
    if np.any(bank.sum(axis=-1) != SHOTS):
        raise ValueError("Missing detector launches")
    diff = (bank[:, :, :, 0] - bank[:, :, :, 1]) / SHOTS
    return diff[:, :, 0] + 1j * diff[:, :, 1]


def run():
    ports, edge = geometry()
    rng = np.random.default_rng(120035)
    rows = []
    for control in ("native", "reverse", "sham", "exchange_only", "commutator"):
        k = 0 if control == "exchange_only" else 8
        for n1, n2 in ((0, 1), (2, 0), (1, 0)):
            words = [
                w
                for w in product(range(3), repeat=4)
                if w.count(1) == n1 and w.count(2) == n2
            ]
            h, q = local_blocks(words, 0, edge, k)
            a = expm(-1j * 2 * pi * 6 / 59 * q / 8)
            b = expm(-1j * 2 * pi * 55 / 59 * h)
            if control == "reverse":
                u = a @ b
            elif control == "sham":
                u = np.eye(len(words))
            elif control == "commutator":
                a = expm(-1j * pi * q / 16)
                b = expm(-1j * pi * h / 2)
                u = a @ b @ a.conj().T @ b.conj().T
            else:
                u = b @ a
            banks = acquire(u, rng)
            est = decode(banks, len(words))
            det = np.linalg.det(est)
            # Two coordinates per entry for every table, each a sample mean of {-1,0,1}.
            eps = sqrt(2 * log(2 * 680 / 0.01) / SHOTS) + 1e-12
            bound = (1 + sqrt(2) * len(words) * eps) ** len(words) - 1
            rows.append(
                dict(
                    control=control,
                    ones=n1,
                    twos=n2,
                    banks=banks,
                    det=[float(det.real), float(det.imag)],
                    joint_radius=bound,
                    max_transition=float(np.max(abs(est - np.diag(np.diag(est))) ** 2)),
                )
            )
    return dict(ports=ports, edge=edge, launches_per_bank=SHOTS, rows=rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    data = run()
    payload = json.dumps(data, separators=(",", ":")) + "\n"
    file = HERE / "records.json"
    if args.write:
        file.write_text(payload)
    elif file.read_text() != payload:
        raise ValueError("Phase readback changed")
    for r in data["rows"]:
        print(
            r["control"],
            r["ones"],
            r["twos"],
            r["det"],
            r["joint_radius"],
            r["max_transition"],
        )
