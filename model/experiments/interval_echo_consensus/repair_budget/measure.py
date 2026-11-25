"""Predict routed strict-write counts and measure them with a separate native instrument."""

import argparse
from fractions import Fraction
from itertools import combinations, product
import json
from math import comb, log, sqrt
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo.candidate import native_law
from model.executable_model import CELL_SEAMS, federation

HERE = Path(__file__).resolve().parent
SHOTS = 2**18
WAITS = (0, 1, 2, 3, 5, 9, 16, 32)


def rates(ports, vacancies, doubles):
    if min(vacancies, doubles) < 0 or ports < 2 or vacancies + doubles > ports:
        raise ValueError("Invalid port inventory")
    choices = comb(ports, 2)
    return [
        Fraction((vacancies - j) * (doubles - j), choices)
        for j in range(min(vacancies, doubles))
    ]


def forecast(ports, vacancies, doubles, t):
    if t < 0:
        raise ValueError("Negative attempt count")
    probs = [float(p) for p in rates(ports, vacancies, doubles)]
    trans = np.diag([1 - p for p in probs] + [1.0])
    for j, p in enumerate(probs):
        trans[j, j + 1] = p
    return np.linalg.matrix_power(trans, t)[0]


def terminal_budget(ports, vacancies, doubles):
    """Exact finite mean and a certified 99% integer attempt cap."""
    ps = rates(ports, vacancies, doubles)
    if not ps:
        return dict(mean=Fraction(0), cap=0)
    m = len(ps)
    a = abs(vacancies - doubles)
    mean = sum(1 / p for p in ps)
    # For q=p_last and 0<q<1, (1-q)^ceil(1/q) <= exp(-1) < 1/2.
    # Choose the smallest integer k with coefficient/2^k <= 1/100.
    coeff = Fraction(m * (a + 2), m + a + 1)
    k = 0
    while 2**k < 100 * coeff:
        k += 1
    q = ps[-1]
    cap = 1 if q == 1 else k * ((q.denominator + q.numerator - 1) // q.numerator)
    return dict(mean=mean, cap=cap)


def native_patch():
    seamset = set(CELL_SEAMS)
    ports = next(
        group
        for group in combinations(range(12), 4)
        if sum(edge in seamset for edge in combinations(group, 2)) == 5
    )
    # Select a physical spanning tree; the two other native edges are writable chords.
    tree = []
    components = [{p} for p in ports]
    for a, b in combinations(ports, 2):
        if (a, b) not in seamset:
            continue
        left = next(c for c in components if a in c)
        right = next(c for c in components if b in c)
        if left != right:
            tree.append((ports.index(a), ports.index(b)))
            left.update(right)
            components.remove(right)
    words = list(w for w in product(range(3), repeat=4) if sum(w) == 4)
    material_index = {w: i for i, w in enumerate(words)}
    n = len(words)
    exchange, contrast = native_law.local_term_parts()
    h = (exchange + 8 * contrast)[:9, :9]
    generator = np.zeros((n, n))
    for i, w in enumerate(words):
        for a, b in tree:
            col = native_law.SRC.index((w[a], w[b]))
            for out, (u, v) in enumerate(native_law.SRC):
                if abs(h[out, col]) < 1e-12:
                    continue
                new = list(w)
                new[a], new[b] = u, v
                generator[material_index[tuple(new)], i] += h[out, col]
    echo = native_law.echo()[:, :9].reshape(-1, 9, 9)
    return ports, tree, words, material_index, generator, echo


def operators(words, index, echo, eta):
    n = len(words)
    result = []
    for a, b in combinations(range(4), 2):
        jumps = []
        for k, (_, sign) in enumerate(native_law.LABELS, 1):
            if abs(sign) != 2:
                continue
            jump = np.zeros((n, n))
            for i, w in enumerate(words):
                inp = native_law.SRC.index((w[a], w[b]))
                amplitude = echo[k, native_law.SRC.index((1, 1)), inp]
                if abs(amplitude) < 1e-12:
                    continue
                out = list(w)
                out[a], out[b] = 1, 1
                jump[index[tuple(out)], i] = amplitude * sqrt(eta)
            jumps.append(jump)
        P = sum(j.T @ j for j in jumps)
        quiet = np.diag(np.sqrt(np.maximum(0, 1 - np.diag(P))))
        # Native admitted columns are exact orthogonal unit columns.
        if eta == 1:
            quiet = np.diag((np.diag(P) < 0.5).astype(float))
        np.testing.assert_allclose(P + quiet.T @ quiet, np.eye(n), atol=5e-15, rtol=0)
        result.append((quiet, jumps))
    return result


def acquisition():
    rng = np.random.default_rng(111303)
    ports, tree, words, index, H, echo = native_patch()
    rows = []
    for mode, eta in [
        ("native", 1.0),
        ("attenuated", 0.5),
        ("quiet", 0.0),
        ("repeated", 1.0),
    ]:
        ops = operators(words, index, echo, eta)
        for source in ("sharp", "coherent"):
            # Independent endpoint ports in the admission domain for both preparations.
            psi = np.zeros(len(words), complex)
            psi[index[(2, 0, 0, 2)]] = 1
            if source == "coherent":
                psi[index[(0, 2, 2, 0)]] = 1j
                psi /= np.linalg.norm(psi)
            banks = [
                np.outer(psi, psi.conj()),
                np.zeros_like(H, dtype=complex),
                np.zeros_like(H, dtype=complex),
            ]
            for t in range(max(WAITS) + 1):
                if t in WAITS:
                    prob = np.array([b.trace().real for b in banks])
                    np.testing.assert_allclose(prob.sum(), 1, rtol=0, atol=3e-12)
                    rounded = np.maximum(np.round(prob, 13), 0)
                    counts = rng.multinomial(SHOTS, rounded / rounded.sum())
                    rows.append(
                        dict(
                            mode=mode, source=source, attempts=t, writes=counts.tolist()
                        )
                    )
                if t == max(WAITS):
                    break
                nxt = [np.zeros_like(b) for b in banks]
                for j, b in enumerate(banks):
                    if j == 2:
                        nxt[j] += b
                        continue
                    # A stage-dependent protected flight runs on all three tree edges.
                    u = (
                        np.eye(len(words))
                        if mode == "repeated"
                        else expm(-1j * (1 + t % 3) * (j + 1) * H / 9)
                    )
                    incoming = u @ b @ u.conj().T
                    menu = ops[:1] if mode == "repeated" else ops
                    for quiet, jumps in menu:
                        nxt[j] += quiet @ incoming @ quiet.T / len(menu)
                        nxt[j + 1] += sum(g @ incoming @ g.T for g in jumps) / len(menu)
                banks = nxt
    return dict(ports=list(ports), tree=tree, shots=SHOTS, banks=rows)


def run():
    data = acquisition()
    tolerance = sqrt(log(2 * 3 * len(data["banks"]) / 0.01) / (2 * SHOTS))
    for row in data["banks"]:
        pred = forecast(4, 2, 2, row["attempts"])
        observed = np.array(row["writes"]) / SHOTS
        row["accepts_native_forecast"] = bool(max(abs(pred - observed)) <= tolerance)
    data["simultaneous_half_width"] = tolerance
    data["forecasts"] = []
    for level in range(3):
        ports = federation(level).n
        for v, d in ((1, 1), (2, 2), (9, 3), (ports // 2, ports // 2)):
            result = terminal_budget(ports, v, d)
            mean = result["mean"]
            # Keep the large harmonic sum exact in code; print its integer ceiling in the receipt.
            data["forecasts"].append(
                dict(
                    ports=ports,
                    vacancies=v,
                    doubles=d,
                    mean_ceil=-(-mean.numerator // mean.denominator),
                    cap99=result["cap"],
                )
            )
    return data


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    result = run()
    payload = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.write:
        (HERE / "counts.json").write_text(payload)
    elif (HERE / "counts.json").read_text() != payload:
        raise ValueError("Count replay changed")
    print(len(result["banks"]), "complete count banks")
