"""Full-menu phase and material records of a populated odd dimer ring."""

import argparse
from functools import cache
from itertools import product
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.sparse import csr_matrix, diags
from scipy.sparse.linalg import eigsh

from model.candidate_laws.interval_echo import candidate
from model.executable_model.carrier import federation
from .. import experiment as parent
from .geometry import layouts

HERE = Path(__file__).resolve().parent
SHOTS = parent.SHOTS
TIMES = 2.0 ** np.arange(-10, 8)
ENERGY_GRID = 2**20


@cache
def ground(n, signed_fields, omitted, kinetic):
    """Certify the Perron minimum of the actual guarded spin block."""
    j = (4 * candidate.KAPPA + (0.5 if kinetic else 0)) / 2
    h = 0.5 if kinetic else 0
    words = np.arange(2**n)
    flips = words[:, None] ^ (1 << np.arange(n))
    z = 2 * ((words[:, None] >> np.arange(n)) & 1) - 1
    diagonal = -j * (z @ signed_fields)
    for k in range(n):
        if k != omitted:
            diagonal += j * z[:, k] * z[:, (k + 1) % n]
    if not kinetic:
        return float(min(diagonal)), 0.0
    rows = np.repeat(words, n)
    hop = csr_matrix(
        (np.full(rows.shape, -h), (rows, flips.ravel())), shape=(len(words), len(words))
    )
    operator = diags(diagonal) + hop
    if n <= 3:
        eig, vectors = np.linalg.eigh(operator.toarray())
    else:
        eig, vectors = eigsh(operator, k=1, which="SA", tol=1e-13, v0=np.ones(2**n))
    # Repair small far-tail entries before using the exact Collatz inequalities.
    mu = float(eig[0])
    v = np.maximum(abs(vectors[:, 0]), 1e-80)
    if min(diagonal) <= mu:
        raise ValueError("Inconsistent variational ground estimate")
    for _ in range(200):
        ratios = diagonal - h * v[flips].sum(axis=1) / v
        if np.ptp(ratios) < 2e-9:
            break
        v = 0.5 * (v + h * v[flips].sum(axis=1) / (diagonal - mu))
    # Certify the floating proposal with integer inequalities on each row.
    # Every component of this particular positive vector is an exact dyadic.
    scale = 2**34
    low = int(np.floor((min(ratios) - 2e-11) * scale))
    high = int(np.ceil((max(ratios) + 2e-11) * scale))
    components = [x.as_integer_ratio() for x in v]
    denominator = max(b for _, b in components)
    weights = [a * (denominator // b) for a, b in components]
    for i, weight in enumerate(weights):
        neighbor_weight = sum(weights[i ^ (1 << k)] for k in range(n))
        four_hw = int(4 * diagonal[i]) * weight - 2 * neighbor_weight
        if not 4 * low * weight <= scale * four_hw <= 4 * high * weight:
            raise ValueError("Failed exact ground-energy enclosure")
    if (high - low) / scale > 5e-8:
        raise ValueError("Ground interval too wide")
    return (high + low) / (2 * scale), (high - low) / (2 * scale)


def canonical(fields, cut):
    n = len(fields)
    if cut >= 0:
        # Omitting link k opens the chain at k+1; then reflection is allowed.
        ix = cut + 1
        f = fields[ix:] + fields[:ix]
        options = (f, f[::-1])
        cut = n - 1
    else:
        options = tuple(f[i:] + f[:i] for f in (fields, fields[::-1]) for i in range(n))
    return min(x for f in options for x in (f, tuple(-v for v in f))), cut


def energies(layout, charges, cut, kinetic):
    board = federation(0)
    pairs = layout["pairs"]
    fields = [0] * len(pairs)
    for site, q in zip(layout["marks"], charges, strict=True):
        fields[site] = q
    f, cut_rep = canonical(tuple(fields), cut)
    eig, err = ground(len(pairs), f, cut_rep, kinetic)
    ex, contrast = candidate.native_law.local_term_parts()
    local = (ex if kinetic else 0) + candidate.KAPPA * contrast
    index = candidate.native_law.SRC.index

    def energy(a, b):
        value = local[index((a, b)), index((a, b))]
        half = round(2 * value) / 2
        if abs(half - value) > 1e-13:
            raise ValueError("Native diagonal changed")
        return half

    offset = 0.0
    for k, (a, b) in enumerate(pairs):
        if board.degree[a] != board.degree[b]:
            raise ValueError("Unequal dimer on-site energies")
        offset += energy(2, 0) + 4 * candidate.KAPPA * (6 - board.degree[a])
        if k != cut:
            offset += (energy(0, 2) + energy(0, 0)) / 2
    for pin, q in zip(layout["pins"], charges, strict=True):
        offset += (energy(q + 1, 0) + energy(q + 1, 2)) / 2
        offset += candidate.KAPPA * (6 - board.degree[pin]) * (q + 1) ** 2
    return float(offset + eig), err


def two_pins(m, a, b):
    out = []
    for x, y in product((-1, 0, 1), repeat=2):
        signs = [0] * m
        signs[a], signs[b] = x, y
        out.append(tuple(signs))
    return out


def public_instruments(layout, rng):
    """Sacrificial local sum reads and independent token ranging with a held-link control."""
    native = candidate.native_law
    source = native.SRC
    u = native.echo()
    ports = []
    for pin, standard in zip(layout["pins"], layout["standards"], strict=True):
        menu = []
        for q in (-1, 0, 1):
            column = u[:, source.index((q + 1, 1))].reshape(-1, len(source))
            p = (abs(column) ** 2).sum(axis=0)
            dyadic = np.rint(2 * p) / 2
            np.testing.assert_allclose(p, dyadic, rtol=0, atol=1e-13)
            menu.append(rng.multinomial(16384, dyadic).tolist())
        ports.append(dict(pin=pin, standard=standard, terminal_words=menu))
    ex, contrast = native.local_term_parts()
    full = expm(-1j * np.pi * (ex + candidate.KAPPA * contrast))
    rails = [source.index((2, 0)), source.index((0, 2))]
    gate = full[np.ix_(rails, rails)]
    np.testing.assert_allclose(gate, [[0, 1], [1, 0]], rtol=0, atol=1e-13)
    routes = {}
    for tag, route in layout["routes"].items():
        arms = []
        for skip in (-1, (len(route) - 2) // 2):
            token = np.zeros(federation(0).n, complex)
            token[route[0]] = 1
            reads = []
            for k, (a, b) in enumerate(zip(route, route[1:])):
                if k != skip:
                    token[[a, b]] = gate @ token[[a, b]]
                prob = abs(token[route[-1]]) ** 2
                sharp = round(prob)
                if abs(sharp - prob) > 1e-12:
                    raise ValueError("Nonsharp ranging transfer")
                reads.append(int(rng.binomial(4096, sharp)))
            arms.append(reads)
        routes[tag] = arms
    return dict(meters=ports, ranging=routes)


def run():
    rng = np.random.default_rng(20112532)
    instrument_rng = np.random.default_rng(520119)
    boards = []
    for layout in layouts():
        n = len(layout["pairs"])
        m = len(layout["marks"])
        triplets = [u + (0,) * (m - 3) for u in product((-1, 0, 1), repeat=3)]
        if m == 4:
            triplets += [u for u in two_pins(m, 0, 3) if u not in triplets]
        pair = two_pins(m, 0, 1)
        banks = []
        for kinetic, cut, menu in (
            (True, -1, triplets),
            (True, 0, pair),
            (True, n - 1, pair),
            (False, -1, triplets),
        ):
            panels = []
            for charges in menu:
                e, err = energies(layout, charges, cut, kinetic)
                grid = round(e * ENERGY_GRID)
                if abs(e - grid / ENERGY_GRID) + err > 0.5 / ENERGY_GRID:
                    raise ValueError("Energy interval crosses the portable Ramsey mesh")
                if abs(e) * TIMES[0] > np.pi - 0.01:
                    raise ValueError("The initial local clock can wrap")
                read = parent.phase_reader(grid / ENERGY_GRID, rng, TIMES)
                panels.append(
                    dict(q=charges, energy_tick=grid, populations=read["counts_xy"])
                )
            banks.append(dict(kinetic=kinetic, cut=cut, panels=panels))
        boards.append(
            dict(
                layout=layout,
                banks=banks,
                instruments=public_instruments(layout, instrument_rng),
            )
        )
    return dict(
        shots=SHOTS, clock=TIMES.tolist(), energy_grid=ENERGY_GRID, rings=boards
    )


def decode(panel):
    z = 2 * np.array(panel["populations"]) / SHOTS - 1
    phases = np.arctan2(z[1], z[0])
    phase = phases[0]
    for p in phases[1:]:
        phase = p + 2 * np.pi * round((2 * phase - p) / (2 * np.pi))
    return phase / TIMES[-1]


def summary(bank):
    e = {tuple(p["q"]): decode(p) for p in bank["panels"]}
    m = len(next(iter(e)))
    neutral = (0,) * m

    def state(*values):
        q = list(neutral)
        for k, v in values:
            q[k] = v
        return tuple(q)

    pairs = {}
    for i in range(m):
        for j in range(i + 1, m):
            corners = [state((i, a), (j, b)) for a, b in product((-1, 1), repeat=2)]
            if all(k in e for k in corners):
                pair = [
                    e[q] - e[state((i, q[i]))] - e[state((j, q[j]))] + e[neutral]
                    for q in corners
                ]
                pairs[f"{i}{j}"] = dict(
                    odd=float(np.dot([1, -1, -1, 1], pair) / 4),
                    even=float(np.mean(pair)),
                    U=pair,
                )
    triple = {}
    for s in product((-1, 1), repeat=3):
        full = state(*enumerate(s))
        if full in e:
            total = e[full] - e[neutral]
            for i in range(3):
                total += e[state((i, s[i]))]
                total -= e[state(*[(j, s[j]) for j in range(3) if i != j])]
            triple[str(s)] = total
    return pairs, triple


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = run()
    out = json.dumps(data, separators=(",", ":"), allow_nan=False) + "\n"
    p = HERE / "records.json"
    if args.write:
        p.write_text(out)
    elif p.read_text() != out:
        raise SystemExit("Static ring records differ")
    for r in data["rings"]:
        for b in r["banks"]:
            a, t = summary(b)
            print(
                len(r["layout"]["pairs"]),
                b["kinetic"],
                b["cut"],
                {k: round(v["odd"], 6) for k, v in a.items()},
                (min(t.values()), max(t.values())) if t else (),
            )
