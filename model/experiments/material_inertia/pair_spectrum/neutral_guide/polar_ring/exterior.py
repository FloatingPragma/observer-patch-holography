"""Enclose the probability of an exterior vacancy on the native stocked material."""

import argparse
from collections import defaultdict
from fractions import Fraction
from hashlib import sha256
import json
from math import factorial, isqrt
from pathlib import Path

import numpy as np

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.row_22_confinement.mobile_string import (
    commission,
)
from model.executable_model import federation
from . import stock
from .measure import microscopic

HERE = Path(__file__).resolve().parent
RESOLUTION = 10**12
DELAYS = (128, 64)
STOCKS = {"host": (0, 1), "reversal": (1, 0), "vacancy": (0, 0), "excess": (1, 1)}


def check_contacts():
    law = candidate.native_law
    write = law._write_pattern()
    reversal, squared = law._label_ops()
    idle = np.array([d == 0 for _, d in law.LABELS], dtype=int)
    columns = np.kron(idle[:, None], np.eye(9, dtype=int))
    # Twice the written component after the echo, with the common normalization omitted.
    written = 2 * columns - write @ (write.T @ columns)
    identity = np.eye(9, dtype=int)
    assert np.array_equal(np.kron(reversal, identity) @ written, written)
    assert not np.any(np.kron(squared, identity) @ written)
    h0, h1 = law.local_term_parts()
    blank = np.eye(law.DIM)[:, :9]
    swap, costs = law._load_ops()
    np.testing.assert_allclose(2 * h0 @ blank, blank @ (identity - swap), atol=1e-14)
    np.testing.assert_allclose(h1 @ blank, blank @ costs, atol=1e-14)
    # The sum of the two symmetric contrast-two inputs writes this connector at 11.
    source = blank[:, law.SRC.index((0, 2))] + blank[:, law.SRC.index((2, 0))]
    memory = (law.echo() @ source).reshape(-1, 9)[:, law.SRC.index((1, 1))]
    for a in (0, 1):
        for b in (0, 1):
            col = np.kron(memory, identity[:, law.SRC.index((a, b))])
            np.testing.assert_allclose(h0 @ col, 0, atol=1e-14)
            np.testing.assert_allclose(h1 @ col, 4 * (1 - a * b) * col, atol=1e-14)
    native = microscopic()
    np.testing.assert_allclose(native["contact_calibration"], [32, 1, 8], atol=1e-12)


def bit_ports(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


class BoardTerms:
    """Twice the binary generator, relative to the completely filled word."""

    def __init__(self, board, dimers, connectors, opened):
        ring = {p for pair in dimers for p in pair}
        edges = {tuple(sorted(edge)) for edge in board.seams}
        blank = {tuple(sorted(pair)) for pair in dimers}
        blank.update(edge for edge in edges if ring.isdisjoint(edge))
        if opened:
            assert (
                tuple(sorted(opened)) in edges and sum(p in ring for p in opened) == 1
            )
            blank.add(tuple(sorted(opened)))
        self.outside = ((1 << board.n) - 1) ^ sum(1 << p for p in ring)
        self.private = sum(
            1 << i for i, degree in enumerate(board.degree) if degree == 5
        )
        self.neighbors = [[] for _ in range(board.n)]
        for a, b in sorted(blank):
            self.neighbors[a].append(b)
            self.neighbors[b].append(a)
        self.connector = {}
        for a, b in connectors:
            assert a not in self.connector and b not in self.connector
            self.connector[a], self.connector[b] = b, a

    def column(self, holes, exchange=True):
        flips = []
        touched = 0
        for a in bit_ports(holes):
            b = self.connector.get(a)
            if b is not None and (a < b or not (holes >> b) & 1):
                touched += 1
            for b in self.neighbors[a]:
                if not (holes >> b) & 1:
                    flips.append(holes ^ (1 << a) ^ (1 << b))
        diagonal = (
            (16 + int(exchange)) * len(flips)
            + 64 * touched
            - 16 * (self.private & holes).bit_count()
        )
        return diagonal, flips if exchange else []

    def powers(self, initial):
        origin = self.column(initial)[0]
        layers = [{initial: 1}]
        for _ in range(5):
            out = defaultdict(int)
            for word, coefficient in layers[-1].items():
                diagonal, flips = self.column(word)
                out[word] += (diagonal - origin) * coefficient
                for other in flips:
                    out[other] -= coefficient
            layers.append({word: value for word, value in out.items() if value})
        return layers


def root_bounds(value):
    scale = 10**18
    floor = isqrt(value.numerator * scale**2 // value.denominator)
    return Fraction(floor, scale), Fraction(floor + 1, scale)


def boundary_words():
    """Read both idle boundary memories with an addressed terminal echo at unit exterior stock."""
    law = candidate.native_law
    basis = np.eye(law.DIM)
    echo = law.echo()
    parked = law.SRC.index((1, 1))
    idle = (echo @ basis[:, parked]).reshape(-1, 9)[:, parked]
    pointer = {}
    for inner in range(3):
        address = law.SRC.index((inner, 1))
        before = np.kron(idle, np.eye(9)[:, address])
        after = echo @ before
        expected = basis[:, address] if inner == 1 else before
        np.testing.assert_allclose(after, expected, atol=1e-14)
        pointer[inner] = "blank" if inner == 1 else "written"
    return {f"{a}{b}": [pointer[a], pointer[b]] for a in range(3) for b in range(3)}


def vacancy_bounds(layers, outside, inverse_delay):
    # The polynomial and the remainder norm squared are rational on every enrolled word.
    unit = 2 * inverse_delay
    denom = factorial(4) * unit**4
    real, imag = defaultdict(int), defaultdict(int)
    for order in range(5):
        factor = factorial(4) // factorial(order) * unit ** (4 - order)
        component = imag if order % 2 else real
        sign = (1, -1, -1, 1)[order % 4]
        for word, value in layers[order].items():
            if word & outside:
                component[word] += sign * value * factor
    probability = Fraction(
        sum(x * x for x in real.values()) + sum(x * x for x in imag.values()), denom**2
    )
    remainder = Fraction(
        sum(x * x for x in layers[5].values()), (factorial(5) * unit**5) ** 2
    )
    low, high = root_bounds(probability)
    error = root_bounds(remainder)[1]
    low = max(Fraction(0), low - error) ** 2
    high = min(Fraction(1), (high + error) ** 2)
    return [int(low * RESOLUTION), -((-high * RESOLUTION) // 1)]


def charged_scale(readings, delay):
    total = {
        name: [sum(arm[delay][name][i] for arm in readings.values()) for i in (0, 1)]
        for name in STOCKS
    }
    lo = Fraction(total["vacancy"][0] - total["excess"][1], 2)
    hi = Fraction(total["vacancy"][1] - total["excess"][0], 2)
    bracket = (
        total["excess"][1] < total["host"][0] <= total["host"][1] < total["vacancy"][0]
    )
    if lo <= 0:
        return {"charged_bracket": bracket, "reversal_within_five_percent": None}
    diff = [
        total["reversal"][0] - total["host"][1],
        total["reversal"][1] - total["host"][0],
    ]
    ratios = [v / scale for v in diff for scale in (lo, hi)]
    return {
        "charged_bracket": bracket,
        "reversal_charge_units": [float(min(ratios)), float(max(ratios))],
        "reversal_within_five_percent": bracket
        and max(map(abs, ratios)) <= Fraction(1, 20),
    }


def assay(board, dimers, connectors, location):
    ring = {p for pair in dimers for p in pair}
    private, gate = sorted(dimers[location], key=lambda p: board.degree[p])
    assert (board.degree[private], board.degree[gate]) == (5, 6)
    host = sum(1 << p for p in ring if board.degree[p] == 5)
    inputs = {}
    for name, (left, right) in STOCKS.items():
        holes = host & ~(1 << private) & ~(1 << gate)
        inputs[name] = holes | ((1 - left) << private) | ((1 - right) << gate)
    arms, cuts = {}, {}
    closed = BoardTerms(board, dimers, connectors, None)
    assert all(
        not (closed.outside >> p) & 1 for site in ring for p in closed.neighbors[site]
    )
    for label, inside in (("private", private), ("gateway", gate)):
        outside = min(p for p in board.neighbours[inside] if p not in ring)
        cuts[label] = [inside, int(outside)]
        terms = BoardTerms(board, dimers, connectors, cuts[label])
        readings = {str(d): {} for d in DELAYS}
        for name, word in inputs.items():
            assert not word & terms.outside
            assert not terms.column(word, exchange=False)[1]
            layers = terms.powers(word)
            for d in DELAYS:
                readings[str(d)][name] = vacancy_bounds(layers, terms.outside, d)
        arms[label] = readings
    # With no exchange all remaining terms are diagonal on these stocks, even on the opening.
    return {
        "ports": cuts,
        "probability_bounds": arms,
        "comparison": {str(d): charged_scale(arms, str(d)) for d in DELAYS},
    }


def acquire():
    check_contacts()
    rings = stock.charts()
    archive = json.loads((HERE / "records.json").read_text())
    assert json.loads(json.dumps(rings)) == archive["boards"]
    strings = json.loads(
        Path(commission.__file__).with_name("commission.json").read_text()
    )
    panels = []
    for level, size in enumerate((20, 80, 320)):
        board = federation(level)
        for row in [r for r in rings if r["cells"] == size]:
            layout = row["stock"]
            trial = {"cells": size, "ring_owners": row["requested_owners"]}
            if layout:
                trial.update(assay(board, layout["ring"], layout["connectors"], 0))
            else:
                trial["refused"] = "registered circuit absent"
            panels.append(trial)
        layout = commission.material_stock(level, 14)
        assert json.loads(json.dumps(layout)) in strings["stocks"]
        dimers = layout["dimers"]
        connectors = [(left[1], right[0]) for left, right in zip(dimers, dimers[1:])]
        panels.append(
            {"cells": size, "open_dimers": 14, **assay(board, dimers, connectors, 6)}
        )
    return {
        "protocol_sha256": sha256((HERE / "EXTERIOR.md").read_bytes()).hexdigest(),
        "probability_denominator": RESOLUTION,
        "closed_boundary_probability": 0,
        "exchange_off_probability": 0,
        "panels": panels,
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    modes = cli.add_mutually_exclusive_group()
    modes.add_argument("--acquire", type=Path, help="write a fresh simulation record")
    modes.add_argument(
        "--terminal-echo",
        action="store_true",
        help="read all nine two-port words with exterior unit standards",
    )
    args = cli.parse_args()
    if args.terminal_echo:
        print(json.dumps(boundary_words(), indent=2))
        raise SystemExit(0)
    result = acquire()
    if args.acquire:
        with args.acquire.open("x") as out:
            json.dump(result, out, indent=2, allow_nan=False)
            out.write("\n")
    elif result != json.loads((HERE / "exterior.json").read_text()):
        raise ValueError("Exterior readings differ from the deposited probabilities")
    for panel in result["panels"]:
        print(
            panel["cells"],
            panel.get("ring_owners", "open"),
            panel.get("comparison", panel.get("refused")),
        )
