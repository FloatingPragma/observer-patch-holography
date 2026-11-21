"""Measure an information projection on the public records of two overlapping patches.

The supplied controller can read and feed back within its own two-cell patch. The
experiment enumerates exact equiprobable branches of the certified sharp-rail
instrument. The comparison ensemble correlates the two private errors analytically.
"""
from collections import Counter
from fractions import Fraction
from itertools import product
import json
from pathlib import Path

from model.executable_model import federation
from .clock import Registers


class Branch:
    """Return one variate on the desired side of each certified fair threshold."""

    def __init__(self, bits):
        self.bits = iter(bits)

    def random(self):
        return 0.25 if next(self.bits) else 0.75


def sample(car, private, bits, order):
    registers = Registers(car, seed=0)
    if any(abs(rate - 0.5) > 1e-12 for rate in registers.rates.values()):
        raise ValueError("the enrolled source is not fair on both sharp inputs")
    draws = [bits[0]]
    for owner in order:
        draws.extend(bits[1 + 2 * owner:3 + 2 * owner])
    registers.rng = Branch(draws)
    for cell in (0, *private):
        registers.loads[12 * cell + 2] = 2

    shared_digit = (3, 2, 4)
    registers.fair(0, {0})
    records = {}
    for owner in order:
        cell = private[owner]
        region = {0, cell}
        offset = 12 * cell
        scratch = (offset, offset + 1, offset + 10)
        retained = (offset + 6, offset + 7, offset + 10)
        registers.initialize_digits([scratch, retained])
        local = {cell}
        # Exchange the common rail with a vacant rail in the private cell. The
        # native route restores every intermediate port. Undo it after feedback
        # so the next observer receives the original physical overlap record.
        route = registers.route(shared_digit[0], offset + 8)
        for a, b in route:
            registers.swap(a, b, region)
        common = registers.flag(offset + 8, offset + 9, local)
        registers.assign(scratch, registers.fair(cell, local), local)
        second = registers.fair(cell, local)
        error = registers.read(scratch, local) & second
        registers.assign(retained, common ^ error, local)
        for a, b in reversed(route):
            registers.swap(a, b, region)
        records[owner] = (registers.read(shared_digit, region), registers.read(retained, region))
        if records[owner][0] != common:
            raise ValueError("transport failed to restore the shared digit")

    if records[0][0] != records[1][0]:
        raise ValueError("the owners disagree on their overlap")
    return records[0][0], records[0][1], records[1][1]


def accessible_tables(counts):
    tables = []
    for private_column in (1, 2):
        marginal = Counter()
        for row, multiplicity in counts.items():
            marginal[row[0], row[private_column]] += multiplicity
        total = sum(counts.values())
        tables.append([str(Fraction(marginal[pair], total)) for pair in product(range(2), repeat=2)])
    return tables


def private_agreement(counts):
    matches = sum(n for (_, left, right), n in counts.items() if left == right)
    return str(Fraction(matches, sum(counts.values())))


def run():
    car = federation(0)
    private = sorted(b // 12 for a, b in car.external if a // 12 == 0)[:2]
    if len(private) != 2:
        raise ValueError("the shared cell needs two incident neighbours")
    counts = Counter()
    for bits in product(range(2), repeat=5):
        forward = sample(car, private, bits, (0, 1))
        backward = sample(car, private, bits, (1, 0))
        if forward != backward:
            raise ValueError("the two owner orders have different retained digits")
        counts[forward] += 1

    # A comparison distribution for the external joint: both owners share one error.
    common_errors = Counter()
    for shared, first, second in product(range(2), repeat=3):
        private_reading = shared ^ (first & second)
        common_errors[shared, private_reading, private_reading] += 1

    # Minimizing KL to a uniform reference maximizes H(X,Y) = H(X) + H(Y xor X | X).
    # The bound H(Y xor X | X) <= H(Y xor X) is sharp only for independent X and error.
    optimum = ["3/8", "1/8", "1/8", "3/8"]
    tables = accessible_tables(counts)
    if tables != [optimum, optimum] or accessible_tables(common_errors) != tables:
        raise ValueError("the accessible tables differ from the unique information projection")

    return {
        "cells": car.cells,
        "owners": [[0, cell] for cell in private],
        "branches_per_order": sum(counts.values()),
        "owner_orders": 2,
        "constraints": {"shared_one": "1/2", "private_mismatch": "1/4"},
        "reference_per_outcome": "1/4",
        "owner_weights": [1, 1],
        "outcome_order": ["00", "01", "10", "11"],
        "accessible_tables": tables,
        "unique_projected_table": optimum,
        "external_private_agreement": private_agreement(counts),
        "same_tables_with_common_error_private_agreement": private_agreement(common_errors),
    }


if __name__ == "__main__":
    receipt = run()
    expected = json.loads(Path(__file__).with_suffix(".json").read_text())
    if receipt != expected:
        raise SystemExit("the physical projection differs from the stored receipt")
    print(json.dumps(receipt))
