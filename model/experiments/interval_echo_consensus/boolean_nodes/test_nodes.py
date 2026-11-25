import json
from itertools import combinations, product
import numpy as np

from .material import Machine, balanced_pointer, slots
from .measure import HERE, compile_formula, offer_at
from model.candidate_laws.interval_echo.candidate import native_law
from model.executable_model import CELL_SEAMS


def test_each_native_query_column_exactly_with_integer_write_pattern():
    w = native_law._write_pattern()
    assert np.array_equal(w.T @ w, 2 * np.eye(native_law.N_SRC, dtype=int))
    pointer = np.tile([pair == (1, 1) for pair in native_law.SRC], native_law.N_LAB)
    for x, y in product((0, 2), repeat=2):
        i = native_law.SRC.index((x, y))
        column = w[:, i]
        occupied = np.flatnonzero(column) % native_law.N_SRC
        half = (x + y) // 2
        assert all(native_law.SRC[j] == (half, half) for j in occupied)
        assert np.all(pointer[column != 0] == (x + y == 2))
        undo = w.T @ (column * (pointer == (x + y == 2)))
        assert np.array_equal(undo, 2 * np.eye(native_law.N_SRC, dtype=int)[:, i])
        assert balanced_pointer("native")[x, y] == int(x + y == 2)
        first = balanced_pointer("native")[x, y]
        q = 1 if first else 2 * (1 - balanced_pointer("native")[x, 2])
        assert q == (x + y) // 2


def test_native_transducer_for_all_inputs_and_old_outputs():
    machine = Machine(0, 3)
    for op, args in (("copy", (0,)), ("true", ()), ("false", ()), ("nand", (0, 1))):
        for x, y, z in product(range(2), repeat=3):
            s = machine.fresh((x, y, z))
            snap = s.copy()
            q, r = machine.attempt(s, 2, (op, args), "native")
            expected = {"copy": x, "true": 1, "false": 0, "nand": int(not (x and y))}[
                op
            ]
            assert s[machine.rows[2][0][0]] == 2 * expected and (r == 1) == (
                expected != z
            )
            assert all(
                a == b
                for i, (a, b) in enumerate(zip(s, snap))
                if i not in machine.rows[2][0]
            )


def test_capacity_layout_and_controls():
    for g in combinations(range(12), 3):
        op, sc, c, e = slots(g)
        assert len(set(op + sc + (c,))) == 5 and set(op + sc + (c,)).isdisjoint(g)
        assert op in CELL_SEAMS and sc in CELL_SEAMS
    for bank in json.loads((HERE / "records.json").read_text()):
        runs = bank["runs"]
        offers = sum(r["attempts"] for r in runs)
        assert sum(sum(x) for x in bank["reads"]) == offers
        ports = 12 * bank["owners"]
        assert bank["routing_swaps"] <= 16 * (2 * ports - 3) * offers
        assert sum(r["writes"] for r in runs) == sum(row[1] for row in bank["reads"])
        assert all(r["all_consistent"] for r in runs) == (bank["mode"] == "native")
        if bank["mode"] == "native" and not bank["reverse"]:
            assert max(r["attempts"] for r in runs) <= 2 * len(runs[0]["initial"])


def test_compiler_is_input_independent_and_emits_repeat_occurrences():
    assert compile_formula(("nand", ("var", 0), ("var", 0)), 1) == [
        ("copy", (0,)),
        ("copy", (0,)),
        ("nand", (1, 2)),
    ]


def test_request_calendar_from_every_phase_and_stale_register_on_small_board():
    machine = Machine(0, 5)
    nodes = compile_formula(("nand", ("var", 0), ("var", 1)), 2)
    for x, y in product(range(2), repeat=2):
        target = (x, y, int(not (x and y)))
        for stale in product(range(2), repeat=3):
            for phase in range(len(nodes)):
                loads = machine.fresh((x, y) + stale)
                changes = []
                for step in range(len(nodes)):
                    _, mismatch = offer_at(machine, loads, 2, nodes, phase + step)
                    changes.append(mismatch == 1)
                assert any(changes) == (stale != target)
                if stale == target:
                    snapshot = loads.copy()
                    for offset in range(2 * len(nodes)):
                        assert (
                            offer_at(
                                machine, loads, 2, nodes, phase + len(nodes) + offset
                            )[1]
                            != 1
                        )
                    assert loads == snapshot


def test_published_tickets_replay_each_write_and_the_final_quiet_scan():
    for bank in json.loads((HERE / "records.json").read_text()):
        for run in bank["runs"]:
            n = len(run["initial"])
            schedule = list(range(n))
            if bank["reverse"]:
                schedule.reverse()
            bits = list(run["initial"])
            for index, ticket in enumerate(run["tickets"]):
                assert 0 <= ticket < 9
                if ticket % 3 == 1:
                    bits[schedule[index % n]] ^= 1
            assert bits == run["final"]
            assert len(run["tickets"]) == run["attempts"]
            if bank["mode"] == "native":
                assert run["settled"] and run["all_consistent"]
                assert all(t % 3 != 1 for t in run["tickets"][-n:])
