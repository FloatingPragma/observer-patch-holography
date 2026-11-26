"""Check the physical code against the component indices of the Lean selection theorem."""

import json
import re

import pytest

from . import material, measure, selection


@pytest.fixture(scope="module")
def bank():
    return measure.run()


def test_whole_receipt_and_controls(bank):
    assert json.loads(json.dumps(bank)) == json.loads(
        (measure.HERE / "results.json").read_text()
    )
    assert len(bank["sharp_flights"]) == 32
    assert all(
        r["sender"] == r["receiver"] and material.decode(r["blocked"]) is None
        for r in bank["sharp_flights"]
    )
    signs = bank["unit_probability_insertion_records"]
    assert {(i, j, base) for i, j, base, *_ in signs} == {
        (i, j, base)
        for i in range(5)
        for j in range(i + 1, 5)
        for base in range(32)
        if not base & ((1 << i) | (1 << j))
    }
    for i, j, base, clear, exterior, ordinary in signs:
        assert clear == "0100"
        assert material.decode(exterior) ^ material.decode(ordinary) == 1 << i
        assert material.decode(ordinary) == base | (1 << j)
    refusals = bank["refused_additions"]
    assert {(j, base) for j, base, _ in refusals} == {
        (j, base) for j in range(5) for base in range(32) if base & (1 << j)
    }
    assert all(raw == "0001" for _, _, raw in refusals)
    with pytest.raises(ValueError):
        selection.census([0] * 32)


def test_every_lean_component_and_ballot(bank):
    lean = (measure.HERE.parents[3] / "Lean/Screen/ExteriorSelection.lean").read_text()

    def vector(name):
        text = re.search(
            rf"def {name} : Fin 10 → [^:]+ :=\s*!\[([^\]]+)\]", lean
        ).group(1)
        return [t.strip() for t in text.split(",")]

    rows = bank["local_census_at_both_stations"]
    native_to_lean = bank["station_to_lean_component_indices"]
    in_lean_order = sorted(enumerate(rows), key=lambda item: native_to_lean[item[0]])
    for actual, formal in (
        ("colour_dim", "colorDim"),
        ("weak_dim", "weakDim"),
        ("charge", "charge"),
    ):
        assert [row[actual] for _, row in in_lean_order] == list(
            map(int, vector(formal))
        )
    c, w, q = (
        list(map(int, vector(name))) for name in ("colorDim", "weakDim", "charge")
    )
    coloured, doublet = (
        [x == "true" for x in vector(name)] for name in ("isColored", "isDoublet")
    )
    conjugate = list(map(int, vector("conj")))
    assert [native_to_lean[r["conjugate"]] for _, r in in_lean_order] == conjugate
    for mask, nonempty, chiral, *sums in bank["ballots_at_both_stations"]:
        chosen = [native_to_lean[k] for k in range(10) if mask & (1 << k)]
        assert bool(nonempty) == bool(chosen)
        assert bool(chiral) == all(conjugate[k] not in chosen for k in chosen)
        assert sums == [
            sum(c[k] * w[k] * q[k] for k in chosen),
            sum(w[k] * q[k] for k in chosen if coloured[k]),
            sum(c[k] * q[k] for k in chosen if doublet[k]),
            sum(c[k] * w[k] * q[k] ** 3 for k in chosen),
        ]
    for parity, mask in enumerate(reversed(bank["survivors"])):
        selected = {
            word
            for j, row in enumerate(rows)
            if mask & (1 << j)
            for word in row["occupations"]
        }
        assert selected == {
            word for word in range(1, 31) if word.bit_count() % 2 == parity
        }
    assert len(bank["held_one_row_changes"]) == 20


def test_block_rotations_leave_the_determinant_line_fixed(bank):
    assert [r["modes"] for r in bank["coherent_block_connectors"]] == [
        [0, 1],
        [1, 2],
        [3, 4],
    ]
    for rotation in bank["rotations_on_each_connector"]:
        vacuum, first, second, full = rotation["probabilities_by_input"]
        assert vacuum == [1, 0, 0, 0] and full == [0, 0, 0, 1]
        assert first == list(reversed(second))
        assert 0 < first[1] < 1 and 0 < first[2] < 1
