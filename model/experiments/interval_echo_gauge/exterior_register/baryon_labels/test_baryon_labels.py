"""Compare both material readers to the formal arrays and every Grassmann coefficient."""

import json
import re
from itertools import permutations

import pytest

from .. import material
from . import apparatus, grammar, measure


@pytest.fixture(scope="module")
def bank():
    actual = measure.run()
    stored = json.loads((measure.HERE / "results.json").read_text())
    assert json.loads(json.dumps(actual)) == stored
    return stored


def test_complete_phase_and_spin_records(bank):
    phases, radius = apparatus.phase_read(bank["phase_banks"])
    assert len(bank["phase_banks"]) == 31 * 3 * 2 * 2
    assert len(phases) == 31 * 3 * 2 and radius < 0.038
    for word in range(1, 32):
        for side in ("sender", "receiver"):
            b, ell = apparatus.charges(word)
            assert [phases[word, t, side] for t in ("3B", "L", "opposite_L")] == [
                b,
                ell,
                -ell,
            ]
            if word < 31:
                for tag in ("3B", "L"):
                    assert phases[word, tag, side] == -phases[31 ^ word, tag, side]
    for side in ("sender", "receiver"):
        assert {
            (material.decode(a), material.decode(b))
            for p in bank["sharp_packets"]
            for a, b in [p[side]]
        } == {(w, s) for w in range(32) for s in (0, 1)}
    assert all(
        material.decode(frame) is None
        for p in bank["sharp_packets"]
        for frame in p["blocked"]
    )
    assert len(bank["controlled_phase_routes"]) == 20


def test_every_lean_row_and_charge_law(bank):
    lean = (grammar.ROOT / "Lean/Screen/BaryonDimensionSix.lean").read_text()
    rows = bank["both_stations"]["fields"]
    for j, name in enumerate(("q6", "b3", "ell", "triality", "weak", "left")):
        vector = re.search(rf"def {name} : Fin 10 → ℤ := !\[([^\]]+)\]", lean)
        assert [r[j] for r in rows] == [int(x) for x in vector.group(1).split(",")]
    ballots = bank["both_stations"]["ballots"]
    assert len({tuple(r[0]) for r in ballots}) == 715
    expected = {
        tuple(grammar.NAMES.index(n) for n in row)
        for row in grammar.formal.enumerate_operator_rows_unchecked(
            grammar.formal.field_table()
        )
    }
    assert {tuple(row) for row, _, tests in ballots if all(tests)} == expected
    assert [sum(all(tests[:k]) for _, _, tests in ballots) for k in range(6)] == [
        715,
        554,
        12,
        12,
        12,
        8,
    ]
    for _, b3, ell, changed in bank["both_stations"]["charge_law_control"]:
        assert b3 == 3 * ell and abs(ell) == 1 and changed == -ell


def test_every_remote_exterior_sign_and_contraction_coefficient(bank):
    signs = grammar.exterior_signs(bank["wedge_reference_readings"], True)
    unsigned = grammar.exterior_signs(bank["wedge_reference_readings"], False)
    for order in permutations(range(4)):
        inversions = sum(a > b for i, a in enumerate(order) for b in order[i + 1 :])
        assert signs[order] == (-1) ** inversions and unsigned[order] == 1
    assert bank["duplicate_source_reads"] == [[j, "0001"] for j in range(4)]
    phases, _ = apparatus.phase_read(bank["phase_banks"])
    _, components = grammar.field_census(bank["sharp_packets"], phases, "receiver")
    # Map each measured sharp word and separate spin tag into the producer's
    # component frame. For the dagger sector identify its conjugate basis.
    generators = {
        token: grammar.formal.component_generator(grammar.NAMES[f % 5], c, w, spin)
        for f, group in components.items()
        for token, c, w, spin in group
    }
    ranks = []
    for record in bank["both_stations"]["contractions"]:
        positive = tuple(grammar.NAMES[i % 5] for i in record["row"])
        name = grammar.formal.POSITIVE_CLASSES[positive]
        official, arrangements = grammar.formal.contraction_space(name)
        assert len(record["pairings"]) == len(official)
        for raw, expected, layout in zip(
            record["pairings"], official, arrangements, strict=True
        ):
            # Conjugation interchanges undotted and dotted index contractions.
            left, right = (
                ("right", "left") if min(record["row"]) >= 5 else ("left", "right")
            )
            assert raw["weak"] == layout["weak_pairs"]
            # In these rows only one spin class ever has multiple possible
            # pairings, so conjugating does not reorder the list of patterns.
            assert (
                raw[left] == layout["left_spinor_pairs"]
                and raw[right] == layout["right_spinor_pairs"]
            )
            mapped = {}
            for a, b, c, d, coefficient in raw["signed"]:
                order = [generators[k] for k in (a, b, c, d)]
                swaps = sum(x > y for i, x in enumerate(order) for y in order[i + 1 :])
                mapped[tuple(sorted(order))] = coefficient * (-1) ** swaps
            assert mapped == expected
        ranks.append(record["signed_rank"])
        if name == "QQUE":
            assert record["unsigned_rank"] == 0
    assert ranks == [1] * 8
