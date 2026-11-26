import json
import numpy as np

from . import measure as m


def test_compiled_signed_flow_and_parked_hypercharge_reference():
    for axis in ((0, 0, 0, 0), (0, 0, -1, 0), (-0.4, 0.7, -1.2, 0.2), (0, 0, 0, -1)):
        for t in (-0.47, 0, 0.63):
            actual = m.execute(m.field_word(axis, t))
            np.testing.assert_allclose(actual, m.exterior(axis, t), atol=3e-13, rtol=0)
            np.testing.assert_allclose(actual[:, 2], [0, 0, 1], atol=3e-13, rtol=0)
    psi = np.array([1, 1j, 0]) / np.sqrt(2)
    u = m.execute(m.field_word((0, 0, 0, 1)))
    np.testing.assert_allclose(
        abs(m.reading_vector(("pair", 0, 1, 0.72)) @ u @ psi) ** 2,
        abs(m.reading_vector(("pair", 0, 1, 0.72)) @ psi) ** 2,
    )
    ref = np.array([1, 0, 1]) / np.sqrt(2)
    read = m.reading_vector(("pair", 0, 2, np.pi / 2))
    assert abs((abs(read @ u @ ref) ** 2 - abs(read @ ref) ** 2)[0]) > 0.2


def test_local_analyzers_and_changed_field():
    for rd in m.READS:
        actual = m.execute(m.reader(rd))
        target = m.reading_vector(rd)
        for a, b in zip(actual, target, strict=True):
            np.testing.assert_allclose(
                np.outer(a.conj(), a), np.outer(b.conj(), b), atol=5e-13
            )
    saved = json.loads((m.HERE / "results.json").read_text())
    assert m.report() == saved
    assert all(
        sum(row["counts"]) == m.SHOTS
        for panel in saved["panels"]
        for row in panel["banks"]
    )
    assert saved["matched_totals"] == [245, 210, 234]
    assert saved["held_changed_law_excluded"] == 42
