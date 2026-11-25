import json

import numpy as np
from scipy.linalg import expm

from . import measure as m


def test_independent_free_step_and_registered_counts():
    saved = json.loads((m.HERE / "results.json").read_text())
    words, h = m.sector(10, 2, 1, 8)
    ix = {w: i for i, w in enumerate(words)}
    independent = np.zeros_like(h)
    for w, i in ix.items():
        for j in range(9):
            a, b = w[j : j + 2]
            swapped = (*w[:j], b, a, *w[j + 2 :])
            independent[i, i] += (a - b) ** 2 * 8 + 0.5
            independent[ix[swapped], i] -= 0.5
    np.testing.assert_array_equal(independent, h)
    f = np.array([[1, -1], [1, 1]]) / np.sqrt(2)
    frames = [f, f]
    start = m.frame_columns(words, m.SOURCE, frames)[0][:, 2]
    target = m.frame_columns(words, m.DAUGHTER, frames)[0][:, 1]
    p = abs(target @ np.linalg.matrix_power(expm(-1j * h / 2), 1024) @ start) ** 2
    sample = next(x for x in saved["raised"]["records"] if x["exposure"] == 512)
    assert abs(p - sample["probability"]) < 1e-10
    assert abs(p - sample["accepted"] / m.SHOTS) < 1e-3
    coeff = m.exact_transition(words, h, start, target)
    assert coeff == saved["raised"]["integer_moments"] == [0] * 7 + [3185]


def test_complete_records_and_control_samples():
    result = m.run()
    archived = json.loads((m.HERE / "results.json").read_text())
    # Long free exposures permit tiny LAPACK differences between machines.
    for bank in ("raised", "ground", "sealed", "reverse"):
        for live, stored in zip(
            result[bank]["records"], archived[bank]["records"], strict=True
        ):
            assert abs(live["probability"] - stored["probability"]) < 2e-8
            np.testing.assert_allclose(live["counts"], stored["counts"], atol=2, rtol=0)
            np.testing.assert_allclose(
                live["receiver_inventory"], stored["receiver_inventory"], atol=2, rtol=0
            )
    for field in ("calibration", "geometry", "pair_alone", "head_only"):
        assert json.loads(json.dumps(result[field])) == archived[field]
    for name in ("raised", "ground", "sealed", "reverse"):
        bank = result[name]
        for row in bank["records"]:
            assert np.sum(row["counts"]) == np.sum(row["receiver_inventory"]) == m.SHOTS
    assert all(row["accepted"] == 0 for row in result["sealed"]["records"])
    assert result["pair_alone"]["all_time_probability"] == 0
    assert (
        result["raised"]["source_mean_variance"]
        == result["raised"]["target_mean_variance"]
        == [66.5, 0.75]
    )
    assert result["raised"]["records"][7]["probability"] > 0.15
    assert result["ground"]["records"][7]["probability"] < 0.004


def test_head_energy_is_observable_but_conserved_on_filled_seams():
    law = m.candidate.native_law
    x, y = law.local_term_parts()
    h = np.rint(8 * (x + 8 * y)).astype(int)
    np.testing.assert_allclose(h, 8 * (x + 8 * y), atol=1e-12)
    filled = law.SRC.index((1, 1))
    b = np.kron(h[filled::9, filled::9], np.eye(9, dtype=int))
    np.testing.assert_array_equal(b @ h, h @ b)
    np.testing.assert_array_equal(b[:, filled::9], h[:, filled::9])
    # Positive control of the filled-head energy analyzer. A terminal echo
    # followed by a physical load read distinguishes the two even-sign states.
    for coin_phase, energy, extremes in [(1, 0, 1), (-1, 32, 0)]:
        head = np.zeros(99)
        for j, (coin, contrast) in enumerate(law.LABELS):
            if abs(contrast) == 2:
                head[(j + 1) * 9 + filled] = 0.5 * coin_phase**coin
        np.testing.assert_allclose(h @ head, 8 * energy * head)
        read = (law.echo() @ head).reshape(11, 9)
        hits = [law.SRC.index(pair) for pair in ((0, 2), (2, 0))]
        np.testing.assert_allclose(
            np.square(abs(read[:, hits])).sum(), extremes, atol=1e-13
        )
