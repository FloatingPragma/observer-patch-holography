"""Independent charge, eigenvalue and interface checks on the emitted packet."""

from hashlib import sha256
import json
import numpy as np
from scipy.linalg import eigh
from . import strand_exterior as ex, weak_strand as ws
from .measure import HERE


def record():
    return json.loads((HERE / "strand_exterior.json").read_text())


def test_minmax_and_local_port_census():
    q = record()
    assert (
        q["protocol_sha256"]
        == sha256((HERE / "STRAND_EXTERIOR.md").read_bytes()).hexdigest()
    )
    for s in q["blank_standards"]:
        assert s["low_count"] == s["cells"] and s["upper_low"] == 0.25
        assert 1.38 < s["next_lower"] < 1.39
        assert (
            s["one_pair_potential_floor"] == 80
            and s["one_vacancy_potential_floor"] == 32
        )
        assert (
            sum(item["cells"] for item in s["isolated_affordable_census"]) == s["cells"]
        )
        for intervals in q["source_intervals"].values():
            for lo, hi in intervals:
                assert 0.25 < lo < hi < s["next_lower"]
                for c in s["isolated_affordable_census"]:
                    assert all(not lo <= v <= hi for v in c["gaps"])
    # An independently projected, exact small board confirms the min-max
    # low subspace. The unrestricted Hamiltonian need not preserve this subspace.
    b = ex.federation(0)
    ports = np.arange(b.n)
    basis = np.eye(b.cells)[ports // 12] / np.sqrt(12)
    a, z = np.asarray(b.seams).T
    projection = basis[a] - basis[z]
    spectrum = eigh(projection.T @ projection / 2, eigvals_only=True)
    assert -1e-12 < spectrum[0] < 1e-12 and spectrum[-1] <= 0.25 + 1e-12
    source = q["blank_standards"][0]["smallest_board_direct"]
    np.testing.assert_allclose(source[2], (5 - np.sqrt(5)) / 2, atol=2e-12)


def test_full_ternary_contrast_and_exchange_on_the_actual_shutter():
    b, d, heads, _ = ws.ring_layout(6, True, False)
    path, cell, internal = ex.geometry(b, d, heads)
    words, opened = ex.mobile_generator(b, d, heads, path, cell, internal, "open")
    _, cut = ex.mobile_generator(b, d, heads, path, cell, internal, "link_cut")
    delta = (opened - cut).tocsr()
    lookup = {frozenset(w): i for i, w in enumerate(words)}
    m = len(words)
    a, z = path[:2]
    for j, w in enumerate(words):
        expected = np.zeros(m)
        if (a in w) != (z in w):
            expected[j] = 8.5
            expected[lookup[frozenset(w) ^ {a, z}]] = -0.5
        np.testing.assert_allclose(delta[:m, j].toarray().ravel(), expected, atol=1e-13)
    assert all(p not in set(sum(d, ())) for p in cell)


def test_all_lines_preserve_tag_probability_and_the_frozen_sign():
    for n, s in record()["panels"].items():
        ix = s["receiving_affordable"]
        for mount, rows in s["settings"].items():
            for page in rows:
                for name, arms in page.items():
                    expected = s["retained_probabilities"][name]
                    for scored in arms.values():
                        np.testing.assert_allclose(
                            sum(sum(w) for w in scored.values()), expected, atol=2e-10
                        )
            last = rows[-1]["native_pulse"]
            diff = sum(np.take(last["retain"]["1"], ix)) - sum(
                np.take(last["reset"]["1"], ix)
            )
            minimum, maximum = s["all_enrolled_state_contrast_bounds"][mount][-1]
            prob = s["retained_probabilities"]["native_pulse"]
            assert (minimum - 1e-11) * prob <= diff <= (maximum + 1e-11) * prob
            if mount == "link_cut":
                assert abs(diff) < 1e-11
            else:
                assert maximum < -0.0001 and -0.0003 < diff < -0.0002
        assert s["shutter_work"]["open"]["native_pulse"]["keep"] > 0.7
        finite = ex.finite_soft_contrast(s, s["settings"]["open"][-1]["native_pulse"])
        assert -0.0003 < finite < -0.0002


def test_replay_compares_degenerate_energies_as_physical_outcomes():
    saved, observed = record(), record()
    arm = observed["panels"]["8"]["settings"]["open"][-1]["native_pulse"]["retain"]
    assert abs(arm["1"][1] - arm["1"][2]) > 1e-5
    arm["1"][1], arm["1"][2] = arm["1"][2], arm["1"][1]
    ex.compare_exterior(saved, observed)
