"""Independent resolved-label calculation and blinded terminal count readback."""

import json

import numpy as np
from scipy.linalg import expm

from . import measure, readback


def test_two_read_histories_in_the_original_coin_and_sign_basis():
    n, connector, receiver = 4, 2, 3
    words = 2**n
    seam = measure.inherited.strand.candidate.native_law
    exchange, contrast = seam.local_term_parts()
    term = exchange + 8 * contrast
    head_labels = [(coin, d) for coin in (0, 1) for d in (2, -2)]
    heads = [1 + seam.LABELS.index(x) for x in head_labels]
    h = np.zeros((4 * words, 4 * words))
    wall = np.zeros(4 * words, dtype=bool)
    for word in range(words):
        bits = [(word >> i) & 1 for i in range(n)]
        ids = [tag * words + word for tag in range(4)]
        wall[ids] = bits[receiver - 1] == bits[receiver]
        for right in range(1, n):
            pair = (bits[right - 1], 1 - bits[right])
            if right == connector:
                indices = [x * 9 + seam.SRC.index(pair) for x in heads]
                h[np.ix_(ids, ids)] += term[np.ix_(indices, indices)]
            elif pair != (1, 1):
                h[ids, ids] += 32
        for tag, i in enumerate(ids):
            for dimer in range(n):
                h[i, tag * words + (word ^ (1 << dimer))] = -0.5
    flight = expm(-1j * h)
    coins = np.repeat([c for c, d in head_labels], words)
    starts = np.zeros((4, words, 2))
    for col, bond in enumerate((1, 0)):
        starts[:, measure.inherited.initial_word(n, bond), col] = 0.5
    psi = flight @ starts.reshape(-1, 2)
    joint = []
    for a in (0, 1):
        between = flight @ ((coins == a)[:, None] * psi)
        for b in (0, 1):
            finish = flight @ (flight @ ((coins == b)[:, None] * between))
            joint.append(
                [np.sum(abs(finish[mask]) ** 2, axis=0) for mask in (wall, ~wall)]
            )
    output = list(measure.histories(8, n, connector, receiver, (1, 0), (2, 4)))
    np.testing.assert_allclose(output[-1][1]["spaced"], joint, atol=4e-13, rtol=0)


def test_replay_terminal_verdicts_using_counts_alone():
    data = json.loads((measure.HERE / "records.json").read_text())
    banks = [
        {k: v for k, v in b.items() if k != "probabilities"} for b in data["banks"]
    ]
    read = readback.decode(banks, data["clocks"], data["shots"], data["sources"])
    np.testing.assert_allclose(
        read["clock_intervals"], data["readings"]["clock_intervals"], atol=1e-10, rtol=0
    )
    for measured, archived in zip(
        read["comparisons"], data["readings"]["comparisons"], strict=True
    ):
        for key, values in measured["source_responses"].items():
            np.testing.assert_allclose(
                values, archived["source_responses"][key], atol=1e-10, rtol=0
            )
    assert all(
        lo < t < hi
        for t, (lo, hi) in zip(measure.HANDS, read["clock_intervals"], strict=True)
    )
    for row in read["comparisons"]:
        bands = row["source_responses"]
        if row["kappa"] == 0:
            assert all(lo <= 0 <= hi for lo, hi in bands.values())
        if row["appointment"] == 2:
            assert row["spaced_minus_once"][0] <= 0 <= row["spaced_minus_once"][1]
    assert len(banks) == 2 * 3 * 4 * 2


def test_pooled_rapid_reads_keep_half_the_connector():
    (attached, cut), _ = measure.generators(8, 4, 2, 3)
    a = attached.toarray()
    b = cut.toarray()
    h = (a + b) / 2
    t = 0.01
    u = expm(h * t)
    rng = np.random.default_rng(701)
    v = rng.normal(size=len(a)) + 1j * rng.normal(size=len(a))
    rho = np.outer(v, v.conj()) / np.vdot(v, v).real
    limit = u @ rho @ u.conj().T
    errs = []
    for steps in (1, 4, 16):
        x, y = expm(a * t / steps), expm(b * t / steps)
        r = rho
        for _ in range(steps):
            r = (x @ r @ x.conj().T + y @ r @ y.conj().T) / 2
        error = np.sum(abs(np.linalg.eigvalsh(r - limit))) / 2
        assert error <= 64 * t**2 / steps
        errs.append(error)
    assert errs[2] < errs[1] < errs[0]


def test_native_magnitude_and_changed_law_against_other_propagators():
    data = json.loads((measure.HERE / "records.json").read_text())
    for bank in data["banks"]:
        if bank["menu"] != "magnitude":
            continue
        t, j = bank["appointment"], bank["source"]
        if bank["kappa"] == 0:
            p = np.sin(t) ** 2 / 2
        else:
            # The dual kernel acts on signed wall density, including the common far-end spectator.
            kernel = measure.inherited.strand.walls(len(data["route"]), 8, t)[7]
            p = (1 - kernel.sum()) / 2 + kernel[j] + kernel[-1]
        np.testing.assert_allclose(bank["probabilities"][0], p, atol=3e-11, rtol=0)
