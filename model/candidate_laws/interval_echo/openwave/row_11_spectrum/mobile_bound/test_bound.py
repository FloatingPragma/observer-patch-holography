import json
from math import comb

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import coo_matrix
from model.executable_model import federation

from . import measure as ex
from . import sector as native


def direct_word_hamiltonian(n, d, kappa):
    _, addresses = native.orbits(n, d)
    words = list(addresses)
    index = {w: j for j, w in enumerate(words)}
    rows, cols, vals = [], [], []
    for j, w in enumerate(words):
        diagonal = 0.0
        for a in range(n):
            b = (a + 1) % n
            if w[a] == w[b]:
                continue
            diagonal += kappa * (w[a] - w[b]) ** 2 + 0.5
            moved = list(w)
            moved[a], moved[b] = w[b], w[a]
            rows.append(index[tuple(moved)])
            cols.append(j)
            vals.append(-0.5)
        rows.append(j)
        cols.append(j)
        vals.append(diagonal)
    return words, coo_matrix(
        (vals, (rows, cols)), shape=(len(words), len(words))
    ).tocsr()


def test_actual_board_and_translation_reduction():
    board = federation(0)
    for d in (1, 3, 5):
        p = native.prepare(d)
        n = len(p["route"])
        assert len(p["route"]) == n == 15
        assert len({i // 12 for i in p["route"]}) == 5
        assert all(
            native.edge(p["route"][a - 1], p["route"][a]) in board.seams
            for a in range(n)
        )
        words, full = direct_word_hamiltonian(n, d, 8)
        reps, addresses = native.orbits(n, d)
        assert len(words) == n * len(reps) == comb(n, d) * comb(d, (d + 1) // 2)
        sources = np.fromiter((addresses[w][0] for w in words), int)
        shifts = np.fromiter((addresses[w][1] for w in words), int)
        for mode in (0, 1, 2):
            h, center, _ = native.sector(n, d, 8, mode)
            phase = np.exp(2j * np.pi * mode * shifts / n) / np.sqrt(n)
            compression = coo_matrix(
                (phase, (np.arange(len(words)), sources)), shape=(len(words), len(reps))
            ).tocsr()
            error = full @ compression - compression @ (
                h + center * native.diags(np.ones(len(reps)))
            )
            assert max(abs(error.data), default=0.0) < 2e-13
            if d > 1:
                core = sources == 0
                cut = full[~core][:, core]
                np.testing.assert_allclose(
                    (cut.T @ cut).toarray(), np.eye(n) * (d + 1) / 4, atol=1e-15
                )


def test_complete_spectrum_and_local_recoil_flight():
    """Dense 91-dimensional blocks retain every low line in the finite physical filter."""
    record = json.loads((ex.HERE / "results.json").read_text())
    trimer = record["cases"][1]
    lines, vecs = [], []
    for m in range(ex.N):
        h, center, _ = native.sector(ex.N, 3, 8, m)
        vals, vec = eigh(h.toarray())
        lines.append(vals)
        vecs.append(vec)
    widths = np.array(lines)
    w = np.abs(np.array(vecs)[:, 0, :]) ** 2
    # Direct multiplication of the physical timing-filter transfer function on every line.
    f = ((1 + np.exp(-2j * widths / center)) / 2) ** ex.FILTER_STEPS
    p_initial = np.sum(w * abs(f) ** 2) / ex.N
    for panel in trimer["localized_flight"]["records"]:
        all_positions = np.fft.ifft(
            np.sum(w * f * np.exp(-1j * widths * panel["time"]), axis=1)
        )
        np.testing.assert_allclose(
            np.abs(all_positions) ** 2,
            np.array(panel["core_by_first_port"]) / ex.SHOTS,
            atol=0.0004,
        )
        assert abs(1 - p_initial - panel["source_refused"] / ex.SHOTS) < 0.00025
    # An explicit localized vacancy projector between the two eigenstates in full word space.
    words, _ = direct_word_hamiltonian(ex.N, 3, 8)
    _, addr = native.orbits(ex.N, 3)
    i = np.array([addr[w][0] for w in words])
    s = np.array([addr[w][1] for w in words])
    phi0 = vecs[0][i, -1] / np.sqrt(ex.N)
    for panel in trimer["recoil_after_rest"]:
        m = panel["mode"]
        phim = vecs[m][i, -1] * np.exp(2j * np.pi * m * s / ex.N) / np.sqrt(ex.N)
        kicked = phi0 * np.array([w[0] == 0 for w in words])
        assert abs(np.vdot(phim, phi0)) < 1e-13
        expected = abs(np.vdot(phim, kicked) * f[m, -1]) ** 2
        assert abs(expected - panel["p_outgoing_bound_given_source"]) < 2e-11


def test_signed_frequencies_source_accounting_and_pentamer():
    data = json.loads((ex.HERE / "results.json").read_text())
    assert len(data["cases"]) == 6
    for c in data["cases"]:
        d, kappa = c["defects"], c["kappa"]
        for bank in c["momentum_metrology"]:
            for record in bank["records"]:
                assert all(sum(x) == ex.SHOTS for x in record["counts"])
            f = ex.fit(bank["records"], c["ideal_lines"][0]["center"])
            actual = native.line(ex.N, d, kappa, bank["mode"])
            assert (
                abs(f["energy"] - actual["center"] - actual["offset"]) < f["radius_99"]
            )
        for r in c["recoil_after_rest"]:
            assert sum(r["counts"]) == ex.SHOTS and r["counts"][-1] > 5000
            assert r["matched_no_probe_outgoing_mode"] == 0
        for r in c["localized_flight"]["records"]:
            assert (
                sum(r["core_by_first_port"]) + r["other_loads"] + r["source_refused"]
                == ex.SHOTS
            )
    one, three, five = data["cases"][:3]
    assert one["measured_gaps"][0] / three["measured_gaps"][0] > 950
    assert one["measured_gaps"][0] / five["measured_gaps"][0] > 900000
    assert (
        three["localized_flight"]["records"][4]["all_outside_initial_owners"]
        > 0.7 * ex.SHOTS
    )
    assert (
        five["localized_flight"]["records"][-1]["all_outside_initial_owners"]
        > 0.42 * ex.SHOTS
    )
    assert five["highest_separated_ceiling"] == 135
    assert five["sharp_core_probability_floor"] == 0.96
    assert data["cases"][-1]["measured_gaps"][0] > 15 * five["measured_gaps"][0]
