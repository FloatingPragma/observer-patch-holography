"""Independent native-word, energy-line, encounter and detector audits."""

import json
from itertools import product
import numpy as np
from scipy.linalg import expm, eigh
from scipy.stats import beta
from model.candidate_laws.contrast_memory_spherical_v1 import law
from . import experiment as ex, converter, calorimeter


def test_full_native_conservation_and_local_extinction():
    energy, contrast = law.local_term_parts()
    echo = law.echo()
    for k in (0, 8):
        h = energy + k * contrast
        for amount in (0, 1, 2):
            counts = np.tile(
                [int(a == amount) + int(b == amount) for a, b in law.SRC], 11
            )
            np.testing.assert_allclose(
                h * (counts[None, :] - counts[:, None]), 0, atol=1e-12
            )
        for a, b in ((0, 2), (2, 0)):
            y = echo[:, law.SRC.index((a, b))]
            assert {i % 9 for i in np.flatnonzero(abs(y) > 1e-10)} == {
                law.SRC.index((1, 1))
            }
            assert {
                law.LABELS[i // 9 - 1][1] for i in np.flatnonzero(abs(y) > 1e-10)
            } == {a - b}
            np.testing.assert_allclose(h @ y, 0, atol=1e-12)


def word_matrix(n, k):
    board, ports, wings = ex.geometry()
    sites = ports + wings[0][: 2 * n] + wings[1][: 2 * n]
    pairs = ex.PAIRS
    states = []
    for a, b in pairs:
        for row in product((0, 1), repeat=2 * n):
            c = [1] * 4
            c[a] = 2
            c[b] = 0
            states.append(tuple(c + [v for bit in row for v in (2 * bit, 2 - 2 * bit)]))
    index = {w: i for i, w in enumerate(states)}
    edges = [(0, 1), (1, 2), (2, 3)]
    for start, inner in ((4, 1), (4 + 2 * n, 2)):
        edges.append((inner, start))
        edges.extend((j, j + 1) for j in range(start, start + 2 * n - 1))
    full = law.seam_term(k)
    H = np.zeros((len(states),) * 2)
    for col, w in enumerate(states):
        H[col, col] = k * sum((6 - board.degree[p]) * a * a for p, a in zip(sites, w))
        for a, b in edges:
            src = law.SRC.index((w[a], w[b]))
            for dst in np.flatnonzero(full[:, src]):
                v = list(w)
                v[a], v[b] = law.SRC[dst]
                row = index.get(tuple(v))
                if row is not None:
                    H[row, col] += full[dst, src]
    return H


def test_complete_independent_word_generator():
    for k in (0, 8):
        h = ex.basis(k, 2)[0].toarray()
        direct = word_matrix(2, k)
        np.testing.assert_allclose(
            direct - h, np.eye(len(h)) * (direct[0, 0] - h[0, 0]), atol=1e-11
        )


def test_first_encounter_total_from_dense_microscopic_propagator():
    case = json.loads((ex.HERE / "records.json").read_text())["cases"][0]
    assert (case["contrast"], case["dimers_per_wing"], case["launch_sign"]) == (
        8,
        2,
        -1,
    )
    H = word_matrix(2, 8)
    # With both incoming matter marks far from the center, its diagonal block
    # contains the free-wing ground product and a constant central cost.
    i = ex.PAIRS.index((3, 0))
    d = 16
    state = np.zeros(len(H), complex)
    state[i * d : (i + 1) * d] = eigh(H[i * d : (i + 1) * d, i * d : (i + 1) * d])[1][
        :, 0
    ]
    U = expm(-0.25j * (H - H[0, 0] * np.eye(len(H))))
    outcomes = []
    for _ in range(16):
        state = U @ state
        for addr in ex.READY:
            sl = slice(addr * d, (addr + 1) * d)
            outcomes.append(float(np.vdot(state[sl], state[sl]).real))
            state[sl] = 0
    np.testing.assert_allclose(
        [*outcomes, np.vdot(state, state).real],
        case["opportunity"]["probabilities"],
        atol=1e-10,
    )


def test_remote_counts_and_causal_cut():
    results = json.loads((ex.HERE / "records.json").read_text())
    # Binomial intervals conditional on the observed number of heralds. Union
    # coverage includes every arm, sign, wing size, time and parity bin.
    tail = 0.01 / (2 * 8 * 4 * 5 * 4)
    for c in results["cases"]:
        for r in c["detectors"]:
            assert sum(r["counts"]) == results["shots"]
        for t in ex.TIMES:
            cut = [
                r for r in c["detectors"] if r["arm"].endswith("_cut") and r["t"] == t
            ]
            np.testing.assert_allclose(
                cut[0]["probabilities"], cut[1]["probabilities"], atol=2e-11
            )
        if c["contrast"] != 8 or c["dimers_per_wing"] != 3:
            continue
        arms = [
            r
            for r in c["detectors"]
            if r["t"] == 32 and r["arm"] in ("converted", "parked")
        ]
        bands = []
        for r in arms:
            k = r["counts"][3]
            n = sum(r["counts"][:4])
            bands.append((beta.ppf(tail, k, n - k + 1), beta.isf(tail, k + 1, n - k)))
        assert bands[0][0] > bands[1][1]
        assert 0.64 < sum(c["wing_energy_above_ground_given_hit"]) < 0.67
    assert all(
        case["extinction_and_either_remote_excitation_zero"]["counts"][1] == 0
        for case in results["uniform"]
    )


def test_autonomous_converter_vacuum_and_resolved_calorimetry():
    result = json.loads((ex.HERE / "converter.json").read_text())
    assert converter.primitive() == 0
    for c in result["studies"]:
        for r in c["records"]:
            assert r["conserved_energy_residual"] < 1e-7
            if c["coupling"] == 0:
                assert all(
                    x["energy_outcomes"]["counts"][-1] == ex.SHOTS
                    for x in r["launches"]
                )
    c = json.loads((ex.HERE / "calorimeter.json").read_text())
    joint = np.array(
        c["constant_converter"]["full_spectrum"]["probabilities"][:-1]
    ).reshape(8, 8)
    assert (joint[2:, :2].sum() + joint[:2, 2:].sum()) / joint.sum() > 0.999
    for r in c["neutral_ground"]["paired_detectors"]:
        assert r["counts"][3] == 0
    for r in c["constant_converter"]["paired_detectors"]:
        assert sum(r["counts"]) == 2**30


def test_finite_receiver_averages_a_local_operator():
    h = ex.wing(8, 3, 1)
    e, v = eigh(h)
    eff, settings = calorimeter.terminal_sensors(h)
    j = 1
    delt = settings[j]["frequency"]
    theta = 2 * np.pi / (calorimeter.BINS * (e[1] - e[0]))
    N = calorimeter.BINS
    x = np.kron(np.eye(4), np.array([[0, 1], [1, 0]]))
    # Evaluate the time sum explicitly on a calibrated upper eigenstate.
    phase = e - e[2] + delt
    local = v.T @ x @ v[:, 2]
    t = np.arange(N) * theta
    amp = np.array(
        [local[i] * np.exp(1j * de * t).mean() for i, de in enumerate(phase)]
    )
    np.testing.assert_allclose(
        v[:, 2] @ eff[j] @ v[:, 2], np.vdot(amp, amp), atol=2e-10
    )
