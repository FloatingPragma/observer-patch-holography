"""Cross-check every receiver bank using raw incidences and native seam columns."""

import json
from itertools import permutations

import numpy as np
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import clock
from .. import experiment as dipole
from . import measure, onset


def direct_local(words, kappa):
    cell = dipole.carrier.cell()
    labels = {word: i for i, word in enumerate(words)}
    h0, h1 = clock.law.local_term_parts()
    seam = h0 + kappa * h1
    h = np.zeros((len(words), len(words)))
    for i, word in enumerate(words):
        loads = np.zeros(12, dtype=int)
        loads[word[0]] = 2
        if len(word) == 2:
            loads[word[1]] = 1
        h[i, i] = kappa * ((6 - cell.degree) @ (loads**2))
        for a, b in cell.seams:
            u, v = loads[[a, b]]
            col = clock.law.SRC.index((u, v))
            h[i, i] += seam[col, col]
            if u != v:
                swap = {a: b, b: a}
                row = labels[tuple(swap.get(x, x) for x in word)]
                h[row, i] += seam[clock.law.SRC.index((v, u)), col]
    np.testing.assert_array_equal(h, h.T)
    return h - h[0, 0] * np.eye(len(h))


def same_sample(stored, probabilities, shots):
    p = np.asarray(probabilities)
    assert min(p) >= -1e-11
    np.testing.assert_allclose(sum(p), 1, atol=1e-11, rtol=0)
    assert sum(stored) == shots
    mean = p * shots
    assert np.all(abs(stored - mean) < 9 * np.sqrt(np.maximum(mean, 0)) + 25)


def test_every_face_and_probe_address():
    result = json.loads((measure.HERE / "records.json").read_text())
    shots = result["shots"]
    points = dipole.geometry()
    board = dipole.carrier.federation(1)
    h = dipole.candidate.generator(board).tocsr()
    h.setdiag(h.diagonal() - 6 * clock.KAPPA)
    initial = np.zeros((board.n, 3))
    initial[:12] = dipole.frame(1)
    propagated = expm_multiply(-64j * h, initial)[12 * 18 : 12 * 19]
    x, y, z = propagated.T
    sources = dict(
        zip(
            measure.NAMES,
            (
                x,
                y,
                (x + y) / np.sqrt(2),
                (x - y) / np.sqrt(2),
                (x + 1j * y) / np.sqrt(2),
                (x - 1j * y) / np.sqrt(2),
                z,
                x * 0,
            ),
            strict=True,
        )
    )
    owner = dipole.carrier.cell()
    for name, psi in sources.items():
        occupancy = np.vdot(psi, psi).real
        adm = psi.copy()
        adm[result["probe_ports"]] = 0
        p_admit = np.vdot(adm, adm).real
        decoded = np.zeros((2, 3))
        for face_index, (face, counts) in enumerate(
            zip(result["face_order"], result["face_counts"][name], strict=True)
        ):
            a, b, c = face
            assert (
                np.cross(points[b] - points[a], points[c] - points[a]) @ points[a] > 0
            )
            j = np.zeros((12, 12), complex)
            for i, k in ((a, b), (b, c), (c, a)):
                assert (min(i, k), max(i, k)) in owner.seams
                j[k, i], j[i, k] = 0.5j, -0.5j
            # The nonzero eigenvalues are ±sqrt(3)/2. Each sign is a polynomial.
            normalized = (2 / np.sqrt(3)) * j
            assert (
                np.max(abs(normalized @ normalized @ normalized - normalized)) < 1e-14
            )
            bins = [
                np.vdot(
                    psi, ((normalized @ normalized + sign * normalized) / 2) @ psi
                ).real
                for sign in (1, -1)
            ]
            same_sample(counts, [*bins, occupancy - sum(bins), 1 - occupancy], shots)
            matched = result["matched_face_counts"][name][face_index]
            mbins = [
                np.vdot(
                    adm, ((normalized @ normalized + sign * normalized) / 2) @ adm
                ).real
                for sign in (1, -1)
            ]
            same_sample(
                matched,
                [*mbins, p_admit - sum(mbins), 1 - occupancy, occupancy - p_admit],
                shots,
            )
            normal = np.cross(points[b] - points[a], points[c] - points[a])
            normal /= np.linalg.norm(normal)
            for k, item in enumerate((counts, matched)):
                decoded[k] += (
                    (item[0] - item[1]) * np.sqrt(3) * normal / (2 * shots * 20)
                )
        for k, field in enumerate(("axial_reading", "matched_axial_reading")):
            np.testing.assert_allclose(
                result[field][name], decoded[k], atol=5.1e-11, rtol=0
            )
    words = list(permutations(range(12), 2))
    a, b = result["probe_ports"]
    direction = (points[b] - points[a]) / np.linalg.norm(points[b] - points[a])
    for law, kappa in (("native", clock.KAPPA), ("zero_contrast", 0)):
        h2 = direct_local(words, kappa)
        h1 = direct_local([(p,) for p in range(12)], kappa)
        for window in measure.TIMES:
            p1, p2 = expm((-1j * window) * h1), expm((-1j * window) * h2)
            for name, psi in sources.items():
                inside = np.vdot(psi, psi).real
                signal = psi.copy()
                signal[[a, b]] = 0
                admit = np.vdot(signal, signal).real
                rec = result["records"][law][str(window)][name]
                changes = {}
                for mode, phase in (("plus", 1j), ("minus", -1j), ("real", 1)):
                    start = np.zeros(12, complex)
                    start[[a, b]] = np.array([1, phase]) / np.sqrt(2)
                    single = abs(p1 @ start) ** 2
                    loaded = p2 @ np.array(
                        [start[probe] * signal[arrival] for probe, arrival in words]
                    )
                    double = np.zeros(12)
                    np.add.at(double, [p for p, _ in words], abs(loaded) ** 2)
                    for tag, local in (("contact", double), ("stored", admit * single)):
                        same_sample(
                            rec["counts"][mode][tag],
                            [inside - admit, *((1 - inside) * single), *local],
                            shots,
                        )
                    changes[mode] = (double - admit * single) @ points
                odd = (changes["plus"] - changes["minus"]) / 2
                assert (
                    abs(rec["odd_parallel"] - odd @ direction)
                    < result["contrast_half_width"]
                )
                if name != "blank" and law == "native":
                    assert rec["odd_parallel"] + result["contrast_half_width"] < 0


def test_initial_current_and_acceleration():
    onset.audit()
