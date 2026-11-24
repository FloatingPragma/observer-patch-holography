"""Independent word projections, propagation and raw-count checks for the dipole arrays."""

import json

import numpy as np
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import candidate
from model.executable_model import carrier
from . import experiment, survey


def microscopic(board, host, defect, kappa):
    law = candidate.native_law
    exchange, contrast = law.local_term_parts()
    h = exchange + kappa * contrast
    ix = [law.SRC.index((defect, host)), law.SRC.index((host, defect))]
    np.testing.assert_allclose(np.delete(h[:, ix], ix, axis=0), 0, atol=2e-14)
    source = np.diag(kappa * (6 - board.degree) * (defect**2 - host**2)).astype(float)
    for seam in board.seams:
        source[np.ix_(seam, seam)] += h[np.ix_(ix, ix)]
    return source - np.eye(board.n) * 6 * kappa


def test_projection_and_all_time_polarization_dark_subspace():
    ports = [1, 4]
    for full in (experiment.one_bridge(ports), carrier.federation(0)):
        for tag, host, defect in (("unit", 0, 1), ("vacancy", 1, 0)):
            for k in (0, 8):
                direct = microscopic(full, host, defect, k)
                np.testing.assert_allclose(
                    direct, experiment.generator(full, tag, k), atol=1e-12
                )
                if full.cells == 2:
                    f = experiment.frame(ports[0])[:, :2]
                    q = np.vstack([f, np.zeros_like(f)])
                    np.testing.assert_allclose(
                        direct @ q, q @ (q.T @ direct @ q), atol=2e-14
                    )


def test_angular_sources_have_countercirculation_with_identical_density():
    for port in range(12):
        poly = experiment.geometry()
        f = experiment.frame(port)
        wave = (f[:, 0] + 1j * f[:, 1]) / np.sqrt(2)

        def currents(s):
            return np.array(
                [np.imag(np.conj(s[i]) * s[j]) for i, j in carrier.CELL_SEAMS]
            )

        current = currents(wave)
        oriented = [
            poly[port] @ np.cross(poly[i], poly[j]) for i, j in carrier.CELL_SEAMS
        ]
        assert current @ oriented > 0.1
        np.testing.assert_allclose(-current, currents(wave.conj()))
        np.testing.assert_allclose(abs(wave) ** 2, abs(wave.conj()) ** 2)
        div = np.zeros(12)
        for j, (a, b) in zip(current, carrier.CELL_SEAMS, strict=True):
            div[a] -= j
            div[b] += j
        np.testing.assert_allclose(div, 0, atol=1e-15)


def test_dense_seam_evolution_predicts_registered_vacancy_quadratures():
    j = json.loads((experiment.HERE / "results.json").read_text())
    row = next(r for r in j["boards"] if r["prefix"] == "level0_full_vacancy_8.0")
    board = carrier.federation(0)
    src = np.zeros((board.n, 3))
    src[:12] = experiment.frame(row["source_port"])
    # Independent native microscopic projection and direct exponential action.
    h = microscopic(board, 1, 0, 8)
    u = expm_multiply(-4j * h, src)
    with np.load(experiment.HERE / "counts.npz") as data:
        for d in row["paths"]:
            sl = slice(12 * d["owner"], 12 * d["owner"] + 12)
            transfer = experiment.frame(d["facing_port"]).T @ u[sl]
            reads = data[d["key"] + "_quadratures"][3]
            measured = 2 * reads / j["shots"] - 1
            np.testing.assert_allclose(
                measured,
                [transfer.real, transfer.imag],
                atol=j["simultaneous_quadrature_error"],
            )
            counts = data[d["key"] + "_detections"]
            assert np.all(counts.sum(axis=-1) == j["shots"])
            assert np.all(counts[:-1, :, :, -1] == 0)


def test_physical_filter_is_coherent_average_including_rejected_trials():
    board = experiment.one_bridge([1, 4])
    h = microscopic(board, 0, 1, 8)
    center = (5 - np.sqrt(5)) / 2
    direct = expm(-0.25j * (h - center * np.eye(24)))
    packet = np.eye(24, dtype=complex)
    total = np.zeros_like(packet)
    for _ in range(256):
        total += packet / 256
        packet = direct @ packet
    w, v = np.linalg.eigh(h)
    avg = np.exp(-0.25j * np.outer(np.arange(256), w - center)).mean(axis=0)
    np.testing.assert_allclose(total, (v * avg) @ v.T, atol=1e-12)
    report = json.loads((experiment.HERE / "results.json").read_text())
    row = next(r for r in report["boards"] if r["prefix"] == "level0_single_unit_8.0")
    a, b = row["source_port"], row["paths"][0]["facing_port"]
    full = experiment.one_bridge([a, b])
    f = experiment.frame(a)
    avg_wave = np.vstack([f, np.zeros_like(f)]).astype(complex)
    evolve = expm(-0.25j * (microscopic(full, 0, 1, 8) - center * np.eye(24)))
    output = avg_wave.copy()
    for _ in range(255):
        avg_wave = evolve @ avg_wave
        output += avg_wave
    prob = np.sum(abs(output / 256) ** 2, axis=0)
    with np.load(experiment.HERE / "counts.npz") as data:
        n = data[row["paths"][0]["key"] + "_detections"][-1, 0, :3, -1]
        np.testing.assert_allclose(1 - prob, n / report["shots"], atol=0.0005)


def test_clock_and_frequency_summary_from_saved_counts():
    report = json.loads((experiment.HERE / "survey.json").read_text())
    full = json.loads((experiment.HERE / "results.json").read_text())
    with np.load(experiment.HERE / "flight.npz") as bank:
        for row in report["boards"]:
            assert not row["survey"]["passing_frequencies"]
            assert not row["survey"]["pass_at_all_three"]
            assert not any(row["survey"]["finite_pass_at_all_three"].values())
            if row["material"] != "unit":
                continue
            primary = next(
                r
                for r in full["boards"]
                if r["cells"] == row["cells"]
                and r["material"] == "unit"
                and r["contrast"]
                and r["bridge_setting"] == "full"
            )
            for d in primary["paths"]:
                source = bank[f"n{row['cells']}_d{d['steps']}_counts"] / report["shots"]
                for i, k in enumerate(report["launches"]):
                    for f, p in enumerate(report["thresholds"]):
                        hits = np.flatnonzero(source[:, i] >= p)
                        t = row["flight"]["readings"][str(d["steps"])][k][f][1]
                        assert t == (float(hits[0] / 2) if len(hits) else None)
    assert survey.EPS < 0.00045
