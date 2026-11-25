from itertools import product
import json

import networkx as nx
import numpy as np

from model.candidate_laws.interval_echo import candidate
from model.executable_model.carrier import federation
from .measure import HERE, TIMES, SHOTS, ENERGY_GRID, decode, energies, summary
from .mixed_clock import (
    BINS,
    GRID,
    LAUNCHES,
    PERIOD,
    REPEATS,
    failure_bound,
    ground_reads,
    native_matrix,
)


def receipt():
    return json.loads((HERE / "records.json").read_text())


def test_routes_sharp_standards_and_full_census():
    board = federation(0)
    graph = nx.Graph(board.seams)
    for entry in receipt()["rings"]:
        r = entry["layout"]
        pairs = r["pairs"]
        n = len(pairs)
        ports = [p for edge in pairs for p in edge] + r["pins"] + r["standards"]
        assert len(set(ports)) == len(ports)
        assert (len(set(p // 12 for row in pairs for p in row)), n) in ((1, 3), (5, 11))
        for k, (a, b) in enumerate(pairs):
            assert graph.has_edge(a, b) and graph.has_edge(b, pairs[(k + 1) % n][0])
            assert board.degree[a] == board.degree[b]
        for k, (pin, standard) in enumerate(
            zip(r["pins"], r["standards"], strict=True)
        ):
            assert graph.has_edge(pin, standard) and graph.has_edge(
                pin, pairs[r["marks"][k]][0]
            )
            meter = entry["instruments"]["meters"][k]
            for q, counts in zip((-1, 0, 1), meter["terminal_words"], strict=True):
                assert sum(counts) == 16384
                assert all(
                    c == 0 or sum(w) == q + 2
                    for c, w in zip(counts, candidate.native_law.SRC, strict=True)
                )
        for code, path in r["routes"].items():
            a, b = map(int, code)
            assert path[0] == r["pins"][a] and path[-1] == r["pins"][b]
            assert all(graph.has_edge(x, y) for x, y in zip(path, path[1:]))
            assert len(path) - 1 == nx.shortest_path_length(graph, path[0], path[-1])
            good, cut = entry["instruments"]["ranging"][code]
            assert good[-1] == sum(good) == 4096 and sum(cut) == 0


def test_word_generator_and_certified_energies():
    small = receipt()["rings"][0]
    for bank in small["banks"]:
        for panel in bank["panels"]:
            args = (small["layout"], panel["q"], bank["cut"], bank["kinetic"])
            h = native_matrix(*args)
            eigen = np.linalg.eigvalsh(h)
            e, err = energies(*args)
            assert abs(eigen[0] - e) < err + 5e-12
            assert (
                abs(panel["energy_tick"] / ENERGY_GRID - eigen[0])
                < 0.5 / ENERGY_GRID + 5e-12
            )
            assert eigen[-1] < 768


def test_outer_clock_signed_pairs_loop_cut_and_three_source_residual():
    eps = np.sqrt(2 * np.log(2 * 30000 / 0.01) / SHOTS)
    angle = np.arcsin(np.sqrt(2) * eps)
    tol = angle / TIMES[-1] + 0.5 / ENERGY_GRID
    assert 3 * angle < np.pi
    for ring in receipt()["rings"]:
        n = len(ring["layout"]["pairs"])
        for bank in ring["banks"]:
            for panel in bank["panels"]:
                d = np.array(panel["populations"])
                assert d.shape == (2, len(TIMES)) and d.min() >= 0 and d.max() <= SHOTS
                e = panel["energy_tick"] / ENERGY_GRID
                z = 2 * d / SHOTS - 1
                for k, y in enumerate((np.cos(TIMES * e), np.sin(TIMES * e))):
                    assert np.max(abs(z[k] - y)) < eps
                assert abs(decode(panel) - e) < tol
        a, t = summary(ring["banks"][0])
        for key in ("01", "12", "02"):
            assert a[key]["odd"] > tol and a[key]["even"] > 4 * tol
        if n == 11:
            assert a["03"]["odd"] < -tol
        assert max(t.values()) > 30 + 8 * tol
        first, _ = summary(ring["banks"][1])
        second, _ = summary(ring["banks"][2])
        assert first["01"]["odd"] < -15 and second["01"]["odd"] > 15
        disabled, _ = summary(ring["banks"][3])
        assert all(abs(v["odd"]) < tol for v in disabled.values())
        # Fix strength on the 0,1 terminal pair and use independent material range.
        paths = ring["layout"]["routes"]
        for held in ("02",) if n == 3 else ("02", "12", "03"):
            ratio = (len(paths["01"]) - 1) / (len(paths[held]) - 1)
            assert abs(a[held]["odd"] - ratio * a["01"]["odd"]) > (1 + ratio) * tol


def test_blind_clock_all_settings_and_no_eigenstate_selection():
    with np.load(HERE / "clock.npz", allow_pickle=False) as f:
        raw = f["exits"]
    assert raw.shape == (27, LAUNCHES, REPEATS)
    assert raw.min() >= 0 and raw.max() < BINS
    measured = ground_reads(raw)
    d = 8 * PERIOD / BINS + 0.5 / GRID
    assert failure_bound() < 0.00002
    layout = receipt()["rings"][0]["layout"]
    expected = np.array(
        [
            np.linalg.eigvalsh(native_matrix(layout, q))[0]
            for q in product((-1, 0, 1), repeat=3)
        ]
    ).reshape((3, 3, 3))
    assert np.max(abs(measured - expected)) < d
    for e in (measured[:, :, 1], measured[:, 1, :], measured[1, :, :]):
        contrast = (e[0, 0] + e[2, 2] - e[0, 2] - e[2, 0]) / 4
        assert contrast > d
    for s in (0, 2):
        triple = (
            measured[s, s, s]
            - measured[s, s, 1]
            - measured[s, 1, s]
            - measured[1, s, s]
            + measured[s, 1, 1]
            + measured[1, s, 1]
            + measured[1, 1, s]
            - measured[1, 1, 1]
        )
        assert triple - 8 * d > 31
    # The clock is detecting several populated eigenlines on unselected launches.
    med = np.median(raw[13], axis=1) * PERIOD / BINS
    assert (med.max() - med.min()) > 30


def test_clock_probabilities_match_coherent_pointer_evolution():
    from scipy.linalg import expm
    from .mixed_clock import clock_distribution

    small = receipt()["rings"][0]["layout"]
    h = native_matrix(small, (1, 0, -1))
    # The forward circular first moment of a length-M Fourier clock equals a
    # physical one-step coherence plus its one wrap-around coherence.
    poles = np.linalg.eigvalsh(h)
    phase = np.exp(2j * np.pi * np.arange(BINS) / BINS)
    mean = 0j
    for v in poles:
        cumulative = clock_distribution(int(round(v * GRID)))
        mean += np.dot(np.diff(cumulative, prepend=0), phase) / len(poles)
    base = 2j * np.pi / PERIOD
    expected = np.trace(
        (1 - 1 / BINS) * expm(base * h) + expm(-(BINS - 1) * base * h) / BINS
    ) / len(poles)
    assert abs(mean - expected) < 4e-10
