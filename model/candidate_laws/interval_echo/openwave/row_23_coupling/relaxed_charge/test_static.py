"""Independent microscopic, spectral, profile and count checks on the signed work."""

import gzip
import json
from functools import cache

import numpy as np
from scipy.linalg import eigvalsh
from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation

from ..free_bath import native
from . import measure, readback


@cache
def archive():
    return json.loads(gzip.decompress(measure.OUTPUT.read_bytes()))


def direct_words(board, k, pins, shutter=None):
    """Act with each native seam on an actual ternary word and retain the two pins."""
    cols = candidate.native_law.SRC
    h0, h1 = candidate.native_law.local_term_parts()
    seam = h0 + k * h1
    assert not np.any(abs(seam[9:, :9]) > 1e-13)
    ports = [v for v in range(board.n) if v not in pins]
    ix = dict(zip(ports, range(len(ports))))
    base = np.ones(board.n, dtype=int)
    for v, q in pins.items():
        base[v] = q
    edges = [
        (a, b)
        for a in range(board.n)
        for b in board.neighbours[a]
        if a < b and (shutter is None or shutter[a] == shutter[b])
    ]
    ai, bi = np.array(edges).T
    diagonal = np.array([seam[cols.index((a, b)), cols.index((a, b))] for a, b in cols])

    def energy(w):
        return ((6 - board.degree) * k) @ (w**2) + diagonal[w[ai] * 3 + w[bi]].sum()

    offset = energy(base)
    out = np.zeros((len(ports), len(ports)))
    for i, v in enumerate(ports):
        word = base.copy()
        word[v] = 0
        out[i, i] = energy(word) - offset
        for a, b in edges:
            if a not in ix or b not in ix or word[a] == word[b]:
                continue
            source = cols.index((word[a], word[b]))
            for target in np.flatnonzero(abs(seam[:9, source]) > 1e-14):
                if target == source:
                    continue
                changed = word.copy()
                changed[a], changed[b] = cols[target]
                vacant = [v for v in (a, b) if changed[v] == 0]
                assert len(vacant) == 1
                out[ix[vacant[0]], i] += seam[target, source]
    return np.array(ports), out


def test_two_charges_and_shutter_against_full_words():
    board = federation(0)
    row = archive()["boards"]["20,8"]
    shell = next(r for r in row["shell_banks"] if r["distance"] == 4)
    home, probe = row["source"], shell["port"]
    sides = native.distances(board, home) <= native.distances(board, probe)
    pins = {home: 0, probe: 2}
    free, whole = direct_words(board, 8, pins)
    projected, ix = native.vacancy_operator(board, 8, pins)
    np.testing.assert_array_equal(ix, free)
    np.testing.assert_allclose(whole, projected.toarray(), atol=6e-11, rtol=0)
    _, shut = direct_words(board, 8, pins, sides)
    np.testing.assert_array_equal(shut[np.ix_(sides[free], ~sides[free])], 0)
    for side, bank in zip((sides[free], ~sides[free]), shell["shutter"]["-1,1"]):
        ground = eigvalsh(shut[np.ix_(side, side)], subset_by_index=[0, 0])[0]
        assert (
            abs(ground - readback.clock_energy(bank["clock"])) < readback.ENERGY_ERROR
        )


def test_complete_ternary_menus_and_independent_clock_lifts():
    data = archive()
    for key, row in data["boards"].items():
        k = key.split(",")[1]
        cal = data["calibration"][k]
        for r in row["shell_banks"]:
            signs = {f"{a},{b}" for a in (-1, 1) for b in (-1, 0, 1)}
            assert set(r["intact"]) == signs
            read = {}
            for key, bank in r["intact"].items():
                assert sum(bank["address"]) == readback.SHOTS
                q = np.array(bank["clock"]) / readback.SHOTS
                phase = np.angle(2 * q[:, 0] - 1 + 1j * (2 * q[:, 1] - 1))
                times = np.array(data["clock_times"])
                # Independently align each observed phase to the microscopic line cluster.
                e0 = bank["energy_lines"][0]
                for t, a in zip(times, phase):
                    angle = a + 2 * np.pi * np.rint((t * e0 - a) / (2 * np.pi))
                    assert abs(angle - t * e0) < readback.ANGLE_ERROR
                read[key] = angle / times[-1]
            result = r["readings"]
            direct = cal["direct_four_sign_bins"][int(r["distance"] == 1)] / 4
            measured = (
                read["1,1"] - read["1,-1"] + read["-1,-1"] - read["-1,1"]
            ) / 4 + direct
            assert abs(measured - result["relaxed"]) < 1e-10
            for pair, interval in result["local_contact_intervals"].items():
                if interval is not None:
                    c = cal["gain"][pair.split(",")[1]]
                    assert interval[0] <= c <= interval[1]


def test_stock_factorization_and_coupling_identifiability():
    for row in archive()["boards"].values():
        for r in row["shell_banks"]:
            result = r["readings"]
            if r["distance"] >= 3:
                assert abs(result["shutter_two_vacancies"]) < 2 * readback.ENERGY_ERROR
            for bank in r["intact"].values():
                assert bank["dimension"] == 12 * row["owners"] - 2
    # The Perron estimator also detects a changed local coefficient with unchanged hopping.
    from scipy.sparse import csr_matrix

    g = csr_matrix([[2.0, -1.0, 0.0], [-1.0, 2.0, -1.0], [0.0, -1.0, 2.0]])
    aperture = np.array([False, True, False])
    for c in (0.5, -3.0):
        h = g.toarray() + np.diag(aperture * c)
        e, v = np.linalg.eigh(h)
        rng = np.random.Generator(np.random.PCG64(2405))
        data = dict(
            multiplicity=1,
            clock=measure.clock_bank(e[:1], rng),
            address=rng.multinomial(readback.SHOTS, v[:, 0] ** 2).tolist(),
        )
        bounds = readback.contact_interval(g, aperture, data)
        assert bounds[0] <= c <= bounds[1] and bounds[1] - bounds[0] < 0.04


def test_two_pins_and_mobile_vacancy_have_a_unique_local_address_read():
    from itertools import combinations, product

    for level in (0, 1):
        board = federation(level)
        code = native.location_instrument(board)
        for tree in code["distinct_local_trees"]:
            for ports in combinations(range(12), 3):
                for marks in product((0, 2), repeat=3):
                    word = np.ones(12, dtype=int)
                    word[list(ports)] = marks
                    recovered = np.zeros(12, dtype=int)
                    for a, b in tree:
                        recovered[b] = recovered[a] - (word[a] - word[b])
                    # Nine ordinary ports calibrate the common additive constant.
                    assert np.array_equal(
                        word, recovered - int(np.median(recovered)) + 1
                    )
