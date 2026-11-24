"""Independent microscopic receiver, stationary support and orbital calculations."""

import json
from itertools import product
from pathlib import Path

import numpy as np
import pytest
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.interval_echo.openwave.quark_winding import registers
from model.candidate_laws.interval_echo.openwave.quark_winding.belt_walls import (
    escape,
    external_field,
    localized,
    orbital_singlet,
)
from model.candidate_laws.interval_echo.openwave.quark_winding.belt_walls.test_walls import (
    compare,
)

HERE = Path(__file__).parent


def stored(name):
    return json.loads((HERE / (name + ".json")).read_text())


def test_remote_against_independent_native_seams():
    n = 4
    kappa = 8
    records = list(product((0, 2), repeat=n))
    ix = {w: i for i, w in enumerate(records)}
    native = sum(
        k * part for k, part in zip((1, kappa), law.local_term_parts(), strict=True)
    )
    for sample, reference in ((1, 0), (0, 1), (1, 1)):
        h = np.zeros((len(records), len(records)))
        for i, w in enumerate(records):
            tips = [(x, 2 - x) for x in w]
            contacts = [(sample, tips[0][0]), (reference, tips[0][1])]
            contacts.extend((a[1], b[0]) for a, b in zip(tips, tips[1:]))
            for a, b in contacts + tips:
                at = law.SRC.index((a, b))
                h[i, i] += native[at, at]
            for j, (a, b) in enumerate(tips):
                h[ix[w[:j] + (b,) + w[j + 1 :]], i] += native[
                    law.SRC.index((b, a)), law.SRC.index((a, b))
                ]
        u = expm(-1j * 5.625 * h)
        start = [ix[(0, 2, 0, 2)], ix[(2, 0, 2, 0)]]
        p = np.sum(abs(u[:, start]) ** 2, axis=1) / 2
        measured = sum((w[-2] - 1) * v for w, v in zip(records, p, strict=True))
        assert external_field.read(n, sample - reference, kappa, [5.625])[
            0
        ] == pytest.approx(measured, abs=2e-10)


def test_remote_counts_and_exact_apertures():
    raw = stored("external_field")
    replay = json.loads(json.dumps(external_field.run()))
    # At zero expectation, binomial samplers can select complementary internal
    # branches when a floating probability crosses 1/2 by rounding.
    for label in ("registered", "exchange_only"):
        for size, panel in replay[label].items():
            reference = raw[label][size]
            assert abs(panel["signal"]["0"]) < 1e-11
            assert abs(reference["signal"]["0"]) < 1e-11
            assert abs(panel["broken"]) < 1e-11
            assert abs(reference["broken"]) < 1e-11
            observed = panel["calibration"]["0"]
            saved = reference["calibration"]["0"]
            assert observed["plus"] + observed["minus"] == external_field.SHOTS
            assert observed["mean"] == pytest.approx(
                2 * observed["plus"] / external_field.SHOTS - 1, abs=5.1e-12, rel=0
            )
            assert observed["plus"] in (saved["plus"], saved["minus"])
            panel["calibration"]["0"] = saved
            assert panel["broken_plus"] in (
                reference["broken_plus"], external_field.SHOTS - reference["broken_plus"]
            )
            panel["broken_plus"] = reference["broken_plus"]
    compare(replay, raw)
    for route in raw["routes"]:
        assert min(route["receiver_distances"]) >= 4
    for panel in raw["registered"].values():
        assert panel["all_visible_bit_error_upper"] < 6e-7
        assert abs(panel["signal"]["2"]) < 0.001 and abs(panel["signal"]["1"]) > 0.9
    for panel in raw["pooled"].values():
        for loads, x in panel.items():
            q = int(loads) - 10 / 3
            lo, hi = x["simultaneous_99_interval"]
            assert lo < q < hi
    for n, rows in raw["middle_ambiguity"].items():
        assert all(
            bit["same_view"] > bit["same_view_walls"]
            for bit in rows["indistinguishable"]
        )
        assert rows["segmented_single_pole_variation"] == (0.5 if n == "3" else 1)
    for w in range(1024):
        if external_field.chart(w):
            for radius in (2, 3):
                weights, windows = external_field.shoulders(w, radius)
                arr = np.array([(w >> i) & 1 for i in range(10)])
                assert weights @ arr - windows == 3 * w.bit_count() - 10


def test_released_support_and_compact_census():
    r = stored("escape")
    compare(escape.run(), r)
    for v, limit in zip(r["sectors"], (0.211, 0.247), strict=True):
        assert v["three_whole_owners_upper"] < limit
        assert v["trial_escape_quadratic"] > 1.6
    c = stored("localized")
    compare(localized.run(), c)
    for shape in c["gate_classes"]:
        sampled = 1
        for row in shape["histograms"]:
            a, b = row["counts"]
            sampled += row["words"] * (1 if a == b else 2)
            if a + b >= 2:
                assert all(x["compact_dimension"] == 0 for x in row["native"])
        assert sampled == 3**9
    # Independent real diagonalization with the full material-to-gateway leakage.
    # Eigenvectors in a degenerate line must be tested as a subspace.
    gates = (0, 1, 2)
    for a, b in ((0, 1), (0, 3), (1, 1), (1, 2)):
        words, exits, hop, contrast = localized.sector(gates, a, b)
        n = len(words)
        cut = np.zeros((len(exits), n))
        for i, via in enumerate(exits):
            cut[i, via] = 1
        h = hop.toarray() + np.diag(16 * (np.array([1, 4, 1]) @ contrast))
        e, v = np.linalg.eigh(h)
        tally = 0
        for block in np.split(np.arange(n), np.where(np.diff(e) > 1e-8)[0] + 1):
            tally += len(block) - np.linalg.matrix_rank(cut @ v[:, block], tol=1e-8)
        assert tally == (5 if a + b == 1 else 0)


def test_native_internal_algebra_and_singlet():
    assert registers.run() == stored("../registers")
    checks = stored("../registers")["8"]
    assert checks["free"]["derived_dimension"] == 3
    assert checks["free"]["brackets_confined_to_zero_contrast_coin"]
    assert checks["free_and_echo"]["derived_dimension"] == 0
    r = stored("orbital_singlet")
    compare(json.loads(json.dumps(orbital_singlet.run())), r)
    assert r["local_port_read_commutant"] == 1
    for board in r["comparisons"]:
        near, far = board["placements"]
        assert max(near["distances"]) < min(far["distances"])
        for location in (near, far):
            assert location["unit_energies"] == {"-1": 35, "1": 67}
            assert location["triple_positive_energy"] == 201
            assert location["neutral_pair_energy"] == 102
            assert location["triple_full_residual"] < 1e-10
            assert location["opposite_pair_full_residual"] < 1e-10
