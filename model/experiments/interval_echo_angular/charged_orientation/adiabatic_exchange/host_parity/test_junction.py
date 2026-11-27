"""Exact detector, holonomy, coverage and passive embedding regressions."""

import json
from itertools import product
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from . import experiment as e
from .passive import commission, local_certificate


def test_native_route_and_source():
    for sign in (1, -1):
        p, wave = e.prepare(sign)
        np.testing.assert_allclose(p, 0.5)
        comparison = e.evolved(wave, e.LEFT + e.RIGHT[::-1])
        np.testing.assert_allclose(comparison, sign * wave, atol=2e-12)
        np.testing.assert_allclose(e.evolved(wave, (3,)), comparison, atol=2e-12)
        # The non-unit labels occupy 0 and 3. The source is a native energy line.
        h = e.word_matrix(8, (3,))
        np.testing.assert_allclose(h @ wave, (32 + (1 - sign) / 2) * wave, atol=2e-12)
    assert e.permutation(e.LEFT + e.RIGHT[::-1]) == (3, 2, 1, 0)


def test_closed_path_parities():
    for length in (2, 4, 6, 8):
        for route in product(e.BONDS, repeat=length):
            tags = tuple(range(4))
            for leaf in route:
                if (tags[0] in (1, 2)) == (tags[leaf] in (1, 2)):
                    break
                tags = e.swap(tags, 0, leaf)
            else:
                if {tags[1], tags[2]} == {1, 2}:
                    assert (tags[1] == 2) == (tags[0] == 3)
    # Enumerate all closed single-core walks through depth eight from each port.
    for origin in range(4):
        paths = [(origin, tuple(range(4)))]
        for _ in range(8):
            followers = []
            for current, tags in paths:
                neighbours = e.BONDS if current == 0 else (0,)
                for destination in neighbours:
                    moved = e.swap(tags, current, destination)
                    if destination == origin:
                        assert moved == tuple(range(4))
                    followers.append((destination, moved))
            paths = followers


def test_physical_reflection():
    flip = np.array([e.WORDS.index(e.swap(w, 1, 2)) for w in e.WORDS])
    branch = np.array([e.WORDS.index(tuple(2 - v for v in w)) for w in e.WORDS])
    for k in (4, 8):
        h = e.word_matrix(k)
        np.testing.assert_array_equal(h[flip][:, flip], h)
        np.testing.assert_array_equal(h[branch][:, branch], h)
        for sign in (1, -1):
            _, psi = e.prepare(sign, k)
            np.testing.assert_allclose(psi[flip], psi, atol=2e-12)
            evolved = expm(-11j * h) @ psi
            np.testing.assert_allclose(evolved[branch], sign * evolved, atol=2e-12)
            for i, w in enumerate(e.WORDS):
                if e.positions(w) == (0, 3):
                    assert branch[i] == flip[i]
                    if sign == -1:
                        assert abs(evolved[i]) < 1e-12


def test_spectator_and_replay():
    logged = json.loads(Path(__file__).with_name("results.json").read_text())
    s = e.spectator()
    np.testing.assert_allclose(s["odd_source_with_spectator"], [1, 0.25])
    np.testing.assert_allclose(
        s["mean_exit_spectrum"], [0, 0, 0, 0, 0.25, 0.25, 0.75, 0.75], atol=1e-13
    )
    for saved, actual in zip(logged["free"], e.release(), strict=True):
        np.testing.assert_allclose(saved["law"], actual["law"], atol=6e-12)
        assert all(sum(bank) == e.SHOTS for bank in saved["counts"])
        for bank in saved["counts"]:
            np.testing.assert_allclose(
                np.array(bank) / e.SHOTS, actual["law"], atol=0.001
            )
        if actual["host_sign"] == -1:
            assert all(bank[3] == 0 for bank in saved["counts"])
    np.testing.assert_allclose(logged["free_parked"], [1, 0, 0, 0, 0, 0], atol=1e-12)


def test_unshuttered_mounts():
    assert local_certificate() < 1e-13
    for level in range(3):
        case = commission(level, 3)
        assert case["enrolled"] and len(case["free_held"]) == 20
        assert len(case["quiet_contacts"]) == 13
        assert len(commission(level, 4)["quiet_contacts"]) == 14


def test_charged_vacancies_and_winding():
    from . import vacancies as v

    report = v.run()
    assert len(report["free"]) == 10 and all(
        len(b["cut"]) == 11 for b in report["boards"]
    )
    for sign, single_dark in ((1, 0.25), (-1, 0.75)):
        row = report["addressed"][str(sign)]
        np.testing.assert_allclose(
            [
                row["dark"],
                row["marked_one_winding_dark"],
                row["corrected_dark"],
                row["twice_dark"],
            ],
            [(1 - sign) / 2, single_dark, 0, 0],
            atol=1e-12,
        )
    for r in report["free"]:
        assert r["joint"][0] > 1e-5 and r["joint"][2] < 1e-22
        assert all(b[2] == 0 and sum(b) == v.SHOTS for b in r["counts"])
    for t in (3.17, 9.71):
        psi = expm(-1j * t * v.hamiltonian(v.WORDS, v.RING, 5.3)) @ v.inject(-1, 8)[1]
        np.testing.assert_allclose(psi[list(v.TARGET)], 0, atol=1e-12)
