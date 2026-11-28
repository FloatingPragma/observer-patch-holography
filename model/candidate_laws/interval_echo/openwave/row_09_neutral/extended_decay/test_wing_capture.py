"""Check the routed capture, receiver traces and filled-host obstruction."""

from hashlib import sha256
from itertools import combinations
import json
import numpy as np
from scipy.linalg import eigh
from . import filled_outlet
from . import open_ports, wing_capture, wing_tracer
from .measure import HERE
from ..passive_decay.measure import hamiltonians


def record(module):
    data = json.loads((HERE / f"{module}.json").read_text())
    assert (
        data["protocol_sha256"]
        == sha256((HERE / f"{module.upper()}.md").read_bytes()).hexdigest()
    )
    return data


def test_coherent_swap_preserves_the_full_tagged_state_with_an_external_reference():
    hardware = wing_capture.capture_hardware()
    assert hardware == record("wing_capture")["hardware"]
    dimers = hardware["guide_wing"]
    terminals = dimers[1]
    outside = [5, 10]
    words = list(combinations([*terminals, *outside], 2))
    _, actions = wing_tracer.tracer_maps(words, terminals)
    _, frame = eigh(hamiltonians(4)[1])
    ground = frame[:, 0]
    rng = np.random.default_rng(86311)
    raw = rng.normal(size=(2 * len(words), 3)) + 1j * rng.normal(
        size=(2 * len(words), 3)
    )
    raw /= np.linalg.norm(raw)
    tag = actions["any_excited"] @ raw
    # Explicitly perform E<->E', D<->E, E<->D', D<->E on the six-index
    # amplitude (D,E,D',E',outside,reference), keeping the reference arbitrary.
    after = np.zeros((2 * len(words), 3, 4), complex)
    for q in outside:
        rows = [
            d * len(words)
            + next(j for j, w in enumerate(words) if set(w) == {terminals[1 - e], q})
            for d in (0, 1)
            for e in (0, 1)
        ]
        before = tag[rows].reshape(2, 2, 3)
        registers = np.einsum("der,ab->deabr", before, ground.reshape(2, 2))
        for a, b in ((1, 3), (0, 1), (1, 2), (0, 1)):
            registers = registers.swapaxes(a, b)
        expected = np.einsum("de,abr->deabr", ground.reshape(2, 2), before)
        np.testing.assert_allclose(registers, expected, atol=4e-13)
        after[rows] = registers.reshape(4, 4, 3).transpose(0, 2, 1)
    actual = after.reshape(-1, 4)
    actual = actual @ actual.conj().T
    marginal = np.zeros_like(actual)
    for tag in wing_tracer.EXTRACTED:
        branch = (actions[tag] @ raw).ravel()
        marginal += np.outer(branch, branch.conj())
    np.testing.assert_allclose(actual, marginal, atol=4e-13)


def test_unmeasured_traces_recover_the_existing_source_and_receiver():
    tagged = record("wing_tracer")
    captured = record("wing_capture")
    original = record("joined_outlet")
    for name, central in tagged["settings"].items():
        unmarked = original["settings"][name]["arms"]
        for index, page in enumerate(
            (
                captured["settings"][name]["162"],
                central,
                captured["settings"][name]["164"],
            )
        ):
            for case, arms in page["cases"].items():
                clock = 1 + index + (3 if case == "decay" else 0)
                for arm, reads in arms.items():
                    np.testing.assert_allclose(
                        reads["untouched"]["core_by_outlet"],
                        unmarked[arm]["core_by_outlet"][clock],
                        atol=5e-11,
                        rtol=0,
                    )


def test_held_calibration_is_sensitive_to_capture_and_the_source_fails():
    evidence = record("wing_capture")
    for name, pages in evidence["settings"].items():
        for page in pages.values():
            r = {}
            for c, arms in page["cases"].items():
                pulse = arms["D_phase" if c == "known" else "core_pulse"]
                r[c] = wing_capture.contrast(
                    pulse, arms["idle"], page["isolated_wing_gaps"]
                )
            if name == "native":
                assert 0.017 < r["known"][4] < 0.022
                assert r["decay"][0] < 0 and r["decay"][5] < 0
                assert 0 < r["decay"][4] < 0.0008
            if name == "exit_cut":
                np.testing.assert_allclose([r[c][4] for c in r], 0, atol=1e-10)
            if name in ("source_cut", "no_reversal"):
                np.testing.assert_allclose(r["decay"][:5], 0, atol=1e-10)


def test_filled_host_has_only_a_vacancy_and_its_pair_floor_is_physical():
    evidence = record("filled_outlet")
    b, _, cell, _, _, contacts, _ = open_ports.aperture()
    h, _, _, inside, words = filled_outlet.system("native")
    assert h.shape == (304, 304)
    for i, word in enumerate(words):
        missing = set(cell) - set(word)
        assert len(missing) == int(inside[i])
    internal = [e for e in contacts if set(e) <= set(cell)]
    assert len(internal) == 30
    # Direct full-word potential difference from the homogeneous host. A blank
    # exchange contributes positive semidefinite (I-SWAP)/2 on every edge.
    costs = []
    for vacancy in cell:
        for excess in set(cell) - {vacancy}:
            load = np.ones(b.n, dtype=int)
            load[vacancy], load[excess] = 0, 2
            local = sum((load[a] - load[z]) ** 2 for a, z in internal)
            private = (load**2 - 1) @ (6 - b.degree)
            costs.append(filled_outlet.KAPPA * (local + private))
    pair = evidence["balanced_pair"]
    assert min(costs) == pair["isolated_potential_floor"] == 72
    assert 75 < pair["isolated_ground"] < 76
    for entry in pair["global_boards"]:
        board = filled_outlet.federation({20: 0, 80: 1, 320: 2}[entry["cells"]])
        load = np.ones(board.n, dtype=int)
        vacancy, excess = entry["vacancy_excess_witness"]
        load[vacancy], load[excess] = 0, 2
        a, z = np.asarray(board.seams).T
        cost = (load**2 - 1) @ (6 - board.degree) + ((load[a] - load[z]) ** 2).sum()
        assert filled_outlet.KAPPA * cost == entry["full_contrast_min"] == 80
    assert evidence["vacancy_calibration"]["native"][6] > 0.015
    assert max(evidence["vacancy_calibration"]["bridge_cut"]) == 0
