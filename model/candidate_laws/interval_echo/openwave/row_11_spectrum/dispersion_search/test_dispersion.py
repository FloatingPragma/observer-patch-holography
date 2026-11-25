import json

import numpy as np
from scipy.linalg import eigh, expm

from . import measure
from .apparatus import Guide


def test_nonadjacent_native_phase_and_interferometer():
    g = Guide(16, 8)
    n = g.n
    a, b = 4, 7
    gates = []
    for j, typ in enumerate(g.word):
        h = np.zeros((n, n))
        h[j : j + 2, j : j + 2] = (
            [[-4, -4], [-4, -4]] if typ == "1" else [[8.5, -0.5], [-0.5, 8.5]]
        )
        gates.append(h)
    shifted = []
    for ref in (0, 1):
        orig = g.sources[ref]
        u = np.eye(n, dtype=complex)
        for k, target in ((1, b), (0, a)):
            for j in range(orig + k, target):
                full = expm(-1j * (np.pi if g.word[j] == "0" else np.pi / 8) * gates[j])
                u = full @ u
        shifted.append(u @ gates[orig] @ u.conj().T)
    sb, sm = [g.transported_sign(ref, a, b) for ref in (0, 1)]
    angle = 0.7381
    tb = (-angle / (8.5 + 0.5 * sb / sm)) % (2 * np.pi)
    tm = (-tb * sb / (8 * sm)) % (np.pi / 4)
    actual = expm(-1j * tm * shifted[1]) @ expm(-1j * tb * shifted[0])
    expected = np.eye(n, dtype=complex)
    expected[a, a] = expected[b, b] = np.exp(1j * angle)
    np.testing.assert_allclose(actual, expected, atol=5e-13)
    vec = np.arange(1, n + 1) * np.exp(0.37j * np.arange(n))
    vec /= np.linalg.norm(vec)
    built = g.prepare(vec)
    probabilities = g.read_pair(built, a, b)
    np.testing.assert_allclose(
        probabilities,
        abs(np.array([built[a] + built[b], built[a] - built[b]])) ** 2 / 2,
    )


def test_replay_disjoint_banks_and_finite_bulk():
    saved = json.loads((measure.HERE / "records.json").read_text())
    assert saved["selected_period"] == 16
    assert [
        sum(s["compatible"] for s in saved["negative_training_joint_holdouts"][m])
        for m in ("relativistic", "quadratic", "lattice", "lattice_shell")
    ] == [1, 1, 17, 17]
    for r in saved["observations"]:
        g = Guide(r["n"], 16)
        assert r["reference_index"] == 2 and len(r["modes"]) == r["n"] // 16 - 2
        for mode in r["modes"]:
            mid, c = measure.signature(g, mode["index"])
            theta = np.arccos(c)
            assert mode["phase_interval99"][0] < theta < mode["phase_interval99"][1]
            delta = g.levels[mode["index"]] - g.levels[r["reference_index"]]
            lo, hi = mode["recoil"]["interval99"]
            assert lo <= delta <= hi
        for receiver in r["source_receiver"]["receivers"]:
            assert receiver["signed_load"] == r["charge"]
            assert all(
                abs(b["local_load"] - b["reciprocal_analyzer"])
                < 2 * measure.PROB_ERROR * measure.SHOTS
                for b in receiver["readings"]
            )
        assert r["charged_chart"]["maximum_sharp_charge_dimension"] == 1
    # Exercise the actual source, momentum, Ramsey and dual-current readers
    # against one frozen independent seed without replaying every long layout.
    ports, board = measure.geometry()
    for length, expected_cut in ((96, 352), (144, 526)):
        path = set(ports[:length])
        assert sum((a in path) != (b in path) for a, b in board.seams) == expected_cut
    repeated = measure.acquire(
        96, 16, -1, ports[:96], board, np.random.default_rng(830021 + 96 - 1)
    )
    frozen = saved["observations"][0]
    # Both the deposited and freshly sampled records must agree individually
    # with the independent law probabilities. A numerical eigensolver can
    # change the random variates consumed by earlier multinomial banks.
    g = Guide(96, 16)
    h = np.diag(g.diag) + np.diag(g.off, 1) + np.diag(g.off, -1)
    energies, modes = eigh(h, driver="evd")
    source_phase_port = 96 - 1 - 24
    for mode in modes.T:
        if mode[source_phase_port] < 0:
            mode *= -1
    for observed in (frozen, repeated):
        verify_count_banks(observed, energies, modes)

    # Refit only the same two short-guide vacancy settings. Check each held
    # conclusion at its original sign, length and momentum, for all four laws.
    expected = saved["negative_training_joint_holdouts"]
    conclusions = {law: {} for law in expected}
    for held in saved["observations"][1:]:
        refitted = measure.fit(repeated, held)
        for law in expected:
            for setting in refitted[law]["held"]:
                key = (setting["charge"], setting["length"], setting["mode"])
                conclusions[law][key] = setting["compatible"]
    for law, settings in expected.items():
        assert conclusions[law] == {
            (s["charge"], s["length"], s["mode"]): s["compatible"] for s in settings
        }


def verify_count_banks(observed, energies, modes):
    """Evaluate the held readings from the two-mode spectral law in a dense basis."""
    assert (observed["charge"], observed["n"], observed["reference_index"]) == (
        -1,
        96,
        2,
    )
    assert [m["index"] for m in observed["modes"]] == [2, 3, 4, 5]
    assert observed["charged_chart"]["reverse_transfer_by_arrived_head"] == [0.0, 1.0]

    def check_count(count, probability):
        lo, hi = measure.interval(count, measure.SHOTS)
        assert lo - 1e-12 <= probability <= hi + 1e-12

    for mode in observed["modes"][1:]:
        delta = energies[mode["index"]] - energies[2]
        recoil = mode["recoil"]
        assert recoil["interval99"][0] <= delta <= recoil["interval99"][1]
        for t, (in_phase, quadrature) in zip(
            measure.WAITS, recoil["counts"], strict=True
        ):
            check_count(in_phase, (1 + np.cos(delta * t)) / 2)
            check_count(quadrature, (1 + np.sin(delta * t)) / 2)

    crossing = observed["source_receiver"]
    assert crossing["mode"] == 5
    delta = energies[5] - energies[2]
    for receiver, (i, port) in zip(
        crossing["receivers"], ((24, 1832), (72, 2078)), strict=True
    ):
        assert (
            receiver["index"],
            receiver["port"],
            receiver["signed_load"],
            receiver["incoming_head"],
        ) == (i, port, -1, 1)
        a, b = modes[i, 2], modes[i, 5]
        incoherent = (a * a + b * b) / 2
        check_count(receiver["incoherent_source"], incoherent)
        for t, row in zip(
            (0, 32768, 131072, 524288), receiver["readings"], strict=True
        ):
            assert row["t"] == t
            probability = incoherent - a * b * np.cos(delta * t)
            for instrument in (
                "local_load",
                "reciprocal_analyzer",
                "arrival_and_native_head",
            ):
                check_count(row[instrument], probability)


def test_exact_charged_energy_fiber_and_spin_controls():
    from . import spin_fibers

    report = json.loads((measure.HERE / "spin_fibers.json").read_text())
    actual = spin_fibers.run()
    for copied, case in zip(report["fixtures"], actual["fixtures"], strict=True):
        assert copied["cyclicity"] == case["cyclicity"]
        assert case["cyclicity"]["span_dimension"] == case["length"]
        assert case["aperture_certificate"]["gcd_degrees"] == [0, 0]
        assert all(
            row["exact_receiver_dimension"] == 1
            and row["source_receiver_singular_values"][1] == 0
            for row in case["pole_readouts"]
        )
        source = case["declared_spin_control"]
        assert source["half_turn_commutator_sign"] == -1
        assert (
            source["known_split_control"]["individual_fixed_energy_spin_dimension"] == 1
        )
        for detector in source["analyzers"]:
            assert min(detector["spin_sensitive_singular_values"]) > 0.04
            assert detector["simultaneous99_success_lower"] > 0.001
