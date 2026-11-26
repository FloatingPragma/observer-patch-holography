"""Physical columns, independent evolution, receiver distributions and charge control tests."""

import json
from functools import cache

import numpy as np
from scipy.sparse.linalg import expm_multiply

from . import measure, readers, ring


@cache
def receipt():
    return json.loads((measure.HERE / "records.json").read_text())


def test_typed_contacts_and_complete_independent_short_flight():
    words, sites = ring.basis(14)
    for value, k in ((1, ring.K), (2, ring.K), (1, 0)):
        h = ring.generator(14, value, k)
        bonds = (sites != np.roll(sites, -1, axis=1)).sum(axis=1)
        diag = (0.5 + k * value**2) * bonds
        np.testing.assert_allclose(
            h.diagonal() - diag, (h.diagonal() - diag)[0], atol=1e-12
        )
        # Independent scalar normalization and Krylov propagation on the
        # complete population space, then compare every word at the held time.
        start = np.zeros(len(words))
        start[np.searchsorted(words, ring.launch(14))] = 1
        out = expm_multiply(-3j * h, start)
        with np.load(measure.HERE / "words.npz") as record:
            name = "native" if k else "contrast_off"
            np.testing.assert_allclose(
                abs(out) ** 2, record[f"law_14_{value}_{name}"][4], atol=1e-8, rtol=1e-6
            )
    for theta in ring.ANGLES:
        z = [readers.phase(a, theta) / readers.phase(0, theta) for a in range(3)]
        np.testing.assert_allclose(z, np.exp(1j * theta * np.arange(3)), atol=1e-12)


def test_routed_words_and_joint_inclusive_counts():
    with np.load(measure.HERE / "words.npz") as bank:
        for row in receipt()["boards"]:
            n, level = row["n"], row["level"]
            assert row["receivers"] == json.loads(
                json.dumps(readers.stations(level, n))
            )
            for amount in (1, 2):
                for sign in (-1, 1):
                    gates, two, _, _ = ring.regional(n, sign, amount)
                    summary = row["release"][f"{amount}_{sign}"]
                    for name in ("native", "contrast_off"):
                        counts = bank[f"counts_{level}_{n}_{amount}_{sign}_{name}"]
                        np.testing.assert_array_equal(counts.sum(axis=1), measure.SHOTS)
                        for mask, tag in ((gates, "joint_gate"), (two, "two_wall")):
                            np.testing.assert_array_equal(
                                counts[:, mask].sum(axis=1),
                                summary[name][tag + "_counts"],
                            )
                            np.testing.assert_allclose(
                                bank[f"law_{n}_{amount}_{name}"][:, mask].sum(
                                    axis=1, dtype=float
                                ),
                                summary[name][tag + "_probabilities"],
                                atol=1e-7,
                            )
                    assert summary["native"]["joint_gate_probabilities"][4] > 0.4
                    assert summary["contrast_off"]["joint_gate_probabilities"][4] < 0.02
                    assert min(summary["native"]["two_wall_probabilities"]) > 1 - n / (
                        4 * (ring.K * amount**2 + 0.5)
                    )


def test_signed_voltage_is_not_a_gain_assigned_to_a_wall():
    for board in receipt()["boards"]:
        assert measure.same_record(
            measure.read_motion(board["held_motion"], board["standards"]),
            board["inferred_motion"],
        )
        for key, scans in board["inferred_motion"].items():
            value, sign = map(int, key.split("_"))
            for s in scans:
                assert abs(s["charge_per_wall_step"] - sign * value / 2) < 0.002
                assert abs(s["ordinary_per_site_step"] - value) < 0.003
                assert abs(s["parked_per_wall_step"]) < 0.001
        for source in board["held_motion"].values():
            for row in source:
                assert (row["final_bond"] - row["first_bond"]) % board["n"] == 2
                assert (row["displaced"] ^ row["initial"]).bit_count() == 2


def test_full_turn_profiles_from_the_native_voltage_not_wall_labels():
    n = 14
    _, sites = ring.basis(n)
    for sign in (-1, 1):
        occupation = sites if sign == 1 else 1 - sites
        for amount in (1, 2):
            gates, _, taper, flat = ring.regional(n, sign, amount)
            for j, arm in enumerate(((-5, -4, -3), (2, 3, 4))):
                shape = np.zeros(n)
                for b in arm:
                    shape[b % n] += 0.5
                    shape[(b + 1) % n] += 0.5
                for weights, expected in (
                    (shape, taper[:, j] / 2),
                    ((shape > 0).astype(int), flat[:, j]),
                ):
                    for theta in (np.pi / 2, np.pi, 2 * np.pi):
                        factor = np.array(
                            [
                                readers.phase(amount, theta * w)
                                / readers.phase(0, theta * w)
                                for w in weights
                            ]
                        )
                        relative = np.prod(np.where(occupation, factor, 1), axis=1)
                        for parity in (0, 1):
                            reference = factor[(np.arange(n) + parity) % 2 == 1].prod()
                            np.testing.assert_allclose(
                                relative / reference,
                                np.exp(1j * theta * expected),
                                atol=2e-12,
                            )
                assert np.all(taper[gates, j] == amount * sign)
            board = receipt()["boards"][0]["release"][f"{amount}_{sign}"]
            for profile in ("flat", "tapered"):
                result = board["two_voltage_probes"][profile][3][0]
                accepted = result["counts"][4:]
                desired = 3 if amount == 1 and profile == "tapered" else 0
                assert sum(accepted) > 0 and sum(accepted) == accepted[desired]


def test_closed_stock_transfer_and_remote_stationarity():
    from . import stock

    result = stock.report()
    assert measure.same_record(
        result, json.loads((measure.HERE / "stock.json").read_text())
    )
    for case in result:
        for trial in case["stocks"]:
            assert (
                trial["ring_after"] - trial["ring_before"]
                == trial["sign"] * trial["amount"]
            )
            assert trial["ring_after"] + trial["reservoir_after"] == trial["conserved"]
