"""Native columns, full-board cuts, time histories and the small stationary witness."""

from fractions import Fraction
from hashlib import sha256
import json
import numpy as np

from model.candidate_laws.interval_echo import observations as obs
from . import local, measure, rate


def archived():
    return json.loads((measure.HERE / "results.json").read_text())


def test_full_native_service_and_routed_binary_rails():
    weights, write, half = local.operators()
    assert (write[1] != 0).sum() == 2 and write[1].sum() == 64
    assert half[1, 1] == half[1, 2] == 32
    assert local.LEFT_CHARGE[1] != local.LEFT_CHARGE[2]
    for source in range(28):
        for key, table in local.samplers().items():
            observed = np.bincount(table[source], minlength=28) / table.shape[1]
            np.testing.assert_array_equal(observed, weights[key][source] / 4096)
    for level in range(3):
        b = obs.federation(level)
        words = obs.gateway_words(b)
        for gateway in words:
            src = (np.arange(12)[None, :] % 3).astype(np.uint8)
            src[0, obs.MATTER] = int(gateway % 2)
            trial = src.copy()
            obs.route(trial, np.array([0]), np.array([gateway]), words)
            expected = src.copy()
            expected[0, [obs.MATTER, gateway]] = src[0, [gateway, obs.MATTER]]
            np.testing.assert_array_equal(trial, expected)
            obs.route(trial, np.array([0]), np.array([gateway]), words)
            np.testing.assert_array_equal(trial, src)
    assert all(x == 0 for x in measure.control_checks().values())


def test_every_sharp_echo_class_preserves_its_local_ledger():
    # Work directly with every full 99-state echo column whose input lies in
    # the closed binary sector, and measure the two conserved quantities.
    ids = np.array(
        [
            (0 if m is None else 1 + local.law.LABELS.index(m)) * 9
            + local.law.SRC.index((a, b))
            for m in local.MEM
            for a in (0, 1)
            for b in (0, 1)
        ]
    )
    u = local.law.echo()
    for i, col in enumerate(ids):
        support = abs(u[:, col]) > 1e-12
        assert not np.any(support[np.setdiff1d(np.arange(99), ids)])
        active = support[ids]
        assert np.all(local.LEFT_CHARGE[active] == local.LEFT_CHARGE[i])
        assert np.all(local.RIGHT_CHARGE[active] == local.RIGHT_CHARGE[i])


def test_all_board_counts_and_the_fresh_stationary_standard():
    r = archived()
    assert (
        r["protocol_sha256"]
        == sha256((measure.HERE / "PROTOCOL.md").read_bytes()).hexdigest()
    )
    for level in range(3):
        b = obs.federation(level)
        ref = r["boards"][str(b.cells)]
        layers, _, patch, signed, _ = measure.hardware(b)
        assert ref["boundary"] == sum(abs(signed))
        assert ref["fresh_exact"] == measure.fresh_prediction(layers, patch, signed)
        for label, arm in ref["trials"].items():
            for key, row in arm.items():
                net = (
                    row["cumulative_births"][0]
                    - row["cumulative_returns"][0]
                    - row["cumulative_export"][0]
                )
                assert abs(net - row["stored_unit"][0]) < 2e-12
                if label == "retained":
                    assert row["charge_residual_abs"] == [0, 0]
                    assert (
                        row["load_contrast"][0] >= ref["echo_only_load_contrast_lower"]
                    )
                    assert (
                        abs(
                            row["original_patch_loss"][0]
                            - ref["boundary"] * row["cut_current"][0]
                        )
                        < 1e-11
                    )
                if label == "fresh":
                    assert row["cumulative_returns"] == [0, 0]
                    for op, target in (
                        ("load_contrast", "load_contrast"),
                        ("signed_unit_head", "directed_unit_head"),
                    ):
                        mean, se = row[op]
                        assert (
                            abs(mean - ref["fresh_exact"][key][target])
                            <= 5 * se + 1e-10
                        )


def test_compact_uniform_classes_and_positive_rate_limit():
    panel = archived()["compact"]
    summary = rate.derive()
    saved = json.loads((measure.HERE / "rate.json").read_text())
    measure.compare(saved, summary)
    fractions = panel["retained"]["exact_uniform_stationary"]
    assert Fraction(fractions["near_minus_far"]) == Fraction(32, 57)
    assert Fraction(fractions["signed_cut"]) == Fraction(41, 57)
    assert panel["retained"]["blank_patches_in_initial_class"] == 1
    for row in summary["rows"].values():
        assert (
            row["recurrent_size"] == 666 and row["blank_patches_in_initial_class"] == 6
        )
        assert Fraction(row["exact_uniform_stationary"]["near_minus_far"]) == 0
        assert Fraction(row["exact_uniform_stationary"]["signed_cut"]) == 0
    for name, row in panel.items():
        birth, returned, exported = row["stationary_birth_return_export"]
        assert birth > 0 and abs(birth - returned - exported) < 1e-12
        if name != "fresh":
            assert exported == 0
    assert [
        summary["rows"][str(d)]["first_cycle_contrast_below_quarter"]
        for d in (1, 2, 8, 32)
    ] == [8, 10, 25, 84]


def test_replicas_and_the_actual_maintenance_supply():
    from . import supply

    s = supply.cost()
    measure.compare(json.loads((measure.HERE / "supply.json").read_text()), s)
    data = archived()
    assert (
        data["replica_protocol_sha256"]
        == sha256((measure.HERE / "REPLICAS.md").read_bytes()).hexdigest()
    )
    assert data["paired_trajectories"] * 2 == data["independent_copies"]
    for name, b in s["boards"].items():
        assert sum(b["count_numerators"]) == b["count_denominator"]
        numerator = sum(i * x for i, x in enumerate(b["count_numerators"]))
        assert numerator / b["count_denominator"] == b["mean_pairs"] == b["cut"] / 2
        assert b["cut"] == data["boards"][name]["boundary"]
    small = data["compact"]
    assert (
        abs(small["retained"]["stationary_cut_replica_product"] - (41 / 57) ** 2)
        < 1e-13
    )
    assert small["retained_exchange"]["stationary_cut_replica_product"] < 1e-25
    assert small["retained_exchange"]["stationary_cut_unit"] > 0.53
