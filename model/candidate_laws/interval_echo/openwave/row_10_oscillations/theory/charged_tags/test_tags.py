import json
from itertools import product
import numpy as np
import pytest
from .measure import (
    HERE,
    LAW,
    STARTS,
    SOURCES,
    WORDS,
    read,
    predictions,
    acquire,
    instruments,
    ANGLES,
)


def test_all_registered_cohorts_replay():
    record = json.loads((HERE / "counts.json").read_text())
    law = json.loads(LAW.read_text())
    fresh = json.loads(json.dumps(acquire(law)))
    # Trace distances are eigenvalue sums and carry library roundoff.
    computed = fresh.pop("computed")
    assert [(r["strength"], r["step"]) for r in computed] == [
        (r["strength"], r["step"]) for r in record["computed"]
    ]
    np.testing.assert_allclose(
        [r["daughter_trace_distance"] for r in computed],
        [r["daughter_trace_distance"] for r in record["computed"]],
        rtol=0,
        atol=1e-12,
    )
    assert fresh == {k: v for k, v in record.items() if k != "computed"}
    law[0]["probabilities"][7] -= 1e-8
    with pytest.raises(ArithmeticError, match="acquisition law"):
        acquire(law)
    rows = read(record)
    for r in rows:
        if r["k"] != 8 or r["step"] != 2:
            continue
        assert r["charge_clicks"] == [STARTS] * 4
        assert r["phase"]["joint"][0] == [STARTS, 0]
        for mode in ("donor_read", "order_read"):
            lo, hi = r["phase"][mode][2]
            assert lo < 0 < hi


def test_all_daughter_statistics_and_joint_sign_conservation():
    for k, step in product((4, 8), range(3)):
        p, distance = predictions(k, step)
        _, convert, _ = instruments(k, step)
        alpha = k * ANGLES[step] / 8
        for sign in (-1, 1):
            a = WORDS.index((1 + sign, 1 - sign, 1))
            b = WORDS.index((1 + sign, 1, 1 - sign))
            exact = np.zeros(27, complex)
            exact[a] = np.exp(-0.5j * alpha) * np.cos(alpha / 2)
            exact[b] = -np.exp(-0.5j * alpha) * np.sin(alpha / 2)
            np.testing.assert_allclose(convert[:, a], exact, atol=1e-12)

        strength = np.sin(step * np.pi / 2 * k / 8 / 2) ** 2
        np.testing.assert_allclose(distance, 1 - strength, atol=1e-12)
        for label in SOURCES:
            sharp = p[label, "sharp"]
            for w, prob in zip(WORDS, sharp, strict=True):
                if sum(w) != 3:
                    assert prob < 1e-25
        for readout in ("sharp", "local_half_swap"):
            marg = [p[s, readout].reshape(9, 3).sum(axis=1) for s in SOURCES[2:]]
            measured = abs(marg[0] - marg[1]).sum() / 2
            assert measured <= distance + 1e-12
            if readout == "local_half_swap":
                np.testing.assert_allclose(measured, distance, atol=1e-12)


def test_additional_internal_tag_survives_at_one_fixed_charge():
    record = json.loads((HERE / "counts.json").read_text())
    b = {tuple(r["key"]): r["counts"] for r in record["same_charge_control"]}
    assert len(b) == 16
    for site, sign in product(("opposite", "second_owner"), (-1, 1)):
        for reverse in (False, True):
            for readout in ("sharp", "half_swap"):
                c = b[site, sign, reverse, readout]
                assert c[-1] == 0 and sum(c) == STARTS
                assert (
                    sum(
                        n
                        for w, n in zip(WORDS, c[:-1], strict=True)
                        if sum(w[:2]) - 2 == sign
                    )
                    == STARTS
                )
        a = np.array(b[site, sign, False, "half_swap"])
        c = np.array(b[site, sign, True, "half_swap"])
        assert (
            max(a) == max(c) == STARTS
            and min(a + c) == 0
            and sum(abs(a - c)) == 2 * STARTS
        )
