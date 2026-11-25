from fractions import Fraction as Fr
import json
from math import exp, log, sqrt

import numpy as np

from .surface import CAP, HERE, OWNERS, SHOTS, branches


def test_measured_modes_accept_and_reset_through_registered_contacts():
    native = branches("native")
    np.testing.assert_allclose(native, [[0.5, 0.25, 0.25]] * 2, atol=3e-13, rtol=0)
    np.testing.assert_allclose(
        branches("parked_filter"), [[0.5, 0.5, 0]] * 2, atol=3e-13, rtol=0
    )
    np.testing.assert_allclose(
        branches("parked_proposal"), [[1, 0, 0], [0, 0.5, 0.5]], atol=3e-13, rtol=0
    )
    accept = (Fr(1, 2), Fr(1, 4))
    leak = Fr(1, 4)
    assert tuple(w / (1 - leak) for w in accept) == (Fr(2, 3), Fr(1, 3))
    assert sum(sum(accept) * leak**t for t in range(CAP)) + leak**CAP == 1


def test_first_member_meets_common_ref_visible_and_native_calibration_fields():
    k = OWNERS
    labels = [(owner, n) for owner in range(k) for n in (0, 1)]
    ref = {(a, n): Fr(2 - n, 3 * k) for a, n in labels}
    assert sum(ref.values()) == 1 and all(r > 0 for r in ref.values())
    beta = log(2)
    log_z = log(3 * k / 2) - 32 * beta
    for (owner, n), r in ref.items():
        assert abs(-log(r) - (beta * (32 + n) + log_z)) < 1e-13
        fibre = [(a, s) for a, s in labels if a == owner]
        assert [ref[z] / sum(ref[y] for y in fibre) for z in fibre] == [
            Fr(2, 3),
            Fr(1, 3),
        ]
    # A different correlated owner/sign law with the same prescribed moment.
    p = np.array([float(ref[z]) for z in labels])
    for b in (0, 0.3, beta, 2, 5 * beta):
        mass = 1 / (1 + exp(b))
        gibbs = np.tile([(1 - mass) / k, mass / k], k)
        other = np.zeros_like(gibbs)
        other[0], other[-1] = 1 - mass, mass

        def kl(law, base):
            nz = law > 0
            return np.sum(law[nz] * np.log(law[nz] / base[nz]))

        assert abs(kl(other, p) - kl(gibbs, p) - kl(other, gibbs)) < 1e-13


def test_all_measured_attempts_timeouts_and_independent_reference():
    receipt = json.loads((HERE / "surface.json").read_text())
    assert (receipt["owners"], receipt["shots_per_owner_and_sign"], receipt["cap"]) == (
        OWNERS,
        SHOTS,
        CAP,
    )

    # For three modes, 20 owners, two old signs and nine bins (1,080
    # means) plus pooled read comparisons: conservative 99% Hoeffding.
    def radius(n):
        return sqrt(log(2 * 1200 / 0.01) / (2 * n))

    for mode, owners in receipt["banks"].items():
        assert [o["owner"] for o in owners] == list(range(OWNERS))
        probs = branches(mode)
        for owner in owners:
            assert sorted(r["old"] for r in owner["rows"]) == [0, 1]
            for row in owner["rows"]:
                history = np.r_[np.array(row["accepted"]).ravel(), row["timeout"]]
                assert sum(history) == SHOTS
                left, expected = 1, []
                for t in range(CAP):
                    pg, pe, pl = probs[row["old"] if t == 0 else 1]
                    expected.extend((left * pg, left * pe))
                    left *= pl
                assert max(abs(history / SHOTS - np.r_[expected, left])) < radius(SHOTS)
    native = receipt["banks"]["native"]
    src = json.loads((HERE / "records.json").read_text())
    board = next(b for b in src["boards"] if b["owners"] == OWNERS)
    panel = next(p for p in board["panels"] if p["mode"] == "native")
    bank = next(row for row in panel["banks"] if row["steps"] == 1)
    raw = np.asarray(bank["counts"]).sum(axis=(0, 1))
    sg, se = raw[-2:]
    for old in (0, 1):
        accepted = sum(
            (np.sum(owner["rows"][old]["accepted"], axis=0) for owner in native)
        )
        g, e = accepted
        # Conditional on each realized admission count, terminal signs
        # are iid: the ideal balanced mixture was launched stratified
        # by sharp A/B and each stratum has the same spectral weights.
        assert abs(e / (g + e) - se / (sg + se)) < radius(g + e) + radius(sg + se)
