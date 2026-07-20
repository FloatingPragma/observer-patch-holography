"""Regression tests for the independent hierarchy interval witness.

The witness must certify the bracket and the negative derivative enclosure on
the declared interval, and every antecedent-perturbation control must fail
closed. The stored receipt must agree with the checked-in claims.
"""
import json
import pathlib
import sys

import mpmath as mp

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "computations"))

import hierarchy_interval_witness as w  # noqa: E402

iv = w.iv


def test_bracket_sign_change_fast():
    b1, b2, b3 = iv.mpf(33) / 5, iv.mpf(1), iv.mpf(-3)
    f_lo = w.phi(iv.mpf(w.I_LO), w.P_STR, b1, b2, b3, mu_iters=48)
    f_hi = w.phi(iv.mpf(w.I_HI), w.P_STR, b1, b2, b3, mu_iters=48)
    assert w.sign_of(f_lo) == 1
    assert w.sign_of(f_hi) == -1


def test_derivative_enclosure_negative_one_piece():
    b1, b2, b3 = iv.mpf(33) / 5, iv.mpf(1), iv.mpf(-3)
    mid = (mp.mpf(w.I_LO) + mp.mpf(w.I_HI)) / 2
    a_piece = iv.mpf([w.I_LO, str(mid)])
    res = w.phi(w.Dual.var(a_piece), w.P_STR, b1, b2, b3, mu_iters=48)
    assert res.d.b < 0
    assert mp.mpf(res.d.a) > -12 and mp.mpf(res.d.b) < -10


def test_perturbed_beta_breaks_bracket():
    b1, b2, b3 = iv.mpf(33) / 5, iv.mpf(1), iv.mpf(-3)
    f_lo = w.phi(iv.mpf(w.I_LO), w.P_STR, b1, b2 + iv.mpf('0.05'), b3, mu_iters=48)
    f_hi = w.phi(iv.mpf(w.I_HI), w.P_STR, b1, b2 + iv.mpf('0.05'), b3, mu_iters=48)
    assert w.sign_of(f_lo) * w.sign_of(f_hi) != -1


def test_receipt_matches_claims():
    receipt = json.loads((HERE / "certificates" / "independent_interval_witness_receipt.json").read_text())
    assert receipt["bracket_sign_change"] is True
    assert receipt["derivative_strictly_negative"] is True
    assert receipt["unique_root_certified"] is True
    assert receipt["controls_all_fail_closed"] is True
    lo, hi = (mp.mpf(x) for x in receipt["derivative_enclosure_union"].strip("[]").split(","))
    assert lo < mp.mpf("-10.995768") and hi > mp.mpf("-10.985284"), (
        "independent enclosure must contain the published enclosure")
    assert hi < 0
