"""Verify stationary exclusions from original counts and the integer propagation."""

import json

import numpy as np
from scipy.linalg import expm
from scipy.sparse import csr_matrix

from openwave.interval_echo.row_19.selected_collar import (
    certified_pointer_flight as integer,
)
from openwave.interval_echo.row_19.selected_collar import (
    equilibrium_feedback as feedback,
)

HERE = feedback.HERE


def read(name):
    return json.loads((HERE / (name + ".json")).read_bytes())


def test_integer_columns_on_a_noncommuting_head():
    # One address and four head values. Compare both integer and sparse
    # floating propagation with an independent dense matrix exponential.
    H = np.array(
        [[2, 1, 0, -1], [1, -3, 2, 1], [0, 2, 4, 0], [-1, 1, 0, -1]], dtype=np.int64
    )
    sparse = csr_matrix(H)
    fixed = integer.integer_columns(sparse, 0, 6)
    floating = feedback.service.short_steps(-0.25j * sparse, np.eye(4), 6)
    for tick, row in enumerate(fixed["ticks"], 1):
        err = np.linalg.norm(floating[tick] - expm(-0.25j * tick * H), ord=2)
        assert err <= row["flight_operator_error"]


def test_full_pointer_disk_including_every_initially_unresolved_setting():
    tape = read("equilibrium_pointer_price")
    known = read("certified_pointer_flight")
    assert len(tape["history_results"]) == len(known["per_history"]) == 144
    revised = 0
    for key, row in known["per_history"].items():
        prior = tape["history_results"][key]["minimum_stationary_drift"]
        assert prior == row["prior_lower"]
        value, _ = integer.lower_drift(tape, key, row["numerical_allowance"])
        assert abs(value - row["lower"]) < 1e-15
        revised += prior <= 0
    assert revised == 4
    assert (
        min(r["lower"] for r in known["per_history"].values())
        == known["uniform_drift_lower"]
        > 8.4109e-9
    )


def test_feedback_rule_from_its_two_registered_branches():
    x = read("equilibrium_feedback")
    with np.load(HERE / "equilibrium_feedback_counts.npz", allow_pickle=False) as z:
        counts = z["terminal"]
        raw = z["reference"]
    panel = feedback.drift_panel(raw, counts)
    assert panel.shape == (24336, 30)
    w = np.array(x["outcomes"]["blind"]["signed_Pauli_witness"])
    drift = np.min(panel @ w) - np.linalg.norm(w) * x["drift_HS_radius"] - 1e-11
    assert abs(drift - x["outcomes"]["blind"]["minimum_drift"]) < 2e-15
    assert drift > 5.5368e-7
    # At the same duration, integrating over the complete midpoint record
    # reproduces a direct end measurement, with no conditional renormalization.
    loc, sign, t, h = 2, 1, 3, 8
    also = counts[loc, sign, t - 1, h].sum(axis=2)
    center, _ = feedback.service.center(also)
    reference, _ = feedback.service.center(raw)
    which = (((loc * 2 + sign) * 12 + t - 1) * 13 + h) * 13 + h
    np.testing.assert_allclose(center - reference, panel[which], atol=5e-16)
