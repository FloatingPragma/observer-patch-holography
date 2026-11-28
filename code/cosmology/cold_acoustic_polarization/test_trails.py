"""Independent checks of the conserved source ratio and its held comparison."""

import hashlib
import json
from fractions import Fraction

import numpy as np

from .boltzmann import compute
from .measure import HERE, Likelihood, whiten
from .trails import RECEIPT


def test_all_native_levels_and_chain_position():
    result = json.loads((HERE / "trail_results.json").read_text())
    source = json.loads(RECEIPT.read_text())
    assert (
        hashlib.sha256(RECEIPT.read_bytes()).hexdigest()
        == result["native_receipt_sha256"]
    )
    assert [b["live_seams"] for b in result["boards"]] == [14, 32, 64]
    with np.load(HERE / "chain_density_columns.npz") as samples:
        for board, native in zip(result["boards"], source["boards"], strict=True):
            m = len(native["path"]) - 1
            q = Fraction(34 * m - 51, 202)
            assert (
                q
                == Fraction(board["ratio_exact"])
                == Fraction(native["net_work"]) / Fraction(native["initial_energy"])
            )
            assert Fraction(native["work_from_sharp_reads"]) == Fraction(17 * m, 4)
            for channel in ("TE", "EE"):
                rows = np.vstack([samples[f"{channel}_{i}"] for i in range(1, 5)])
                w, b, c = rows.T
                rank = np.average(c <= float(q) * b, weights=w)
                assert (
                    rank
                    == board["published_chain"][channel][
                        "fraction_of_weight_at_or_below"
                    ]
                )
                assert board["published_chain"][channel]["within_central_68"] == (
                    m == 32
                )


def test_independent_middle_separation_and_score_difference():
    result = json.loads((HERE / "trail_results.json").read_text())
    original = json.loads((HERE / "results.json").read_text())
    middle = result["boards"][1]["baryon_density_runs"]["0.02237"]
    lk = Likelihood()
    theory, zeq = compute(0.02237, 1037 / 202)
    fit = lk.calibrate(theory[:, 0], True, 0.2)
    np.testing.assert_allclose(zeq, middle["z_equality"], rtol=1e-10)
    predicted = lk.transformed(theory, fit)
    fd = np.tile(lk.eff * (lk.eff + 1) / (2 * np.pi), 3)
    standard = (
        np.ravel(original["selected_bands_uK2"]["predictions"]["standard_control"]) / fd
    )
    # Express both fixed predictions in common conditional coordinates.
    delta = predicted - standard
    delta = whiten(lk.cond, delta[lk.n :] - lk.B @ delta[: lk.n])
    r = lk.y - standard
    residual = whiten(lk.cond, r[lk.n :] - lk.B @ r[: lk.n])
    np.testing.assert_allclose(
        delta @ delta,
        middle["conditional_prediction_squared_separation_from_standard"],
        rtol=1e-5,
    )
    np.testing.assert_allclose(
        delta @ delta - 2 * residual @ delta,
        middle["conditional_P_excess_score_over_standard"],
        rtol=1e-5,
    )
