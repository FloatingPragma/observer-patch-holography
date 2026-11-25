from copy import deepcopy
from fractions import Fraction
import json

import numpy as np
import pytest
from scipy import sparse

from .. import open_pair_check
from ..corridor_pair import neutral_column
from . import experiment as study


def test_exact_projector_obstruction_and_intervention_domain():
    certificate = study.local_certificate()
    assert certificate["exact_null_action"] and certificate["neutral_memory_dimension"] == 18
    assert certificate["native_null_action_max_error"] < 5e-14
    law = study.candidate.native_law
    # An occupied unequal-load contact is outside the erasure domain.
    # Its neutral word is fixed by the echo itself, so the design needs 00.
    occupied = neutral_column((1, 0))
    np.testing.assert_allclose(law.echo() @ occupied, occupied, rtol=0, atol=5e-14)
    h0, _ = law.local_term_parts()
    selected = np.repeat([False, *(d == 1 for _, d in law.LABELS)], law.N_SRC)
    projector = np.diag(selected.astype(float))
    assert np.linalg.norm(h0 @ projector - projector @ h0) > .1


def independent_primary():
    """Raw NPZ pair graph and Chebyshev propagator from the separate row checker."""
    exchange, contrast, launch, _, detector, _ = open_pair_check.board_generator()
    # Add the erased contact's ordinary binary swaps directly, without the
    # intervention producer's Carrier, PairSector or generator functions.
    n, left, right = 240, 152, 166
    stationary = np.array([p for p in range(n) if p not in (left, right)])

    def indices(a, b):
        a, b = np.minimum(a, b), np.maximum(a, b)
        return a * (2 * n - a - 1) // 2 + b - a - 1

    aa, bb = indices(left, stationary), indices(right, stationary)
    rows, cols = np.r_[aa, bb], np.r_[bb, aa]
    off = sparse.csr_matrix((np.full(len(rows), -.5), (rows, cols)), shape=exchange.shape)
    diagonal = np.zeros(exchange.shape[0])
    diagonal[np.r_[aa, bb]] = .5
    changed = exchange + off + sparse.diags(diagonal, format="csr", dtype=np.float64)
    amplitudes = open_pair_check.flight_amplitudes(changed, contrast + 2*diagonal,
                                                  0, (37,), launch, detector)[0]
    home = float(open_pair_check.RETRY) * abs(amplitudes @ open_pair_check.VECTOR)**2
    return np.r_[home, 1-home.sum()]


def test_primary_release_with_independent_operator_and_propagator():
    archived = json.loads(study.ARCHIVE.read_text())
    row = next(r for r in archived["appointments"] if r["kappa"] == 0 and r["exposure"] == 37)
    np.testing.assert_allclose(independent_primary(), row["bins"]["native_one_echo"],
                               rtol=0, atol=study.TOLERANCE)


def test_all_recorded_arms_keep_source_exhaustion_and_the_full_detector():
    archived = json.loads(study.ARCHIVE.read_text())
    assert len(archived["appointments"]) == 9
    assert archived["primary_gate"]["passed"]
    for row in archived["appointments"]:
        for arm, bins in row["bins"].items():
            assert len(bins) == 31 and min(bins) >= 0 and abs(sum(bins)-1) < study.TOLERANCE
            assert bins[-1] >= archived["source_exhaustion"]
        assert row["bins"]["native_no_echo"] == row["bins"]["native_two_echo"]
        assert row["bins"]["native_no_echo"] == row["bins"]["permanently_cut_one_echo"]
    at_pin = next(row for row in archived["appointments"] if row["kappa"] == 8 and row["exposure"] == 37)
    assert .065 < at_pin["total_variation"] < .066
    short = archived["short_exposure_certificate"]
    retry = sum(Fraction(364, 729)**i for i in range(8))
    assert Fraction(short["difference"]) == Fraction(49, 729)*retry


def test_replay_rejects_changed_memory_rule_bins_and_nonfinite_metrics(monkeypatch):
    law = study.candidate.native_law
    echo = law.echo()
    wrong = echo.copy()
    wrong[0, 0] += .125
    monkeypatch.setattr(law, "echo", lambda: wrong)
    with pytest.raises(ArithmeticError, match="Native implementation"):
        study.local_certificate()
    archived = json.loads(study.ARCHIVE.read_text())
    for corruption in ("bin permutation", "nonfinite metric", "missing arm"):
        wrong = deepcopy(archived)
        row = wrong["appointments"][0]
        if corruption == "bin permutation":
            row["bins"]["native_one_echo"].reverse()
        elif corruption == "nonfinite metric":
            row["total_variation"] = float("nan")
        else:
            del row["bins"]["native_two_echo"]
        with pytest.raises(ArithmeticError):
            study.compare(archived, wrong)
