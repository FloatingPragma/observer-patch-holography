"""Independent regional entropy, geometric admission, and changed source controls."""
from copy import deepcopy
from decimal import Decimal as D, localcontext
from fractions import Fraction as F
import json
import math

import numpy as np
import pytest

from . import experiment
from . import finite_rods
from ..proper_balls.measure import PI
from .memories import LOG_ODDS, P, SENSITIVITY, THETA, Pair


def _density(theta):
    rotate = np.eye(4)
    c, s = math.cos(theta), math.sin(theta)
    rotate[np.ix_([0, 3], [0, 3])] = ((c, -s), (s, c))
    return (rotate*np.array([5, 1, 1, 1])/8) @ rotate.T


def _entropy(state):
    eigen = np.linalg.eigvalsh(state)
    return -sum(eigen*np.log(eigen))


def test_joint_and_restricted_matrices_independently_measure_entropy_and_modular_change():
    ref = _density(THETA)
    values, vectors = np.linalg.eigh(ref)
    K = (-vectors*np.log(values)) @ vectors.T
    for d in (-0.06, 0.0125, 0.09):
        state = _density(THETA+d)
        pair = Pair(d)
        endpoint = np.trace(state.reshape(2, 2, 2, 2), axis1=1, axis2=3)
        assert math.isclose(endpoint[1, 1], P+pair.probability_change)
        assert abs(_entropy(endpoint)-_entropy(np.diag((1-P, P)))-pair.cut_entropy_change) < 3e-15
        assert abs(_entropy(state)-_entropy(ref)) < 3e-15
        assert abs(np.trace((state-ref)@K)-pair.whole_modular_change) < 3e-15


def test_unforced_memory_rotation_at_fixed_source_rejects_the_full_operator_premise():
    d = 1e-5
    slope = (Pair(d).cut_entropy_change-Pair(-d).cut_entropy_change)/(2*d)
    assert abs(slope-SENSITIVITY*LOG_ODDS) < 1e-9
    assert 0.22 < slope < 0.23
    # The two cut outcomes of a given stock have identical free sheet stress.
    # No scalar action operator can remove both vacuum modular eigenvalues.
    K = -np.log([1-P, P])
    assert math.isclose(max(abs(K-K.mean())), LOG_ODDS/2)


def test_source_blind_radial_brackets_at_every_refinement():
    # Exact squared flat crossing radii for the seven shapes on the rational ray.
    crosses = [F(1)]
    crosses += [1/(1-F(16,25)*F(x,21)) for x in (1,4,16)]
    crosses += [1/(1-F(8,9)*F(x,42)) for x in (9,25,36)]
    ordered = sorted(crosses)
    assert all(right/left > F(161,159)**2 for left,right in zip(ordered,ordered[1:]))
    assert ordered[-1]*F(161,159)**2 < 16*ordered[0]
    assert 36*math.pi/(30*32**2) < 1/160


@pytest.fixture(scope="module")
def receipt():
    return json.loads((experiment.HERE/"results.json").read_text())


def test_finite_source_uses_the_tensor_decoder_and_both_controls(receipt):
    for field in ("ideal_tomography", "finite_tomography"):
        def arms(stock):
            return [row for row in receipt[field] if row["stock"] == stock]
        held = arms("held")
        assert [r["inverse_ell"] for r in held] == list(experiment.SCALES)
        for key in ("largest_metric_action_error","largest_entropy_action_error","largest_balance_over_ell4"):
            sequence = [r[key] for r in held]
            assert sequence[1] < sequence[0]/5 and sequence[2] < sequence[1]/5
            assert sequence[-1] < 0.001
        expected = (19/33,0,0,0,0,0,0,0,0,-19/33)
        np.testing.assert_allclose(held[-1]["from_action"], expected, atol=2e-15)
        assert arms("half_defect")[-1]["largest_balance_over_ell4"] > 5
        assert arms("blind_pulse")[-1]["largest_balance_over_ell4"] > 11
    for field in ("rows", "finite_rows"):
        rows = receipt[field]
        for r in rows:
            ell = 1/r["inverse_ell"]
            for stock in ("ordinary", "held"):
                trial = r["trials"][stock]
                assert trial["relative_entropy"] >= 0
                # Numerical coefficients are deliberately loose; the uniform
                # eighth-order bound is proved for all finite source boards.
                assert trial["relative_entropy"] < 1e5*ell**8
                assert abs(trial["vacuum_modular_change"]-trial["entropy_change"]-trial["relative_entropy"]) < 1e-18
        assert len(rows) == 21


def test_counted_holdout_predicts_one_geometric_response(receipt):
    records = receipt["counted_endpoints"]
    low, high = records["held_minus_ordinary_entropy_99"]
    assert 0 < records["opposite_half_strength_area"] < low
    assert low < records["opposite_finite_rod_area"] < high
    assert len(set(records["one_counts"].values())) == 3


def _decimal_fraction(x):
    return D(x.numerator)/D(x.denominator)


def _rest_series(K, ell):
    # The rest-frame angular average of (1-mu²)^n equals product(2i/(2i+1)).
    # Integrate the sinc series radially for area, volume and diamond weight.
    moments = [D(1)]
    for n in range(1,18):
        moments.append(moments[-1]*(-K)*D(2*n)/D(2*n+1)/D((2*n)*(2*n+1)))

    def data(r):
        terms = [c*r**(2*n) for n,c in enumerate(moments)]
        v = r**3*sum(t/D(2*n+3) for n,t in enumerate(terms))
        a = r**2*sum(terms)
        w = r**4*sum(t/D((2*n+3)*(2*n+5)) for n,t in enumerate(terms))
        return v,a,w

    r = ell
    for _ in range(7):
        v,a,_ = data(r)
        r -= (v-ell**3/3)/a
    _,a,w = data(r)
    return 4*PI*(a-ell**2), 4*PI*w


def _decimal_cut_entropy(work):
    p = D(3)/8
    angle = work*(2*PI)/(D(5)/3).ln()/(D(3).sqrt()/4)
    x = 2*angle
    sine = sum((-1)**n*x**(2*n+1)/D(math.factorial(2*n+1)) for n in range(12))
    cosine = sum((-1)**n*x**(2*n)/D(math.factorial(2*n)) for n in range(12))
    shifted = p + D(3).sqrt()*sine/8 + (1-cosine)/8
    def h(q):
        return -q*q.ln()-(1-q)*(1-q).ln()
    return h(shifted)-h(p)


def test_normalized_rest_residuals_against_high_precision_entropy_and_integrals(receipt):
    with localcontext() as ctx:
        ctx.prec = 65
        for inverse in (32,512):
            ell = D(1)/inverse
            reference = F(1,inverse)
            for stock, rho in (("ordinary",F(11,33)), ("held",F(19,33))):
                density = _decimal_fraction(rho)
                K = 2*PI*density
                flat = finite_rods.read_ball(0, reference, (1,0,0,0))
                rod = finite_rods.read_ball(finite_rods.angle_gain(reference)*rho, reference, (1,0,0,0))
                inputs = {"rows": _rest_series(K,ell), "finite_rows": (
                    4*PI*_decimal_fraction(rod.area_over_4pi-flat.area_over_4pi),
                    4*PI*_decimal_fraction(rod.weighted_over_4pi))}
                for field,(area,w) in inputs.items():
                    row, = [r for r in receipt[field] if r["shape"] == "rest" and r["inverse_ell"] == inverse]
                    t = row["trials"][stock]
                    entropy = _decimal_cut_entropy(density*w)
                    bulk = entropy/ell**4-8*PI**2*density/15
                    assert math.isclose(t["bulk_over_ell4"], float(bulk), rel_tol=1e-12, abs_tol=1e-20)
                    assert abs(t["geometric_over_ell4"]-float(area/ell**4+4*PI*K/15)) < 2e-14
                    assert abs(t["balance_over_ell4"]-float((area+entropy)/ell**4)) < 2e-14


@pytest.mark.parametrize("field", ["bulk_over_ell4", "geometric_over_ell4", "balance_over_ell4"])
def test_small_residual_replay_detects_changed_fine_scale_ordinary_source(receipt, field):
    row, = [r for r in receipt["rows"] if r["shape"] == "rest" and r["inverse_ell"] == 512]
    changed = deepcopy(row)
    changed["trials"]["ordinary"][field] += 1e-10
    with pytest.raises(ValueError, match=field):
        experiment._compare_boards([row], [changed])
