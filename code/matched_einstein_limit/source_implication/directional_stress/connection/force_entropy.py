"""Different entropy with all twelve mean ranges and both native marked forces fixed."""

import argparse
import json

import numpy as np
import sympy as sp

from ...full_algebra import certificate as c
from ...full_algebra.reuse import invisible
from . import equal_radar
from .force_gradients import force_gradients

HERE = equal_radar.HERE


def range_gradients(basis):
    nums, unit, error = c.exact_effects()
    pairs = c.rational_effects(nums, unit)
    panel = sp.Matrix(
        [[c.inner(block, a) for block in pair for a in basis] for pair in pairs]
    )
    return panel, sp.Rational(error)


def corrected_kernel_direction(ranges, range_error, forces, force_error, *, score=None, proposal=None):
    """Enclose an exact tangent holding all twelve ranges and both forces fixed."""
    a = ranges.col_join(forces)
    assert a.shape == (14, 18)
    if proposal is None:
        row_array = np.array(a, dtype=float)
        scales = np.linalg.norm(row_array, axis=1)
        scaled = row_array / scales[:, None]
        _, _, right = np.linalg.svd(scaled, full_matrices=True)
        v = right[14:].T @ (right[14:] @ np.array(score, dtype=float).ravel())
        d = invisible.dyadic((v / np.linalg.norm(v))[:, None])
        k = invisible.dyadic(np.linalg.pinv(scaled) / scales[None, :])
    else:
        d = sp.Matrix(proposal)
        # Fixed rational control proposals use a platform-independent right inverse.
        right_inverse = a.T * (a * a.T).inv(method="DM")
        k = right_inverse.applyfunc(lambda x: sp.Rational(round(x * 2**38), 2**38))
    d -= k * (a * d)
    n = invisible.frobenius
    # Rounding and physical range propagation give ||Rhat-R|| <30*range_error.
    # The two force rows each have Euclidean error <=force_error.
    op_err = c.iv.sqrt(c.ri((30 * range_error) ** 2 + 2 * force_error**2))
    defect = n(sp.eye(14) - a * k) + op_err * n(k)
    assert defect.b < c.ri(".000002")
    change = n(k) * (n(a * d) + op_err * n(d)) / (1 - defect)
    assert change.b < c.ri(".000002")
    return d, change, defect, op_err


def certify(ranges, range_error, forces, force_error):
    c.iv.dps = 32
    source, _, _, corr, _ = equal_radar.tangent_source()
    frame = [c.prior.sym_matrix(x) for x in c.iso.basis(False)]
    info = source["fixed_state_centers"]
    cond_radius = sp.Rational(source["conditional_HS_error_upper"])
    g, conditional_scores, centers = [], [], []
    lower = (sp.Rational(".20"), sp.Rational(".21"))
    for lo, row in zip(lower, info["even_then_odd"]):
        state = sp.eye(4) / 4 + sum(
            (sp.Rational(v, info["denominator"]) * b for v, b in zip(row, frame)),
            sp.zeros(4),
        )
        assert c.global_survey.positive_definite(state - (lo + cond_radius) * sp.eye(4))
        centers.append(state)
        grad, log_err = c.prior.entropy_coordinates(state, 0, lo, corr)
        g.extend(grad)
        # Nine orthonormal coordinates. Each nominal log coefficient has error
        # <=||B_i||_1*log_err=2*log_err. Projection of the true log-state shift
        # has L2 norm <= radius/lambda_min by the resolvent bound.
        conditional_scores.append(6 * log_err + cond_radius / lo)
    g = sp.Matrix(g)
    d, change, defect, op_err = corrected_kernel_direction(
        ranges, range_error, forces, force_error, score=g
    )
    n = invisible.frobenius
    sizes = [n(d[i : i + 9, :]) for i in (0, 9)]
    norm = n(d) + change
    assert norm.b < c.ri("1.00001")
    score = c.ri(g.dot(d))
    err = n(g) * change + sum(
        c.ri(e) * (size + change) for e, size in zip(conditional_scores, sizes)
    )
    slope = score + c.iv.mpf([-1, 1]) * err
    assert slope.a > c.ri(".000009"), str(slope)
    t = sp.Rational("1e-7")
    # Each conditional center has spectrum above the independently verified
    # lower bounds above. Minimum parity weight .134 includes both branches.
    eig = c.ri(".134") * c.ri(min(lower)) - c.ri(t) * norm
    assert eig.a > c.ri(".026")
    gain = c.ri(t) * slope.a - c.ri(t * t) * norm.b**2 / (2 * c.ri(".026"))
    assert gain.a > c.ri("7e-13")
    # The independent source certificate gives a >.20 Pauli-lottery slack for
    # each sector at every allowed weight. A perturbation of any correlation
    # coordinate vector costs at most 12 |t| ||D||_HS in either lottery.
    slack = c.ri(source["lottery_slack_lower"]) - 12 * c.ri(".01") * norm
    assert slack.a > c.ri(".079")
    labels = ("XX", "XY", "XZ", "YX", "YY", "YZ", "ZX", "ZY", "ZZ")
    probe = max(range(9), key=lambda j: abs(d[j] + d[j + 9]))
    response = (
        c.ri(2 * (d[probe] + d[probe + 9]))
        + c.iv.mpf([-1, 1]) * 2 * c.iv.sqrt(2) * change
    )
    assert 0 not in response
    # Direct positive binary measurement on the midpoint native arms, with an
    # independently chosen stationary-center weight inside its certified bracket.
    p = sp.Rational(".1343")
    parity_weight = (1 - p, p)
    pop = sum(
        w * c.inner(2 * corr[probe], state) for w, state in zip(parity_weight, centers)
    )
    magnitude = 2 * (d[probe] + d[9 + probe])
    rates = [
        float((1 + pop + x * magnitude) / 2)
        for x in (-sp.Rational(".01"), sp.Rational(".01"))
    ]
    nshots = 2**22
    counts = np.random.default_rng(251126).binomial(nshots, rates)
    observed = (int(counts[1]) - int(counts[0])) / nshots
    halfwidth = 2 * np.sqrt(np.log(4 / 0.01) / (2 * nshots))
    assert (
        abs(observed - (rates[1] - rates[0])) < halfwidth and observed + halfwidth < 0
    )
    assert float((c.ri(".01") * response).b) < -2 * halfwidth
    denom = sp.ilcm(*(v.q for v in forces))
    return dict(
        affine_readout_count=14,
        rank=14,
        correlation_slice_dimension=18,
        exact_kernel_dimension=4,
        matrix_error_norm_upper=str(op_err.b),
        right_inverse_defect_upper=str(defect.b),
        true_kernel_correction_HS_upper=str(change.b),
        nominal_entropy_slope_nats=str(score),
        signed_entropy_slope_nats=str(slope),
        native_feasible_segment=["-.01", ".01"],
        lottery_slack_lower_at_endpoints=str(slack.a),
        finite_step=str(t),
        entropy_increase_nats_lower=str(gain.a),
        independent_joint_read=dict(
            product=labels[probe],
            slope=str(response),
            plus_rates=rates,
            trials_per_endpoint=nshots,
            plus_counts=counts.tolist(),
            observed_plus_difference=observed,
            simultaneous_99pct_halfwidth=halfwidth,
        ),
        direction_center_even_then_odd=[float(v) for v in d],
        force_gradient_denominator=int(denom),
        force_gradient_numerators=[
            [int(v * denom) for v in forces[i, :]] for i in range(2)
        ],
        force_gradient_error_per_row=str(force_error),
    )


def run():
    corr = [c.prior.sym_matrix(x) for x in c.iso.basis(True)]
    return certify(*range_gradients(corr), *force_gradients())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (HERE / "force_entropy.json").write_text(data)
    else:
        print(data, end="")
