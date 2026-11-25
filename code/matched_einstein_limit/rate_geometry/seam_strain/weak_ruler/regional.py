"""Entropy of both seam heads, bridge strain and all enclosed spent-pump modes."""

import argparse
from fractions import Fraction
from functools import lru_cache
import json
from pathlib import Path

import numpy as np
import sympy as sp

from . import aperture, threshold
from .. import live_modular, live_source
from ..source import certificate as cert, measure as source

HERE = Path(__file__).parent


def select(wave, ports):
    # Restore the port axis without changing the exact fixed-point limbs.
    return wave.reshape(2, -1, 4, 16)[:, ports].reshape(-1, 16)


def square_root(inside_hi, inside_lo, colerr):
    """Rational center of sqrt(V_in† V_in), with a verified operator error."""
    wave = inside_hi.astype(object) * (1 << 36) + inside_lo
    re, im = np.split(wave, 2, axis=1)
    gram = (
        sp.Matrix(re.T @ re + im.T @ im) + sp.I * sp.Matrix(re.T @ im - im.T @ re)
    ) / (1 << 164)
    vals, vec = np.linalg.eigh(np.array(gram, dtype=complex))
    proposal = (vec * np.sqrt(vals)) @ vec.conj().T
    root = cert.prior.sym_matrix((proposal + proposal.conj().T) / 2)
    _, w, _, radius = live_modular.spectral_center(root)
    assert min(w.diagonal()) - radius > sp.Rational(9, 10)
    # The isometry has eight columns. Restricting output ports keeps the
    # operator error below sqrt(8)*colerr and the true norm at most one.
    gram_err = 6 * colerr + 9 * colerr**2
    residual = cert.prior.one_norm((gram - root * root).applyfunc(sp.expand)) + gram_err
    # The Sylvester equation for the difference of positive square roots.
    return root, residual * sp.Rational(10, 9)


def logarithm(reference, err):
    center, rounding = live_modular.rounded(reference)
    v, w, defect, radius = live_modular.spectral_center(center)
    uncertainty = radius + err + rounding
    lower = min(w.diagonal()) - uncertainty
    assert lower > sp.Rational(1, 10**8)
    logs = [sp.Rational(float(np.log(float(x)))) for x in w.diagonal()]
    for x, y in zip(w.diagonal(), logs, strict=True):
        cert.prior.log_bound(x, y)
    log = (-v * sp.diag(*logs) * v.H).applyfunc(sp.expand)
    logerr = (
        uncertainty / lower
        + sp.Rational(1, 10**12)
        + 3 * defect * max(abs(x) for x in logs)
    )
    return log, logerr, lower


@lru_cache(maxsize=2)
def certify(radius):
    board = source.formation.federation(1)
    ports = aperture.enroll(board, radius)[0]
    if not set(source.formation.region(board)) <= set(ports):
        raise ValueError(
            "the region must contain the physical supports of both heads and the bridge"
        )
    exterior = np.setdiff1d(np.arange(board.n), ports)
    hi, lo, colerr, _ = live_source.pump_center()
    j, jerr = square_root(select(hi, ports), select(lo, ports), colerr)
    outside = live_source.rail_grams(select(hi, exterior), select(lo, exterior))
    vacuum_err = 6 * colerr + 9 * colerr**2
    (formed, product), ferr, _ = live_source.formation_center()

    def compress(head, weights):
        h = sp.kronecker_product(sp.diag(*weights) / 3, head)
        return (
            (j * h * j).applyfunc(sp.expand),
            live_source.rational_read(outside, head, weights),
        )

    ref = compress(formed, (2, 1))
    channel_err = (2 * jerr + jerr**2, vacuum_err)
    logs = [
        logarithm(block, ferr + (1 + ferr) * chan)
        for block, chan in zip(ref, channel_err, strict=True)
    ]
    rotations = live_modular.rotation_centers(formed, ferr)
    signals = {}
    for tag, head, weights, inp_err, norm in [
        ("ordinary", formed, (-1, 1), ferr, 1),
        ("decorrelated", product - formed, (2, 1), 3 * ferr, 2),
        *[(k, h, (2, 1), e, 1) for k, (h, e) in rotations.items()],
        *[
            (a + b, h, (2, 1), 0, 1)
            for (a, b), h in zip(
                ((a, b) for a in "XYZ" for b in "XYZ"),
                live_source.qubit_hidden.correlation_axes(),
                strict=True,
            )
        ],
    ]:
        blocks = compress(head, weights)
        mean = error = 0
        for block, chan, (log, logerr, _) in zip(
            blocks, channel_err, logs, strict=True
        ):
            b, round_err = live_modular.rounded(block)
            err = inp_err + chan * (norm + inp_err) + round_err
            # Each arm of the regional channel contracts trace norm on
            # Hermitian inputs. True ||-log(reference block)|| < 19.
            mean += sp.re(sp.trace(log * b))
            error += norm * logerr + (19 + logerr) * err
        signals[tag] = cert.enclosed(mean, error)
    dz = threshold.z_intervals()
    alpha = [
        -Fraction(x) / Fraction(y) for x in dz["decorrelated"] for y in dz["right"]
    ]
    matched = [
        Fraction(x) + a * Fraction(y)
        for x in signals["decorrelated"]
        for a in alpha
        for y in signals["right"]
    ]
    assert min(matched) > Fraction(19, 100)
    return dict(
        radius=radius,
        ports=ports,
        occupied_and_vacuum_trace=[float(sp.re(sp.trace(x))) for x in ref],
        minimum_reference_eigenvalue_lower=[cert.enclosed(x[2], 0)[0] for x in logs],
        log_operator_error_upper=[cert.enclosed(x[1], 0)[1] for x in logs],
        regional_entropy_slopes=signals,
        fixed_metric_pairing_slope=cert.enclosed(
            (max(matched) + min(matched)) / 2, (max(matched) - min(matched)) / 2
        ),
    )


def census_transfers(row, slopes):
    from mpmath import iv

    # Use the independently enclosed numerator and the enclosed ordinary
    # boundary denominator. Every measured change includes its full interval.
    slopes = {k: iv.mpf(v) for k, v in slopes.items()}
    result = {}
    for side, side_row in row["certificate"]["fixed_census"].items():
        config = side_row["matched"]

        def matched_k(tag):
            return slopes[tag] + iv.mpf(config[tag]["angle"]) * slopes[side]

        for orientation in ("outward", "inward"):
            ordinary = iv.mpf(config["ordinary"][orientation])
            if ordinary.a <= 0 <= ordinary.b:
                continue
            gain = -matched_k("ordinary") / ordinary
            gap = matched_k("decorrelated") + gain * iv.mpf(
                config["decorrelated"][orientation]
            )
            result[side + "_" + orientation] = dict(
                ordinary_gain=aperture.encode(gain), held_residual=aperture.encode(gap)
            )
    return result


def lottery_response(slopes, fronts):
    from mpmath.ctx_iv import MPIntervalContext

    ctx = MPIntervalContext()
    ctx.dps = 40
    scale = fronts["native_lotteries"]["common_scale"]
    margin = fronts["maximum_coordinate_correction_upper"]
    total = ctx.mpf(0)
    for i, (name, v) in enumerate(
        zip(fronts["basis"], fronts["blind_direction_center"], strict=True)
    ):
        # The stored centers are dyadic numbers. Include a binary64 ulp for
        # their serialization as well as the certified true coordinate shift.
        delta = ctx.mpf([-1, 1]) * (
            np.spacing(abs(v)) + (margin if i in fronts["correction_columns"] else 0)
        )
        total += ctx.mpf(slopes[name]) * (ctx.mpf(v) + delta)
    return aperture.encode(total / scale)


def match_all_three(aperture_row, region):
    """Transfer the detector's exact three-reading kernel to the enlarged matter algebra."""
    from mpmath.ctx_iv import MPIntervalContext

    ctx = MPIntervalContext()
    ctx.dps = 40
    kernel = aperture_row["three_equal_readings"]
    slopes = region["regional_entropy_slopes"]
    change = ctx.mpf(slopes["decorrelated"])
    for k, v in kernel["control_slopes"].items():
        change += ctx.mpf(v) * ctx.mpf(slopes[k])
    # In the positive lottery the unnormalized decorrelated arm has weight
    # one. Its reported probability is exactly the reciprocal normalization.
    weight = kernel["finite_native_lotteries"]["positive"]["decorrelated"]
    gap = change * ctx.mpf(weight)
    t = ctx.mpf(2) ** -42
    lam = min(region["minimum_reference_eigenvalue_lower"])
    # On the common support each mixed reference has eigenvalues >= (1-t)lam.
    # The log resolvent bounds |S''| by 4/((1-t)lam) for either arm.
    remainder = 4 * t * t / ((1 - t) * lam)
    interval = gap * t + ctx.mpf([-1, 1]) * remainder
    assert interval.a > 0
    return dict(
        fixed_readings=kernel["fixed"],
        regional_entropy_slope=aperture.encode(change),
        lottery_modular_gap=aperture.encode(gap),
        common_admixture=float(t),
        finite_entropy_difference=aperture.encode(interval),
    )


def assemble(panels):
    record = json.loads((HERE / "aperture.json").read_text())
    fronts = json.loads((HERE / "observability.json").read_text())
    return dict(
        regions=panels,
        fixed_census_and_both_currents={
            q: match_all_three(rows[1], panels[0])
            for q, rows in record["panels"].items()
        },
        matched_front_native_lottery_entropy_slopes=[
            lottery_response(row["regional_entropy_slopes"], fronts) for row in panels
        ],
        fixed_census_transfers={
            q: census_transfers(rows[1], panels[0]["regional_entropy_slopes"])
            for q, rows in record["panels"].items()
        },
    )


def run():
    return assemble([certify(r) for r in (1, 2)])


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    result = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if p.parse_args().write:
        (HERE / "regional.json").write_text(result)
    else:
        print(result, end="")
