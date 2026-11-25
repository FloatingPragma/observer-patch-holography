"""Curvature and fixed-volume tests on the six-station native radar survey."""

from itertools import combinations, permutations
from pathlib import Path
import argparse
import json

import numpy as np
import sympy as sp
from scipy.optimize import brentq

from ..full_algebra import certificate as cert
from ..full_algebra import native_precision as native
from ..null_optics.certificate import bounded_sources

HERE = Path(__file__).resolve().parent
iv = cert.iv


def surveyed_ranges():
    """Enclose all six-station mean ranges on the certified equilibrium."""
    board = native.iso.pairing.formation.federation(1)
    source = native.iso.pairing.formation.region(board)
    ports = (*source, 5, 6)
    columns = np.array([4 * p + i for p in ports for i in range(4)])
    weights = native.prior.weights()
    blocks = []
    for odd in (False, True):
        pooled = 0
        for sign in (-1, 1):
            h = native.generator(board, source, sign, odd)
            panels = native.kernel_panels(h, columns, weights)
            pooled += np.array(native.prior.direction_effects(panels, 12), dtype=object)
        blocks.append(pooled)
    effects = cert.rational_effects(np.stack(blocks, axis=1), 8 * native.prior.FRAC)
    center, state_error = bounded_sources(sp.eye(4) / 4)["state"]
    error = 3 * state_error + sp.Rational(native.prior.error_effect())
    means = np.full((6, 6), iv.mpf(0), dtype=object)
    for (i, j), effect in zip(permutations(range(6), 2), effects, strict=True):
        mean = sum(cert.inner(h, x) for h, x in zip(effect, center))
        means[i, j] = cert.ri_band(mean, error)
    return ports, means


def squared_distances(means):
    """Keep the original five clock gains fixed when enrolling the held receiver."""
    log_ratios = np.full((6, 6), iv.mpf(0), dtype=object)
    for i, j in permutations(range(6), 2):
        log_ratios[i, j] = iv.ln(means[i, j] / means[j, i])
    gains = [sum(log_ratios[i, :5]) / 5 for i in range(5)]
    gains.append(sum(log_ratios[5, i] + gains[i] for i in range(5)) / 5)
    squared = np.full((6, 6), iv.mpf(0), dtype=object)
    for i, j in permutations(range(6), 2):
        squared[i, j] = means[i, j] * means[j, i] * iv.exp(-gains[i] - gains[j])
    return squared


def eigenvalue_floor(squared, floor):
    """Sylvester test on 1-perp with a row-sum allowance for interval errors."""
    n = len(squared)
    scale = 10**8
    centers = sp.Matrix(
        n, n, lambda i, j: sp.Rational(round(float(squared[i, j].mid) * scale), scale)
    )
    err = cert.ri("0.00010")
    assert all(
        abs(v - cert.ri(c)).b < err for v, c in zip(squared.flat, centers, strict=True)
    )
    contrast = sp.eye(n)[:, :-1] - sp.Matrix(n, n - 1, lambda i, j: int(i == n - 1))
    test = (
        contrast.T
        * (-centers / 2 - (floor + sp.Rational(n - 1, 20000)) * sp.eye(n))
        * contrast
    )
    assert all(test[:k, :k].det() > 0 for k in range(1, n))


def det_interval(matrix):
    n = len(matrix)
    value = iv.mpf(0)
    for order in permutations(range(n)):
        odd = sum(a > b for i, a in enumerate(order) for b in order[i + 1 :]) % 2
        term = iv.mpf(-1 if odd else 1)
        for i, j in enumerate(order):
            term *= matrix[i, j]
        value += term
    return value


def spherical_test(squared):
    dist = np.array([[float(iv.sqrt(v).mid) for v in row] for row in squared])
    speed = brentq(
        lambda z: np.linalg.eigvalsh(np.cos(z * dist[:5, :5]))[0], 0.45, 0.78
    )
    # Frozen, outward brackets of radius 1e-4 about the numerical inverse radius.
    scale = 10000
    tick = int(speed * scale)
    lo, hi = (
        cert.ri(sp.Rational(tick - 1, scale)),
        cert.ri(sp.Rational(tick + 2, scale)),
    )

    def gram(z):
        return np.array([[iv.cos(z * iv.sqrt(v)) for v in row] for row in squared])

    lower, upper = gram(lo)[:5, :5], gram(hi)[:5, :5]
    assert det_interval(lower).a > 0 and det_interval(upper).b < 0
    z = iv.mpf([lo.a, hi.b])
    both = gram(z)
    assert max(z * iv.sqrt(v) for v in squared[:5, :5].flat) < iv.pi
    assert max(both[i, j].b for i, j in permutations(range(5), 2)) < 0
    # At a determinant zero the 4x4 positive leading minor gives a rank-4 PSD Gram.
    assert all(det_interval(both[:j, :j]).a > 0 for j in range(1, 5))
    fail = both.sum() / 6
    assert fail.b < 0
    # For smaller positive curvature the centered six-station Gram has five
    # positive directions. For larger curvature its uniform direction is negative.
    anchor_sq = cert.ri(9)
    max_row_error = max(
        sum(abs(row[j] - anchor_sq) for j in range(6) if j != i).b
        for i, row in enumerate(squared)
    )
    alpha_cut = cert.ri("0.65")
    low_margin = anchor_sq / 2 - alpha_cut**2 * anchor_sq**2 / 24 - max_row_error / 2
    high_rayleigh = gram(alpha_cut).sum() / 6
    assert low_margin.a > 0 and high_rayleigh.b < 0
    return dict(
        inverse_radius_bracket=[(tick - 1) / scale, (tick + 2) / scale],
        all_positive_curvatures=dict(
            inverse_radius_split=0.65,
            centered_lower_coefficient=str(low_margin),
            uniform_upper_at_split=str(high_rayleigh),
        ),
        curvature_bracket=str(z * z),
        calibrated_sixth_rayleigh=str(fail),
    )


def hyperbolic_test(squared):
    """Exclude every negative space-form curvature for these six minimizing ranges."""
    distances = np.array([[iv.sqrt(v) for v in row] for row in squared])
    # In the Poincare ball, sinh(mu*d/2) equals a Euclidean chord times
    # endpoint factors. Those factors cancel from all three Ptolemy products.
    best_gap = 0
    for sites in combinations(range(6), 4):
        i, j, k, last = sites
        pairings = [((i, j), (k, last)), ((i, k), (j, last)), ((i, last), (j, k))]
        bounds = [sum(distances[e] for e in edges) for edges in pairings]
        top = max(range(3), key=lambda n: float(bounds[n].mid))
        gap = min((bounds[top] - v).a for n, v in enumerate(bounds) if n != top)
        if gap > best_gap:
            best_gap, chosen, longest = gap, sites, pairings[top]
    lower = min(distances[e].a for e in longest)
    large = 1
    while (1 - iv.exp(-large * lower)) ** 2 * iv.exp(large * best_gap / 2) <= 2:
        large += 1
    assert (1 - iv.exp(-large * lower)) ** 2 * iv.exp(large * best_gap / 2) > 2

    # Near zero the five-point tangent Gram lower bound strictly dominates
    # the Taylor remainder, giving four negative directions, one too many.
    small = cert.ri("0.2")
    upper_d2 = cert.ri("9.03")
    remainder = 4 * small**2 * upper_d2**2 * iv.exp(small * iv.sqrt(upper_d2)) / 24
    assert remainder.b < cert.ri("3.49")

    charts = list(combinations(range(6), 5))
    count = 0
    min_width = 1

    def cover(start, stop):
        nonlocal count, min_width
        mu = iv.mpf([cert.ri(start).a, cert.ri(stop).b])
        # Multiply the whole Gram by exp(-3*mu) to keep cancellation controlled.
        g = np.array(
            [
                [(iv.exp(mu * (v - 3)) + iv.exp(-mu * (v + 3))) / 2 for v in row]
                for row in distances
            ]
        )
        for sites in charts:
            minor = g[np.ix_(sites, sites)]
            midpoints = np.array([[float(x.mid) for x in row] for row in minor])
            try:
                inverse = np.linalg.inv(midpoints)
            except np.linalg.LinAlgError:
                continue
            x = np.array([[iv.mpf(float(v)) for v in row] for row in inverse])
            defect = np.eye(5, dtype=int) - x @ minor
            if max(sum(abs(v) for v in row).b for row in defect) < 1:
                count += 1
                min_width = min(min_width, stop - start)
                return
        if stop - start < sp.Rational(1, 20000):
            raise ArithmeticError(
                f"Unresolved hyperbolic Gram interval {start}, {stop}"
            )
        mid = (start + stop) / 2
        cover(start, mid)
        cover(mid, stop)

    cover(sp.Rational(1, 5), sp.Integer(large))
    return dict(
        small_inverse_radius_max=0.2,
        certified_subintervals=count,
        minimum_subinterval_width=str(min_width),
        large_inverse_radius_min=large,
        ptolemy_station_indices=list(chosen),
        ptolemy_largest_pairing_indices=longest,
        pairing_sum_gap_lower=float(best_gap.a) - 1e-12,
    )


def volume_interventions():
    """Repeat one physical decorrelation with three entropy-preserving switch settings."""
    from .. import signed_transfers as transfers
    from ..full_algebra import station
    from ..dressed_charge.measure import read_source

    state, _ = read_source()
    r = transfers.right_coin(state)
    left = np.trace(state.reshape(2, 2, 2, 2, 2), axis1=2, axis2=4)
    product = np.stack([np.kron(block, r) for block in left])
    axes = transfers.phase.AXES[1:]
    bloch = np.array([np.trace(r @ a).real for a in axes])
    axis = sum(v * a for v, a in zip(bloch, axes)) / np.linalg.norm(bloch)
    effects = station.radar_menus()[6]
    baseline_area, baseline_volume = station.geometry(station.pair(effects, state))
    t = 1e-7
    mixed = state + t * (product - state)

    def survey(angle):
        u = transfers.rotation(axis, angle)
        rotated = u @ mixed @ u.conj().T
        area, volume = station.geometry(station.pair(effects, rotated))
        return rotated, area - baseline_area, volume - baseline_volume

    match = brentq(lambda angle: survey(angle)[2], -0.00006, 0.0, xtol=1e-14)
    rows = []
    for angle in (0.0, match, 2 * match):
        x, da, dv = survey(angle)
        ds = station.entropy(x) - station.entropy(state)
        assert np.max(abs(transfers.right_coin(x) - r)) < 1e-14
        assert (
            np.max(abs(np.trace(x.reshape(2, 2, 2, 2, 2), axis1=2, axis2=4) - left))
            < 1e-14
        )
        assert abs(np.trace((x - state)[1])) < 1e-14
        rows.append(
            dict(angle=angle, entropy_change=ds, area_change=da, volume_change=dv)
        )
    assert rows[0]["volume_change"] < 0 < rows[2]["volume_change"]
    return dict(unpairing_probability=t, numeric_center=rows)


def source_volume():
    """Exact local freedom of surveyed volume at the same energy and entropy."""
    signed = json.loads((HERE.parent / "full_algebra/certificate.json").read_text())[
        "signed_transfers"
    ]

    def enclosed(key):
        lo, hi = signed[key].strip("[]").split(",")
        return iv.mpf([lo.strip(), hi.strip()])

    det = (
        enclosed("ordinary_energy")
        * enclosed("pairing_entropy")
        * enclosed("right_adjuster_volume")
    )
    assert det.b < cert.ri("-0.0000017")
    return str(det)


def run():
    iv.dps = 32
    ports, means = surveyed_ranges()
    squared = squared_distances(means)
    eigenvalue_floor(squared[:5, :5], sp.Rational("3.49"))
    assert max(v.b for v in squared[:5, :5].flat) < cert.ri("9.03")
    # If |sectional curvature| R^2 <= 1/2, sinc >= 11/12 and sinhc <= 12/11.
    error = sp.Rational(903, 100) * (sp.Rational(12, 11) ** 2 - 1)
    assert 2 * error < sp.Rational("3.49")
    return dict(
        cells=80,
        stations=list(map(int, ports)),
        deadline=6,
        directed_mean_range_bounds=[[str(v) for v in row] for row in means],
        projected_squared_distance_bounds=[[str(v) for v in row] for row in squared],
        reference_centered_gram_floor=3.49,
        max_squared_reference_distance=9.03,
        max_flat_squared_error_if_curvature_radius_squared_at_most_half=str(error),
        required_flat_squared_error_lower=1.745,
        positive_space_form=spherical_test(squared),
        negative_space_form=hyperbolic_test(squared),
        volume_interventions=volume_interventions(),
        energy_entropy_volume_jacobian=source_volume(),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    result = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if parser.parse_args().write:
        (HERE / "curved.json").write_text(result)
    else:
        print(result, end="")
