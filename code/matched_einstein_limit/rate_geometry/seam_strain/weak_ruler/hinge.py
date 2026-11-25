"""Intrinsic hinge curvature of four Euclidean tetrahedra around a radar apex.

Four outer stations form the labeled boundary and the fifth is the apex of
the filled star used in live_survey. Each radial edge lies in three tetrahedra.
Parallel-transport angle on a normal cross section gives its deficit, 2 pi
minus their three interior dihedral angles. In three dimensions the
distributional *scalar* curvature integral is twice the length-weighted sum
of interior deficits. Boundary dihedral bending has a different support.
This convention follows Regge (1961), doi:10.1007/BF02733251.

The apex link is a triangulated two-sphere with four spherical triangles.
Their spherical excess gives Omega, and summing their angles gives the
identity sum(delta_i)=4 pi-Omega. Before the outer faces, every radial hinge
has length r inside the ball and V=Omega r^3/3. Thus I_R=2r(4 pi-Omega),
and the mean scalar curvature times r^2 is 6(4 pi-Omega)/Omega. A nonzero
coefficient makes this reconstruction singular on approach to the apex.
When the intrinsic cone radius adjusts to keep that volume fixed,
dI_R/ds=-4(2 pi+Omega)/(Omega r) times dA/ds.

Four native controls are the rail probability, partial re-pairing and two
marginal rotations. Two implicit-function tests hold the whole filled volume
and integrated bulk scalar exactly fixed. The third fixed observable is the
exterior area in one test and the boundary bending in the other. A nonzero
control minor and a negative re-pairing slope give admissible one-sided
curves for negative left rotations. At these faithful references the entropy
derivative equals the independently measured modular derivative. It is
nonzero for the retained heads and strain, and on the area-fixed curve also
for the twelve-port material algebra. Identifying either algebra with the
matter in a proper ball requires a separate physical measurement. These
intrinsic spatial comparisons do not specify a spacetime null direction.
"""

import argparse
from fractions import Fraction
from itertools import combinations
from functools import lru_cache
import json
from math import nextafter, inf, comb, isqrt, factorial
from pathlib import Path

import numpy as np

from .. import live_survey
from . import flight

HERE = Path(__file__).parent
ARMS = ("ordinary", "decorrelated", "hidden")
ROTATIONS = ("left", "right")
EDGES = tuple(combinations(range(5), 2))
BASES = tuple(combinations(range(4), 3))
METRICS = (
    "bulk_scalar",
    "filled_volume",
    "spherical_link",
    "boundary_bending",
    "outer_area",
)
FIRST_FIVE = (0, 1, 2, 3, 4)


def bend(d):
    """Interior deficits, exterior boundary angles, and integrated quantities."""
    angles = {e: 0 for e in EDGES}
    for triangle in BASES:
        tet = (*triangle, 4)
        for a, b in combinations(tet, 2):
            c, k = (v for v in tet if v not in (a, b))

            # Project both other vertices onto the plane normal to the hinge.
            # Their included angle is the interior dihedral at this edge.
            def q(i, j):
                return (d[a, i] ** 2 + d[a, j] ** 2 - d[i, j] ** 2) / 2

            base = d[a, b] ** 2
            x, y = q(b, c), q(b, k)
            nn = (d[a, c] ** 2 - x * x / base) * (d[a, k] ** 2 - y * y / base)
            cs = (q(c, k) - x * y / base) / np.sqrt(nn)
            if nn.real <= 0 or abs(cs.real) >= 1:
                raise ValueError("degenerate tetrahedral hinge")
            angles[(a, b)] += np.arccos(cs)
    deficits = np.array([2 * np.pi - angles[(i, 4)] for i in range(4)])
    outer = tuple(combinations(range(4), 2))
    boundary = np.array([np.pi - angles[e] for e in outer])
    area, volume, omega = live_survey.star(d)
    return (
        np.array(
            [
                2 * d[4, :4] @ deficits,
                volume,
                omega,
                2 * np.array([d[e] for e in outer]) @ boundary,
                area,
            ]
        ),
        deficits,
        boundary,
    )


def derive(d, slope):
    return bend(d + 1e-19j * slope)[0].imag / 1e-19


def transfers(changes, k):
    out = {}
    for side in ROTATIONS:
        vr = changes[side][1]
        if abs(vr) < 1e-15:
            raise ValueError("cannot match volume on this rotation")
        records = {}
        for arm in ARMS:
            alpha = -changes[arm][1] / vr
            ds = changes[arm] + alpha * changes[side]
            records[arm] = dict(
                rotation=float(alpha),
                bulk_scalar=float(ds[0]),
                filled_volume=float(ds[1]),
                boundary_bending=float(ds[3]),
                solid_angle=float(ds[2]),
                outer_area=float(ds[4]),
                modular=float(k[arm] + alpha * k[side]),
            )
        ordinary = records["ordinary"]
        # Compare dI = gain*d<K_ref>, both at the independently held volume.
        gain = ordinary["bulk_scalar"] / ordinary["modular"]
        held = {}
        for arm in ARMS[1:]:
            row = records[arm]
            pred = gain * row["modular"]
            held[arm] = dict(
                predicted_bulk_scalar=pred,
                measured_minus_predicted=row["bulk_scalar"] - pred,
                geometry_ratio=row["bulk_scalar"] / ordinary["bulk_scalar"],
                modular_ratio=row["modular"] / ordinary["modular"],
            )
        out[side] = dict(ordinary_gain=gain, sources=records, held=held)
    return out


def panel(record, modular, divisor, continuous):
    rows = flight.pooled(record, divisor, continuous=continuous)
    ix = np.ix_(FIRST_FIVE, FIRST_FIVE)
    d = rows["reference"][ix]
    vals, deficits, boundary = bend(d)
    delta = {k: derive(d, rows[k][ix]) for k in (*ARMS, *ROTATIONS)}
    # Independent higher-order symmetric real finite differences check sign,
    # scale and the complex derivative on every source direction.
    fd_errors = {}
    for k, tangent in delta.items():
        step = 0.001 / max(abs(rows[k][ix]).ravel())
        v = step * rows[k][ix]
        fd = (
            8 * (bend(d + v)[0] - bend(d - v)[0])
            - (bend(d + 2 * v)[0] - bend(d - 2 * v)[0])
        ) / (12 * step)
        fd_errors[k] = float(max(abs(tangent - fd)))
    omega = vals[2]
    coef = 6 * (4 * np.pi - omega) / omega
    q = {k: sum(b) / 2 for k, b in modular["modular_slopes"].items()}
    matched = transfers(delta, q)
    safe = live_survey.admission(d)["safe_cone_radius"]
    cone = {}
    for arm in ARMS:
        g = delta[arm][2]
        cone[arm] = dict(
            dI_per_r=float(-4 * (2 * np.pi + omega) * g / (3 * omega)),
            dA_per_r_squared=float(g / 3),
            modular=q[arm],
        )
    gain = cone["ordinary"]["dI_per_r"] / q["ordinary"]
    for arm in ARMS[1:]:
        cone[arm]["ordinary_prediction_per_r"] = gain * q[arm]
        cone[arm]["residual_per_r"] = cone[arm]["dI_per_r"] - gain * q[arm]
    return dict(
        divisor=divisor,
        continuous=continuous,
        reference=dict(zip(METRICS, map(float, vals), strict=True)),
        interior_deficits=deficits.tolist(),
        boundary_exterior_angles=boundary.tolist(),
        reference_mean_bulk_scalar=float(vals[0] / vals[1]),
        cone_radius_upper=safe,
        cone_mean_scalar_times_r_squared=float(coef),
        cone_volume_matched=cone,
        spherical_link_identity_error=float(abs(sum(deficits) + omega - 4 * np.pi)),
        source_slopes={
            k: dict(zip(METRICS, map(float, v), strict=True)) for k, v in delta.items()
        },
        maximum_real_derivative_discrepancy=fd_errors,
        matched=matched,
    )


def uniform():
    """Exact conditional shape bounds composed with native all-source clock bounds."""
    r = json.loads((HERE / "bounds.json").read_text())["original_local_six"]
    result = {}
    for mode in ("continuous", "spacing_1/16"):
        u = Fraction(3)
        ls = [Fraction(r["by_charge"][q][mode]["range_lower"]) for q in ("-1", "1")]
        edge_lower = (sum(ls) * 10**12 // 2) * Fraction(
            1, 10**12
        )  # Every center and source obeys each charge-specific bound.
        delta = u * u - edge_lower * edge_lower
        # Gram about the apex is within 3*delta in operator norm of the regular
        # face Gram, whose smallest eigenvalue is u²/2. Its face height squared
        # is 1/(ones^T G^-1 ones) >= lambda_min/3.
        eigen = u * u / 2 - 3 * delta
        if eigen <= 0:
            raise ValueError(
                "distance bounds do not admit uniform nondegenerate tetrahedra"
            )
        qlo = (2 * edge_lower * edge_lower - u * u) / 2
        qhi = (2 * u * u - edge_lower * edge_lower) / 2
        # For an arbitrary edge the other two projected face vectors have
        # dot product >= qlo-qhi²/l² and length product <=u².
        cos = (qlo - qhi * qhi / (edge_lower * edge_lower)) / (u * u)
        if not 0 < cos < 1:
            raise ValueError("distance bounds do not force acute interior dihedrals")
        # acos(x)<=pi/2-x for 0<=x<=1; each radial edge has three dihedrals.
        # Replace pi by its upper rational bound for the decreasing ratio.
        twopi = Fraction(44, 7)
        lower = 6 * (twopi + 12 * cos) / (twopi - 12 * cos)
        result[mode] = dict(
            uniform_edge_lower=str(edge_lower),
            uniform_edge_upper=str(u),
            minimum_tetrahedron_gram_eigenvalue_lower=str(eigen),
            safe_cone_radius_squared_lower=str(eigen / 3),
            dihedral_cosine_lower=str(cos),
            solid_angle_upper=str(twopi - 12 * cos),
            cone_mean_scalar_times_r_squared_lower=str(lower),
            lower_decimal=nextafter(float(lower), -inf),
        )
    return dict(
        stations=r["stations"][:5],
        source_scope="every retained density; pooled charge and orientation",
        constraints="Same glued three-dimensional star construction, using five of the previously bounded six enrolled stations.",
        bounds=result,
    )


UNIT = 1 << 52
DIVISOR = 16
N = 6 * DIVISOR
DEGREE = 22


def instrument_constants(h):
    L = Fraction(int(abs(h).sum(axis=1).max()), 40 * DIVISOR)
    assert L <= Fraction(3, 2)
    damp_exact = [
        sum(Fraction(-1, 2 * DIVISOR) ** j / factorial(j) for j in range(n))
        for n in (16, 17)
    ]
    damp = round(sum(damp_exact) * UNIT / 2)
    assert max(abs(x - Fraction(damp, UNIT)) for x in damp_exact) < Fraction(1, UNIT)
    # Each integer nearest rounding produces a column error <=sqrt(2D)/(2U).
    # Taylor's remainder for an exact unitary step has no exponential prefactor.
    term = total = 0
    for j in range(1, DEGREE + 1):
        term = L * term / j + Fraction(1 + isqrt(2 * h.shape[0]), 2 * UNIT)
        total += term
    # All exact previous steps contract. While total column error stays <1
    # the input norm is <=2. The computed pre-click norm is <=3. Damping
    # contributes <=3/U for its factor and <=2/U for 16 rounded components.
    col = N * (
        2 * L ** (DEGREE + 1) / factorial(DEGREE + 1) + total + Fraction(5, UNIT)
    )
    assert col < 1
    d = 3 * col
    branch = 2 * d + d * d
    p = 1 - Fraction(damp, UNIT) ** 2
    # True two-leg amplitude is a contraction. Replacing both pre-click
    # blocks costs <=branch; their Gram costs <=2branch+branch^2. Each
    # click supplies sqrt(p). |p_true^2-p^2|<=4/U.
    eps = p * p * (2 * branch + branch * branch) + Fraction(4, UNIT)
    return damp, eps * comb(N, 3) / (2 * DIVISOR), p * branch + Fraction(2, UNIT)


def pre_click(h, rows, b, damp):
    """Integral pre-click blocks for all senders, using complex transposition."""
    w = np.zeros((h.shape[0], 16), dtype=np.int64)
    w[rows[b], :8] = np.eye(8, dtype=np.int64) * UNIT
    norm = int(abs(h).sum(axis=1).max())
    hits = np.zeros((N + 1, len(rows), 16, 16), dtype=np.int64)
    for t in range(1, N + 1):
        term = w.copy()
        for j in range(1, DEGREE + 1):
            den = 40 * DIVISOR * j
            assert norm * int(abs(term).max()) + den < 2**63
            v = h @ term
            term[:, :8] = (v[:, 8:] + den // 2) // den
            term[:, 8:] = (-v[:, :8] + den // 2) // den
            assert int(abs(w).max()) + int(abs(term).max()) < 2**63
            w += term
        # U D U ... D U is complex symmetric. A launch at b thus returns
        # the transpose of every pre-click amplitude for an arbitrary a→b.
        re = w[rows, :8].swapaxes(1, 2)
        im = w[rows, 8:].swapaxes(1, 2)
        hits[t, :, :8, :8] = hits[t, :, 8:, 8:] = re
        hits[t, :, 8:, :8] = im
        hits[t, :, :8, 8:] = -im
        w[rows[b]] = (w[rows[b]].astype(object) * damp + UNIT // 2) // UNIT
    return hits


def rebate(f, g):
    """Sum (N-s-j) (G_j F_s)* (G_j F_s) without selecting a return."""
    f, g = f.astype(object), g.astype(object)
    total = np.zeros((16, 16), dtype=object)
    weight = np.zeros_like(total)
    prefix = np.zeros_like(total)
    # At k=N-s, weight=sum_{j=1}^k (k-j) G_j*G_j.
    for k in range(1, N):
        weight += prefix
        prefix += g[k].T @ g[k]
        x = f[N - k]
        total += x.T @ weight @ x
    return np.stack([total[:8, :8], total[8:, :8]])


@lru_cache(maxsize=1)
def exact_ranges():
    accumulated = np.zeros((10, 2, 8, 8), dtype=object)
    radius = 0
    for q in (-1, 1):
        h, allrows = flight.hardware(q)
        rows = allrows[:5]
        damp, err, _ = instrument_constants(h)
        terms = [pre_click(h, rows, b, damp) for b in range(5)]
        gain = (UNIT**2 - damp**2) ** 2
        for e, (a, b) in enumerate(EDGES):
            f, g = terms[b][:, a], terms[a][:, b]
            accumulated[e] -= gain * (rebate(f, g) + rebate(g, f))
        radius += err / 2
    den = 8 * DIVISOR * UNIT**8
    accumulated[:, 0] += 3 * den * np.eye(8, dtype=np.int64).astype(object)
    return accumulated, den, radius


def sources():
    """Exact pump centers and trace-norm radii for the six native readings."""
    import sympy as sp
    from .. import live_source, live_modular

    (formed, prod), ferr, _ = live_source.formation_center()
    _, _, col, grams = live_source.pump_center(True)
    e = 4 * col + 4 * col**2

    def f(x, weights=(2, 1)):
        return live_source.rational_read(grams, x, weights)

    ref = f(formed)
    out = dict(
        reference=(ref, ferr + e),
        ordinary=(f(formed, (1, 2)) - ref, 2 * (ferr + e)),
        decorrelated=(f(prod) - ref, 3 * ferr + 2 * e),
    )
    _, _, _, rw = live_source.exact_witness(grams, col, formed, ferr)
    axes = live_source.qubit_hidden.correlation_axes()
    direction = sum(
        (
            sp.Rational(a) * axes[j]
            for a, j in zip(rw["coefficients"], (8, 2, 7, 6), strict=True)
        ),
        sp.zeros(4),
    )
    de = sp.Rational(rw["l2_error"])
    out["hidden"] = (f(direction), de + (1 + de) * e)
    for tag, (rho, re) in live_modular.rotation_centers(formed, ferr).items():
        out[tag] = (f(rho), re + (1 + re) * e)
    return out


def sqrt_upper(x):
    scale = 2**85
    return Fraction(1 + isqrt(x.numerator * scale**2 // x.denominator), scale)


def edge_errors(edges, den, uniform_error):
    """Bound whole-history error before squaring its measured rebate amplitude.

    The orthogonal sum of conditional returns weighted by sqrt((T-s-j)/2)
    has Gram 3I-E. Each history amplitude has error <=p*(2d+d²)+2/U.
    Cauchy on that direct sum bounds its operator error by the square root
    of the sum of weights times this error. Multiplying by the ACTUAL
    computed rebate norm preserves sensitivity on the many rarely hit edges.
    """
    raw = max(instrument_constants(flight.hardware(q)[0])[2] for q in (-1, 1))
    delta = sqrt_upper(Fraction(comb(N, 3), 2 * DIVISOR)) * raw
    ident = 3 * den * np.eye(8, dtype=np.int64).astype(object)
    out = []
    for re, im in edges:
        norm = Fraction(int(np.max(np.sum(abs(ident - re) + abs(im), axis=1))), den)
        out.append(min(uniform_error, 2 * sqrt_upper(norm) * delta + delta**2))
    return out


def paired_inputs(edges, den, err, states):
    """Outward intervals about the known reference trace or zero tangent trace."""
    import sympy as sp
    from mpmath import iv
    from .. import live_bounds

    cols = {}
    norms = {}
    for tag, (x, e) in states.items():
        rr, ii = x.as_real_imag()
        cols[tag] = sp.Matrix([*rr, *ii])
        norms[tag] = sum(abs(v) for v in cols[tag]) + e
    out = {tag: [] for tag in states}
    ident = 3 * den * np.eye(8, dtype=int).astype(object)
    for (re, im), radius_e in zip(edges, err, strict=True):
        re = re - ident
        amp = sp.Rational(int(np.max(np.sum(abs(re) + abs(im), axis=1))), den)
        row = sp.Matrix([[*re.ravel().tolist(), *im.ravel().tolist()]])
        for tag, (_, e) in states.items():
            m = (3 if tag == "reference" else 0) + (row * cols[tag])[0] / den
            rad = amp * e + sp.Rational(radius_e) * norms[tag]
            ri = live_bounds.rational_interval
            # Source formation error can contain an exact algebraic norm.
            # Round its full expression upward, and any algebraic center
            # downward, on a rational mesh before passing it to MP intervals.
            mesh = 2**90
            middle = Fraction(int(sp.floor(m * mesh)), mesh)
            radius = Fraction(int(sp.ceiling(rad * mesh)) + 1, mesh)
            out[tag].append(ri(middle) + ri(radius) * iv.mpf([-1, 1]))
    return out


def radial_jets(lengths, variations):
    """Interior radial edge lengths and the sums of their three dihedral angles."""
    from mpmath import iv
    from .. import live_bounds

    Jet = live_bounds.Jet
    keys = live_bounds.KEYS
    zero = iv.mpf(0)
    hess = np.full((5, 5), zero, object)
    d = np.empty((5, 5), object)
    for e, (i, j) in enumerate(EDGES):
        d[i, j] = d[j, i] = Jet(
            lengths[e], np.array([variations[k][e] for k in keys], object), hess
        )

    def inverse(jet):
        return jet.unary(1 / jet.v, -1 / jet.v**2, 2 / jet.v**3)

    radial = [None] * 4
    # All three hinges of each face tetrahedron have their base at the apex.
    for face in BASES:
        for i in face:
            j, k = [v for v in face if v != i]

            def g(a, b):
                return (d[4, a] * d[4, a] + d[4, b] * d[4, b] - d[a, b] * d[a, b]) * (
                    iv.mpf(1) / 2
                )

            x, y = g(i, j), g(i, k)
            r2 = inverse(d[4, i] * d[4, i])
            nx = d[4, j] * d[4, j] - x * x * r2
            ny = d[4, k] * d[4, k] - y * y * r2
            c = (g(j, k) - x * y * r2) * inverse((nx * ny).sqrt())
            sine = iv.sqrt(1 - c.v**2)
            angle = c.unary(iv.atan2(sine, c.v), -1 / sine, -c.v / sine**3)
            radial[i] = angle if radial[i] is None else radial[i] + angle
    return [(d[4, i], radial[i]) for i in range(4)]


def bulk_jet(lengths, variations):
    from mpmath import iv

    return sum(length * (-angle + 2 * iv.pi) * 2
               for length, angle in radial_jets(lengths, variations))


def total_slopes(lengths, variations):
    """Schlaefli identity d(I+boundary)=2 sum_e bending_e d length_e."""
    from mpmath import iv

    values = dict(zip(EDGES, lengths, strict=True))
    angles = {e: (2 if 4 in e else 1) * iv.pi for e in EDGES}
    for face in BASES:
        tet = (*face, 4)
        for a, b in combinations(tet, 2):
            other = [v for v in tet if v not in (a, b)]

            def d(i, j):
                return values[tuple(sorted((i, j)))]

            def g(i, j):
                return (d(a, i) ** 2 + d(a, j) ** 2 - d(i, j) ** 2) / 2

            c, k = other
            base = d(a, b) ** 2
            x, y = g(b, c), g(b, k)
            cosine = (g(c, k) - x * y / base) / iv.sqrt(
                (d(a, c) ** 2 - x**2 / base) * (d(a, k) ** 2 - y**2 / base)
            )
            angles[a, b] -= iv.atan2(iv.sqrt(1 - cosine**2), cosine)
    return {
        k: 2 * sum(angles[e] * v for e, v in zip(EDGES, dx, strict=True))
        for k, dx in variations.items()
    }


def determinant(a, b, c):
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def certificate(samples):
    from mpmath import iv
    from .. import live_bounds

    centers, den, err = exact_ranges()
    saved = iv.dps
    try:
        iv.dps = 35
        radii = edge_errors(centers, den, err)
        inputs = paired_inputs(centers, den, radii, sources())
        lengths = inputs.pop("reference")
        R = bulk_jet(lengths, inputs)
        total = total_slopes(lengths, inputs)
        A, V, omega = live_bounds.geometry_jets(lengths, inputs)
        ep = live_bounds.endpoints
        idx = {k: i for i, k in enumerate(live_bounds.KEYS)}
        q = {k: iv.mpf(b) for k, b in samples["modular_slopes"].items()}
        methods = {}
        for side in ROTATIONS:
            r = idx[side]
            if V.g[r].a <= 0 <= V.g[r].b:
                methods[side] = dict(volume_slope_unresolved=ep(V.g[r]))
                continue
            matched = {}
            tot = {}
            probe = "left" if side == "right" else "right"
            for tag in (*ARMS, probe):
                i = idx[tag]
                alpha = -V.g[i] / V.g[r]
                matched[tag] = (R.g[i] + alpha * R.g[r], q[tag] + alpha * q[side])
                tot[tag] = total[tag] + alpha * total[side]
            ko = matched["ordinary"][1]
            if ko.a <= 0 <= ko.b:
                methods[side] = dict(ordinary_modular_unresolved=ep(ko))
                continue
            gain = matched["ordinary"][0] / ko
            out = {}
            for tag, (meas, kval) in matched.items():
                out[tag] = dict(
                    bulk_scalar=ep(meas),
                    modular=ep(kval),
                    bulk_plus_boundary=ep(tot[tag]),
                )
                if tag != "ordinary":
                    predicted = gain * kval
                    out[tag].update(
                        predicted=ep(predicted),
                        residual=ep(meas - predicted),
                        total_ordinary_residual=ep(
                            tot[tag] - tot["ordinary"] / ko * kval
                        ),
                    )
            rank = determinant(
                *[(matched[tag][0], tot[tag], matched[tag][1]) for tag in ARMS]
            )
            methods[side] = dict(
                ordinary_calibrated_gain=ep(gain),
                sources=out,
                determinant_bulk_total_modular=ep(rank),
            )
            io = matched["ordinary"][0]
            ic = matched["decorrelated"][0]
            to = tot["ordinary"]
            tc = tot["decorrelated"]
            minor = io * tc - ic * to
            if minor.a > 0 or minor.b < 0:
                witnesses = {}
                for third in ("hidden", probe):
                    it = matched[third][0]
                    tt = tot[third]
                    det = determinant(
                        *[
                            (matched[t][0], tot[t], matched[t][1])
                            for t in ("ordinary", "decorrelated", third)
                        ]
                    )
                    witnesses[third] = dict(
                        ordinary_coefficient=ep((ic * tt - tc * it) / minor),
                        decorrelation_coefficient=ep((to * it - io * tt) / minor),
                        modular_derivative_at_fixed_geometry=ep(det / minor),
                    )
                methods[side]["implicit_geometry_constraints"] = dict(
                    volume_and_two_curvatures_control_minor=ep(V.g[r] * minor),
                    volume_matched_curvature_minor=ep(minor),
                    witnesses=witnesses,
                )

            if side == "right":
                native = ("ordinary", "decorrelated", "left")
                area = {t: A.g[idx[t]] - V.g[idx[t]] / V.g[r] * A.g[r] for t in native}
                (region,) = (
                    v
                    for v in json.loads((HERE / "regional.json").read_text())["regions"]
                    if v["radius"] == 1
                )
                rk = {
                    t: iv.mpf(region["regional_entropy_slopes"][t])
                    for t in (*native, side)
                }
                regional = {t: rk[t] - V.g[idx[t]] / V.g[r] * rk[side] for t in native}
                io, ic, il = [matched[t][0] for t in native]
                ao, ac, al = [area[t] for t in native]
                minor = io * ac - ic * ao
                if minor.a > 0 or minor.b < 0:
                    kd = (
                        determinant(
                            *[(matched[t][0], area[t], matched[t][1]) for t in native]
                        )
                        / minor
                    )
                    kr = (
                        determinant(
                            *[(matched[t][0], area[t], regional[t]) for t in native]
                        )
                        / minor
                    )
                    methods[side]["equal_volume_outer_area_bulk"] = dict(
                        control_minor=ep(V.g[r] * minor),
                        ordinary_per_left_radian=ep((ic * al - ac * il) / minor),
                        partial_repairing_per_left_radian=ep(
                            (ao * il - io * al) / minor
                        ),
                        retained_entropy_per_left_radian=ep(kd),
                        twelve_port_entropy_per_left_radian=ep(kr),
                        regional_ports=region["ports"],
                        regional_scope="Both original heads, strain and enclosed spent-pump modes on these twelve ports. Their identification with a proper ball of the five-station star is open.",
                    )

        return dict(
            divisor=DIVISOR,
            pooled_range_effect_operator_radius=nextafter(float(max(radii)), inf),
            edge_operator_radii=[nextafter(float(r), inf) for r in radii],
            reference_bulk_scalar=ep(R.v),
            filled_volume=ep(V.v),
            solid_angle=ep(omega.v),
            primary_cone_mean_times_r_squared=ep(6 * (4 * iv.pi - omega.v) / omega.v),
            raw_volume_slopes={k: ep(V.g[i]) for k, i in idx.items()},
            matched=methods,
        )
    finally:
        iv.dps = saved


def finite_curvature(base, change):
    """Curvature of the perturbed lengths with enough digits for tiny secants."""
    import mpmath as mp

    length = mp.matrix(5)
    for a, b in EDGES:
        length[a, b] = length[b, a] = mp.mpf(float(base[a, b])) + mp.mpf(
            float(change[a, b])
        )

    def gram(i, j):
        return (length[4, i] ** 2 + length[4, j] ** 2 - length[i, j] ** 2) / 2

    defect = [2 * mp.pi for _ in range(4)]
    for face in BASES:
        for i in face:
            j, k = [t for t in face if t != i]
            g = gram(i, i)
            x, y = gram(i, j), gram(i, k)
            c = (gram(j, k) - x * y / g) / mp.sqrt(
                (gram(j, j) - x * x / g) * (gram(k, k) - y * y / g)
            )
            defect[i] -= mp.acos(c)
    return 2 * sum(length[4, i] * a for i, a in enumerate(defect))


def finite_matches(first_panel):
    """Rotate onto the measured volume level set using actual positive preparations."""
    from scipy.optimize import brentq
    import mpmath as mp
    from .. import live_matching, live_source

    source = flight.live_radar.source
    board = source.formation.federation(1)
    rho, product = source.form(board, source.formation.region(board), -1)
    g = 2 * live_survey.rotations(rho)["right"]
    states = live_source.postpump(True)
    ref = states["reference"]
    channel = live_source.head_channel(True)
    center, den, _ = exact_ranges()
    es = (np.array(center[:, 0], float) + 1j * np.array(center[:, 1], float)) / float(
        den
    )
    base = live_survey.square(live_survey.read(es, ref), 5)
    # Exact changes are traceless; remove the common 3I timeout charge first.
    contrasts = es - 3 * np.eye(8)[None]

    def distance_delta(dx):
        return live_survey.square(live_survey.read(contrasts, dx), 5)

    vals, axes = np.linalg.eigh(ref)
    log = -(axes * np.log(vals)) @ axes.conj().T
    panels = []
    with mp.workdps(42):
        zero = np.zeros((5, 5))
        v0 = live_matching.area_volume(base, zero)[1]
        i0 = finite_curvature(base, zero)
        for t in (1e-6, 1e-7, 1e-8):
            arms = {}
            for tag in ("ordinary", "decorrelated"):
                weights = (1 - t) * np.array(source.RAILS["reference"]) + t * np.array(
                    source.RAILS[tag]
                )
                head = rho if tag == "ordinary" else (1 - t) * rho + t * product
                # Unit conjugation about G/2 changes the incoming head by
                # (cos(theta)-1) C + sin(theta) S. The common rail pump is linear.
                c = live_source.apply_channel(
                    channel, (head - g @ head @ g) / 2, weights
                )
                s = live_source.apply_channel(
                    channel, -0.5j * (g @ head - head @ g), weights
                )
                d0 = t * (states[tag] - ref)

                def change(theta):
                    return d0 - 2 * np.sin(theta / 2) ** 2 * c + np.sin(theta) * s

                def volume(theta):
                    return float(
                        live_matching.area_volume(base, distance_delta(change(theta)))[
                            1
                        ]
                        - v0
                    )

                initial = (
                    first_panel["matched"]["right"]["sources"][tag]["rotation"] * t
                )
                bracket = sorted((0, 2 * initial))
                theta = brentq(volume, *bracket, xtol=1e-19)
                dx = change(theta)
                ivalue = finite_curvature(base, distance_delta(dx))
                spectrum = np.linalg.eigvalsh(ref + dx)
                arms[tag] = dict(
                    rotation=theta,
                    volume_error=volume(theta),
                    bulk_scalar_secant=float((ivalue - i0) / t),
                    modular_secant=float(np.trace(log @ dx).real / t),
                    entropy_secant=float(
                        (vals @ np.log(vals) - spectrum @ np.log(spectrum)) / t
                    ),
                )
            a, b = (arms[tag] for tag in ("ordinary", "decorrelated"))
            gain = a["bulk_scalar_secant"] / a["modular_secant"]
            panels.append(
                dict(
                    admixture=t,
                    arms=arms,
                    ordinary_gain=gain,
                    held_residual=b["bulk_scalar_secant"] - gain * b["modular_secant"],
                )
            )
    return panels


def three_matches(first_panel, use_area=False):
    """Finite checks of the exact local implicit family on native one-sided preparations."""
    import mpmath as mp
    from scipy.optimize import root
    from .. import live_matching, live_source

    source = flight.live_radar.source
    board = source.formation.federation(1)
    formed, product = source.form(board, source.formation.region(board), -1)
    generators = live_survey.rotations(formed)
    ref = live_source.postpump(True)["reference"]
    channel = live_source.head_channel(True)
    centers, den, _ = exact_ranges()
    es = (np.array(centers[:, 0], float) + 1j * np.array(centers[:, 1], float)) / float(
        den
    )
    base = live_survey.square(live_survey.read(es, ref), 5)
    dif = es - 3 * np.eye(8)[None]
    vals, axes = np.linalg.eigh(ref)
    log = -(axes * np.log(vals)) @ axes.conj().T
    zeros = np.zeros((5, 5))
    ref_weights = np.array(source.RAILS["reference"])
    rail_change = np.array(source.RAILS["ordinary"]) - ref_weights
    # Schlaefli gives the exact derivative of I+B along a length segment.
    # Four-point Gauss integration evaluates its very small finite change
    # directly; it never subtracts two large nearly equal Regge actions.
    pts, weights = np.polynomial.legendre.leggauss(4)

    def change(o, c, theta, turn):
        head = formed + c * (product - formed)
        delta = c * (product - formed)
        for side, angle in (("left", turn), ("right", theta)):
            g = 2 * generators[side]
            x = -(np.sin(angle / 2) ** 2) * (head - g @ head @ g) - 0.5j * np.sin(
                angle
            ) * (g @ head - head @ g)
            head = head + x
            delta += x
        w = ref_weights + o * rail_change
        return (
            live_source.apply_channel(channel, delta, w)
            + o * live_source.apply_channel(channel, formed, rail_change)
        ), w

    charts = {}
    for tag, record in first_panel["source_slopes"].items():
        charts[tag] = [
            record["filled_volume"],
            record["bulk_scalar"],
            record["outer_area"]
            if use_area
            else record["bulk_scalar"] + record["boundary_bending"],
        ]
    J = np.array([charts[k] for k in ("ordinary", "decorrelated", "right")]).T
    initial = -np.linalg.solve(J, charts["left"])
    with mp.workdps(42):
        a0, v0 = live_matching.area_volume(base, zeros)
        i0 = finite_curvature(base, zeros)

        def residual(dx):
            dl = live_survey.square(live_survey.read(dif, dx), 5)
            area, volume = live_matching.area_volume(base, dl)
            if use_area:
                return np.array(
                    [
                        float(volume - v0),
                        float(finite_curvature(base, dl) - i0),
                        float(area - a0),
                    ]
                )
            total = 0.0
            for t, w in zip(pts, weights, strict=True):
                _, rad, rim = bend(base + dl * ((t + 1) / 2))
                total += w * (
                    dl[4, :4] @ rad
                    + np.array([dl[i, j] for i, j in combinations(range(4), 2)]) @ rim
                )
            return np.array(
                [
                    float(volume - v0),
                    float(finite_curvature(base, dl) - i0),
                    total,
                ]
            )

        result = []
        for turn in (-1e-4, -1e-5, -1e-6):

            def loss(x):
                dx, _ = change(*(turn * x), turn)
                return np.linalg.solve(J, residual(dx) / turn)

            fit = root(loss, initial, tol=1e-10)
            o, c, theta = fit.x * turn
            dx, w = change(o, c, theta, turn)
            geometry = residual(dx)
            if max(abs(geometry)) > 2e-17 or not (
                0 < c < 1 and 0 < w.min() <= w.max() < 1
            ):
                raise ArithmeticError(
                    "three finite geometric readings were not matched on a positive preparation"
                )
            spectrum = np.linalg.eigvalsh(ref + dx)
            result.append(
                dict(
                    left_rotation=turn,
                    right_rotation=float(theta),
                    rail_lottery=w.tolist(),
                    ordinary_coordinate=float(o),
                    partial_repairing=float(c),
                    geometry_error=dict(
                        zip(
                            (
                                "volume",
                                "bulk_scalar",
                                "outer_area" if use_area else "bulk_plus_boundary",
                            ),
                            map(float, geometry),
                            strict=True,
                        )
                    ),
                    minimum_retained_eigenvalue=float(spectrum.min()),
                    modular_change=float(np.trace(log @ dx).real),
                    entropy_change=float(
                        vals @ np.log(vals) - spectrum @ np.log(spectrum)
                    ),
                )
            )
    return result


def run():
    record = json.loads((HERE / "flight.json").read_text())
    modular = json.loads((HERE.parent / "live_modular.json").read_text())
    panels = [
        panel(record, modular, d, c) for c in (False, True) for d in flight.DIVISORS
    ]
    return dict(
        finite_appointment_certificate=certificate(modular),
        finite_positive_right_matches=finite_matches(panels[0]),
        finite_equal_volume_bulk_and_boundary=three_matches(panels[0]),
        finite_equal_volume_outer_area_and_bulk=three_matches(panels[0], use_area=True),
        primary_stations=list(flight.STATIONS[:5]),
        apex=flight.STATIONS[4],
        hinges="four interior radial edges",
        source="Coupled four-tick vacancy-pumped state with live strain and both original heads.",
        modular="Independent retained eight-state reference negative log. Ordinary and decorrelated are source slopes about the same reference; hidden is nominal XZ+ZX and its modular interval encloses the independently refined tangent.",
        normalization="Twice sum of length times interior dihedral deficit. Boundary bending is reported separately.",
        uncertainty="The 1/16 primary transfer uses enclosed native effects and source states. Continuous pooled primary panels and finite secants are numerical checks. The held five-station cone bounds are exact and uniform.",
        panels=panels,
        uniform_held_star=uniform(),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    report = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if parser.parse_args().write:
        (HERE / "hinge.json").write_text(report)
    else:
        print(report, end="")
