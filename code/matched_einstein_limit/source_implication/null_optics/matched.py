"""Right-marginal rotation to keep the cap-six tetrahedral volume fixed."""

import numpy as np

from .. import signed_transfers


WEIGHTS = (1000, -883, -117, 50, -117, 67, 62, -112, 50, -112, -888, 1000)


def center(states):
    x = states["state"]
    right = signed_transfers.right_coin(x)
    pauli = signed_transfers.phase.AXES[1:]
    bloch = np.array([np.trace(right @ p).real for p in pauli])
    h = np.kron(
        np.eye(2),
        sum(v * p for v, p in zip(bloch, pauli)) / (2 * np.linalg.norm(bloch)),
    )
    rotation = -1j * (h @ x - x @ h)
    effects = signed_transfers.station.radar_menus()[6]
    station = signed_transfers.station
    g = station.range_gradient(station.pair(effects, x))[1]

    def slope(y):
        return g @ station.pair(effects, y)

    result = {}
    for name in ("ordinary", "unpair"):
        angle = -slope(states[name]) / slope(rotation)
        fixed = states[name] + angle * rotation
        assert abs(slope(fixed)) < 1e-12
        result["matched_" + name] = fixed
    return result


def certified_geometry(kernels, sources, clock_effect, clock_error):
    import sympy as sp
    from ..full_algebra import certificate as cert
    from ..isometric_collar import certify as PRIOR
    from .measure import PAIRS

    iv = cert.iv
    iv.dps = 30
    ri = cert.ri
    band = cert.ri_band
    eps = sp.Rational(PRIOR.error_effect())
    nums = np.array(
        [
            sum(
                np.array(PRIOR.direction_effects(kernels[b, q], 12), dtype=object)
                for q in (1, -1)
            )
            for b in (False, True)
        ],
        dtype=object,
    ).swapaxes(0, 1)
    ge = cert.rational_effects(nums, 8 * PRIOR.FRAC)
    mid, rho_error = sources["state"]

    def inner(e, d):
        return sum(cert.inner(a, b) for a, b in zip(e, d))

    ranges = [band(inner(e, mid), eps + 3 * rho_error) for e in ge]
    _, vol = cert.enc.geometry_jet(ranges)
    # The denominator comes from the existing rigorous full-algebra certificate
    # on exactly this center and radius; the code there encloses it more tightly.
    # Recheck that bound here via commutators on the same rational two-way effects.
    d = sp.Rational("2e-7")
    pauli = [
        sp.Matrix([[0, 1], [1, 0]]),
        sp.Matrix([[0, -sp.I], [sp.I, 0]]),
        sp.diag(1, -1),
    ]
    gens = [sp.kronecker_product(sp.eye(2), p) for p in pauli]
    b = [band(sum(cert.inner(x, g) for x in mid), rho_error) for g in gens]
    norm = iv.sqrt(sum(x**2 for x in b))
    assert norm.a > ri(".026")
    axis = [a / norm for a in b]
    generators = [
        [(-sp.I * (g * x - x * g) / 2).applyfunc(sp.expand) for x in mid] for g in gens
    ]

    def rotated(e, unc):
        mu = sum(a.trace() for a in e) / 8
        square = sum(cert.inner(a - mu * sp.eye(4), a - mu * sp.eye(4)) for a in e)
        err = ri(d) * iv.sqrt(ri(square)) + ri(unc) * 2
        # Dimensionless generators have operator norm one, so commutation with
        # g/2 costs at most one in Hilbert-Schmidt norm. The effects are centered.
        return sum(
            component * (ri(inner(e, gen)) + iv.mpf([-1, 1]) * err)
            for component, gen in zip(axis, generators)
        )

    vr = sum(g * rotated(e, eps) for g, e in zip(vol.g, ge))
    assert vr.b < ri("-.000025")

    def match(direction, delta_err, weights, all_positive):
        sampled = [band(inner(e, direction), 3 * delta_err + 2 * eps) for e in ge]
        angle = -sum(g * u for g, u in zip(vol.g, sampled)) / vr
        rows = []
        for q in (1, -1):
            for sender in range(4):
                picks = [
                    k for k, (i, j) in enumerate(PAIRS) if i == sender and weights[k]
                ]
                if not picks:
                    continue
                mass = sum(abs(weights[k]) for k in picks)
                pair = []
                for parity in (False, True):
                    a = sum(
                        (
                            weights[k] * clock_effect(kernels[parity, q], *PAIRS[k], 12)
                            for k in picks
                        ),
                        sp.zeros(4),
                    )
                    pair.append((a + a.H) / 2)
                unc = mass * clock_error
                raw = band(inner(pair, direction), 6 * mass * delta_err + 2 * unc)
                result = raw + angle * rotated(pair, unc)
                if all_positive:
                    result = q * result
                    assert result.a > 0
                else:
                    assert result.a > 0 or result.b < 0
                rows.append(
                    dict(
                        charge=q,
                        sender_index=sender,
                        lower=float(result.a),
                        upper=float(result.b),
                    )
                )
        return angle, rows

    delta, delta_err = sources["unpair"]
    angle, rows = match(delta, delta_err, [sp.Rational(k, 1000) for k in WEIGHTS], True)
    checker = [1, -1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1]
    _, pair_splits = match(delta, delta_err, checker, False)
    total = sum(mid, sp.zeros(4))
    oa, ord_splits = match([total / 2 - a for a in mid], 2 * rho_error, checker, False)
    return dict(
        weights=list(WEIGHTS),
        paired_four_edge_splits=pair_splits,
        ordinary_four_edge_splits=ord_splits,
        ordinary_angle_interval=[float(oa.a), float(oa.b)],
        scale=1000,
        angle_interval=[float(angle.a), float(angle.b)],
        volume_rotation_interval=[float(vr.a), float(vr.b)],
        minimum_signed_lower=min(row["lower"] for row in rows),
        rows=rows,
    )
