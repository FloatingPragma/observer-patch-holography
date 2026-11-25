"""Signed optical changes on the enclosed equilibrium, with exact native aperture effects."""

import json
import sympy as sp

from ...full_algebra import certificate as cert
from ...full_algebra.native_precision import exact_effects
from ...null_optics.certificate import bounded_sources
from . import measure


def run():
    cert.iv.dps = 30
    ri, band, iv = cert.ri, cert.ri_band, cert.iv
    src = bounded_sources(sp.eye(4) / 4)
    x, eps = src["state"]
    unpair, unpair_err = src["unpair"]
    ordinary = [sum(x, sp.zeros(4)) / 2 - b for b in x]
    product = [a + d for a, d in zip(x, unpair)]

    def inner(a, b):
        return sum(cert.inner(p, v) for p, v in zip(a, b))

    def row_bound(e):
        return max(cert.prior.one_norm(b) for b in e)

    # The independently certified angle bounds are widened to exact integers.
    named = json.loads(
        (measure.HERE.parents[1] / "null_optics/certificate.json").read_text()
    )["matched_pairing"]
    angle = {}
    for key, field, low, high in (
        ("ordinary", "ordinary_angle_interval", -997, -983),
        ("unpairing", "angle_interval", -297, -291),
    ):
        interval = named[field]
        assert low < interval[0] < interval[1] < high
        angle[key] = iv.mpf([low, high])
    pauli = [
        sp.Matrix([[0, 1], [1, 0]]),
        sp.Matrix([[0, -sp.I], [sp.I, 0]]),
        sp.diag(1, -1),
    ]
    gens = [sp.kronecker_product(sp.eye(2), a) for a in pauli]
    bloch = [band(inner([a, a], x), eps) for a in gens]
    length = iv.sqrt(sum(a**2 for a in bloch))
    axis = [a / length for a in bloch]
    unit_rots = [
        [(-sp.I * (a @ v - v @ a) / 2).applyfunc(sp.expand) for v in x] for a in gens
    ]

    def rotated(e):
        return sum(
            component * band(inner(e, d), row_bound(e) * eps)
            for component, d in zip(axis, unit_rots)
        )

    def raw(e, d, error):
        return band(inner(e, d), row_bound(e) * error)

    nums, unit, error = exact_effects()
    ge = cert.rational_effects(nums, unit)
    ranges = [band(inner(e, x), 3 * eps + sp.Rational(error)) for e in ge]
    _, v = cert.enc.geometry_jet(ranges)
    dirs = [[cert.prior.sym_matrix(b) for b in d] for d in measure.directions()]
    vs = [
        sum(g * band(inner(e, d), 2 * sp.Rational(error)) for g, e in zip(v.g, ge))
        for d in dirs
    ]
    vn = iv.sqrt(sum(a**2 for a in vs))
    phase = [-vs[1] / vn, vs[0] / vn]
    assert vs[0].b < 0 and vs[1].a > 0

    # The volume adjuster preserves both modular terms by unitary invariance.
    # Bound -log chi and -log pi against the exact central ordinary perturbation.
    def entropy_slope(blocks, state_error):
        center = sp.S.Zero
        log_error = sp.S.Zero
        for b, d in zip(blocks, ordinary):
            val, err = cert.prior.entropy_coordinates(
                b, state_error, sp.Rational(".019"), [d]
            )
            center += val[0]
            log_error = max(log_error, err)
        # The ordinary direction has trace norm <=2 and error <=2 eps.
        # All true block eigenvalues are in [.019,1], so ||log|| < 4.
        return band(center, 2 * log_error + 8 * eps)

    ordinary_full = entropy_slope(x, sp.Rational("2e-7"))
    ordinary_collar = ordinary_full - entropy_slope(product, 2 * eps)
    assert ordinary_collar.a > 0
    board = measure.front.optics.pairing.formation.federation(1)
    route = measure.front.optics.pairing.formation.region(board)
    ops = measure.operators(board, route)
    output = []
    for raw_blocks in ops:
        e = [cert.prior.sym_matrix(h) / (4 * measure.NORM) for h in raw_blocks]
        f = raw(e, x, eps)
        dl_scale = -((2304 * f) ** (-ri(1) / 6)) / (6 * f)
        calc = {
            "ordinary": raw(e, ordinary, 2 * eps) + angle["ordinary"] * rotated(e),
            "unpairing": raw(e, unpair, unpair_err) + angle["unpairing"] * rotated(e),
            "phase": sum(w * ri(inner(e, d)) for w, d in zip(phase, dirs)),
        }
        dl = {tag: dl_scale * z for tag, z in calc.items()}
        assert dl["ordinary"].a > 0 and dl["unpairing"].a > 0 and dl["phase"].b < 0
        # Prior inward/volume-matched entropy and hidden-phase bounds in nats.
        charges = {
            "unpairing": iv.mpf([".231079", ".231133"]),
            "phase": iv.mpf([".0003224", ".0003260"]),
        }
        fractions = {
            kind: {
                tag: str(dl[tag] * ord_charge / (dl["ordinary"] * charges[tag]))
                for tag in charges
            }
            for kind, ord_charge in (
                ("full", ordinary_full),
                ("collar", ordinary_collar),
            )
        }
        output.append(
            dict(slope={k: str(z) for k, z in dl.items()}, transfer_ratio=fractions)
        )
    return dict(
        phase_volume_basis=[str(z) for z in vs],
        phase_weights=[str(z) for z in phase],
        ordinary_full_modular=str(ordinary_full),
        ordinary_collar_modular=str(ordinary_collar),
        apertures=output,
    )


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
