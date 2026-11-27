"""Interval submersion of two early tetrahedra on the exact late-range fiber."""

import argparse
import json

import numpy as np
from scipy.linalg import qr
import sympy as sp

from .... import certificate as old
from ..probe import density_center
from . import kernels

iv = old.iv
ri = old.ri
EPS = sp.Rational(9, 10**16)


def abs_bound(a):
    return max(sum(abs(x) for x in a.row(i)) for i in range(a.rows))


def upward(x):
    return float(np.nextafter(float(x.b), np.inf))


def geometry_inputs():
    x, *_ = density_center()
    center = [old.prior.sym_matrix(b) for b in x]
    rho_error = sp.Rational("0.00000020000001")
    base = [old.prior.sym_matrix(v) for v in old.iso.basis(True)]
    with np.load(kernels.HERE / "deficits.npz") as raw:
        data = raw["numerators"]
    maps, jets = [], []
    for cap, source in zip(kernels.CAPS, data, strict=True):
        ranges, rows = [], []
        for blocks in source:
            herm = []
            for re, im in blocks:
                mat = (sp.Matrix(re) + sp.I * sp.Matrix(im)) / kernels.GRID
                herm.append((mat + mat.H) / 2)
            row = [old.inner(h, b) for h in herm for b in base]
            rows.append(row)
            # The full state has unit trace. Only the deficit can change.
            shift = sum(old.inner(a, b) for a, b in zip(herm, center, strict=True))
            width = ri(rho_error) * iv.sqrt(
                ri(sum(old.inner(a, a) for a in herm))
            ) + ri(EPS)
            ranges.append(ri(sp.Rational(cap, 2) + shift) + iv.mpf([-1, 1]) * width)
        maps.append(sp.Matrix(rows))
        jets.append(old.enc.geometry_jet(ranges))
    return maps, jets


def complement(r):
    # Solve first on the rational measured center, then correct that solution
    # on the same native physical rows. The free coordinates stay fixed.
    piv = qr(np.array(r, dtype=float), mode="economic", pivoting=True)[2][:12].tolist()
    free = [j for j in range(18) if j not in piv]
    inv = r[:, piv].inv(method="DM")
    z = sp.zeros(18, 6)
    for j, k in enumerate(free):
        z[k, j] = 1
    projection = -inv * r[:, free]
    for i, k in enumerate(piv):
        z[k, :] = projection[i, :]
    assert r * z == sp.zeros(12, 6)
    # Each Pauli-product column has trace norm 2, hence map-entry error 2 EPS.
    defect = abs_bound(inv) * 12 * 2 * EPS
    assert defect < sp.Rational(1, 1000)
    inv_bound = abs_bound(inv) / (1 - defect)
    return piv, z, inv, inv_bound, defect


def evaluate():
    iv.dps = 45
    maps, jets = geometry_inputs()
    piv, z, inv, inv_bound, defect = complement(maps[-1])
    derivatives = []
    for m, (area, vol) in zip(maps[:-1], jets[:-1], strict=True):
        reduced = m * z
        entries = np.empty((12, 6), dtype=object)
        carries = [
            (sum(abs(v) for v in m[i, piv] * inv) + defect) / (1 - defect)
            for i in range(12)
        ]
        for j in range(6):
            zsize = sum(abs(x) for x in z.col(j))
            for i in range(12):
                # Enclose M_P R_P^{-1} jointly before bounding its error.
                # A row times the inverse contains much tighter cancellations
                # than the norms of the row and inverse considered separately.
                err = 2 * EPS * zsize * (1 + carries[i])
                entries[i, j] = old.ri_band(reduced[i, j], err)
        derivatives.extend([area.g @ entries, vol.g @ entries])
    d = np.array(derivatives, dtype=object)
    # A numerical proposal for four columns is harmless: every residual
    # with this particular dyadic inverse is checked using outward intervals.
    mid = np.array([[float(v.mid) for v in row] for row in d])
    scale = np.max(abs(mid), axis=1)
    keep = qr(mid / scale[:, None], mode="economic", pivoting=True)[2][:4].tolist()
    # Left preconditioning is better scaled after normalizing the four rows.
    row_scale = np.array([2 ** round(np.log2(1 / s)) for s in scale])
    scaled = d * row_scale[:, None]
    a = mid[:, keep] * row_scale[:, None]
    proposed = np.linalg.inv(a)
    u = np.array(
        [[ri(sp.Rational(round(x * 2**30), 2**30)) for x in row] for row in proposed],
        dtype=object,
    )
    product = scaled[:, keep] @ u
    error = max(
        upward(sum(abs(product[i, j] - int(i == j)) for j in range(4)))
        for i in range(4)
    )
    assert error < 1, error
    return dict(
        late_range_rank=12,
        fixed_marginal_coordinates=18,
        fiber_dimension=6,
        early_geometry_rank=4,
        early_rows=["area 2", "volume 2", "area 4", "volume 4"],
        right_inverse_residual_infinity_upper=error,
        range_inverse_infinity_upper=float(np.nextafter(float(inv_bound), np.inf)),
        late_pivot_indices=piv,
        early_free_pivots=keep,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    report = evaluate()
    path = kernels.HERE / "certificate.json"
    if args.write:
        with path.open("x") as out:
            out.write(json.dumps(report, indent=2) + "\n")
    else:
        saved = json.loads(path.read_text())
        assert saved["early_geometry_rank"] == report["early_geometry_rank"]
        assert saved["right_inverse_residual_infinity_upper"] < 1
        print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
