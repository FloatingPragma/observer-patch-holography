"""Transfer the counted connected charge to the area of the calibrated radar survey."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
import sympy as sp

from matched_einstein_limit.source_implication.full_algebra import certificate as cert
from matched_einstein_limit.source_implication.full_algebra.reuse.retained_equilibrium.nested_ranges import (
    kernels,
)
from . import measure


def run():
    prior_precision = cert.iv.dps
    try:
        cert.iv.dps = 34
        return enclose_transfer()
    finally:
        cert.iv.dps = prior_precision


def enclose_transfer():
    # All input matrix components below are exact binary rationals. The
    # cap-six propagation error plus 4x4 entry rounding is < 1e-15.
    iv = cert.iv
    bank = kernels.HERE / "deficits.npz"
    with np.load(bank) as data:
        e = cert.rational_effects(data["numerators"][-1].astype(object), kernels.GRID)
    root = Path(cert.__file__).parent
    record = json.loads((root / "results.json").read_bytes())
    centers = [
        cert.prior.sym_matrix(block)
        for block in np.array(record["x_real"]) + 1j * np.array(record["x_imag"])
    ]
    centers = [(a + a.H) / 2 for a in centers]
    norm = sum(a.trace() for a in centers)
    centers = [a / norm for a in centers]
    recipe = json.loads((measure.HERE / "results.json").read_bytes())["source"]
    step = sp.Rational(recipe["step"])
    basis = measure.joint_basis()
    D = [
        sum(
            (
                sp.Rational(t) * cert.prior.sym_matrix(op[b])
                for t, op in zip(recipe["basis_coefficients"], basis, strict=True)
            ),
            sp.zeros(4),
        )
        * 2
        * step
        for b in (0, 1)
    ]
    hs = iv.sqrt(cert.ri(sum(cert.inner(a, a) for a in D)))
    assert hs.b < 0.016
    assert all(a.trace() == 0 for a in D)
    for a in D:
        for i, j in np.ndindex(2, 2):
            assert sum(a[2 * b + i, 2 * b + j] for b in (0, 1)) == 0
            assert sum(a[2 * i + b, 2 * j + b] for b in (0, 1)) == 0

    def ranges(blocks, error):
        return np.array(
            [
                cert.ri_band(
                    sum(cert.inner(a, b) for a, b in zip(pair, blocks, strict=True)),
                    sp.Rational(error),
                )
                for pair in e
            ],
            dtype=object,
        )

    r = ranges(centers, "1e-15")
    dr = ranges(D, "1e-15")  # ||D||_1 <= sqrt(8) ||D||_HS < 1
    area, vol = cert.enc.geometry_jet(r)
    dV = vol.g @ dr
    dA = area.g @ dr
    low, _ = cert.enc.geometry_jet(r - dr / 2)
    high, _ = cert.enc.geometry_jet(r + dr / 2)
    finite = high.v - low.v
    proof = json.loads((root / "certificate.json").read_bytes())
    assert proof["gradient_norm_upper"] < 8e-8
    # The certificate bounds the residual of S+eta*A-lambda*V at
    # this numerical center, with its nominal multiplier below 52.
    err = cert.ri("8e-8") * hs + 52 * abs(dV)
    gain = iv.mpf(proof["gain_interval"])
    reader = measure.HERE / "reference_read.json"
    sample = json.loads(reader.read_bytes())["cases"][str(2**38)]["readings"][
        "equilibrium"
    ]
    # I⊗log(rho_R), log(rho_L)⊗I and the parity scalar have zero trace
    # against D. Thus the connected charge is the full entropy derivative.
    # Cover binary64 preparation rounding far below the statistical errors.
    measured = iv.mpf(sample["interval"]) + iv.mpf([-1e-11, 1e-11])
    prediction = (-measured + iv.mpf([-1, 1]) * err) / gain

    def bounds(x):
        return [
            float(np.nextafter(float(x.a), -np.inf)),
            float(np.nextafter(float(x.b), np.inf)),
        ]

    flat = []
    for tape in (r - dr / 2, r + dr / 2):
        square = measure.station.square(np.array([float(x.mid) for x in tape]))
        geometry = measure.station.distance_check(square)
        d = np.array(geometry["projected"])
        gram = (d[0, 1:, None] ** 2 + d[None, 0, 1:] ** 2 - d[1:, 1:] ** 2) / 2
        sites = np.vstack((np.zeros(3), np.linalg.cholesky(gram)))
        fitted = np.linalg.norm(sites[:, None] - sites[None, :], axis=2)
        flat.append(
            dict(
                stations=sites.tolist(),
                smallest_gram_eigenvalue=np.linalg.eigvalsh(gram)[0],
                maximum_symmetric_range_residual=abs(fitted - d).max(),
                directed_log_error_floor=geometry["necessary_log_error_floor"],
            )
        )
    return dict(
        reader_sha256=sha256(reader.read_bytes()).hexdigest(),
        propagation_bank_sha256=sha256(bank.read_bytes()).hexdigest(),
        equilibrium_certificate_sha256=sha256(
            (root / "certificate.json").read_bytes()
        ).hexdigest(),
        derivative_direction="from lower to upper source endpoint, full displacement",
        volume_derivative=bounds(dV),
        ordinary_gain=bounds(gain),
        stationarity_plus_volume_error=float(np.nextafter(float(err.b), np.inf)),
        area_derivative_from_source=bounds(prediction),
        area_derivative_from_radar=bounds(dA),
        finite_endpoint_area_change=bounds(finite),
        flat_realizations_of_symmetric_ranges=flat,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    path = measure.HERE / "source_area.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), result)
    print(result)
