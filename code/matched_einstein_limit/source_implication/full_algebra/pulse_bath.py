"""Connected probe relaxation and its fixed-volume response on the enrolled 80-cell source."""

import argparse
import json
from pathlib import Path

import numpy as np
import sympy as sp

from . import certificate as cert
from . import native_precision, station


def commutator_spectrum(h):
    """Exact Gram eigenvalues in a traceless Hilbert--Schmidt orthonormal Pauli basis."""
    blocks = h.tobsr(blocksize=(4, 4)).data
    reduced = blocks - blocks[:, 0, 0, None, None] * np.eye(4, dtype=int)
    matrices, counts = np.unique(reduced.reshape(-1, 16), axis=0, return_counts=True)
    # Each doubled basis matrix has Gaussian-integer entries. Use integer dot products.
    paulis = 2 * station.iso.basis(False)
    re, im = paulis.real.astype(int), paulis.imag.astype(int)
    assert np.array_equal(paulis, re + 1j*im)
    numerator = np.zeros((15, 15), dtype=np.int64)
    imaginary = np.zeros_like(numerator)
    for flat, count in zip(matrices, counts):
        a = flat.reshape(4, 4)
        u = (a @ re - re @ a).reshape(15, 16)
        v = (a @ im - im @ a).reshape(15, 16)
        assert max(abs(u).max(), abs(v).max()) < 64 and len(blocks) < 2**30
        numerator += count*(u @ u.T + v @ v.T)
        imaginary += count*(u @ v.T - v @ u.T)
    assert not np.any(imaginary)
    eigenvalues = (sp.Matrix(numerator)/4).eigenvals()
    assert min(eigenvalues) >= 256
    return {str(value): int(n) for value, n in sorted(eigenvalues.items())}


def generator_checks():
    board = station.iso.pairing.formation.federation(1)
    route = station.iso.pairing.formation.region(board)
    data = []
    for charge in (-1, 1):
        for odd in (False, True):
            h = native_precision.generator(board, route, charge, odd)
            assert h.shape == (4*board.n, 4*board.n) and (h-h.T).nnz == 0
            row_bound = int(abs(h).sum(axis=1).max())
            assert row_bound < 48
            data.append(dict(charge=charge, odd=odd, integer_row_bound=row_bound,
                             commutator_gram=commutator_spectrum(h)))
    # Changing the contrast gain preserves the graph and both parity projections.
    # In the even sector all changed-law marked blocks are scalar in the heads.
    for charge in (-1, 1):
        local = station.iso.ranging.parity_blocks(charge, 0, 1)
        np.testing.assert_array_equal(local, local[..., 0, 0, None, None]*np.eye(2))

    return dict(zero_gain_even_head_channel='identity', ports=board.n, dwell='1/16', lazy_probability='1/2',
                squared_norm_loss_floor=str(sp.Rational(1, 64*board.n)), sectors=data)


def matched_readout():
    cert.iv.dps = 30
    nums, unit, error = native_precision.exact_effects()
    eps = sp.Rational(error)
    effects = cert.rational_effects(nums, unit)
    gain, target = cert.calibration(effects, eps)
    even = np.array([cert.ri_band(pair[0].trace()/4, eps) for pair in effects], dtype=object)
    odd = np.array([cert.ri_band(pair[1].trace()/4, eps) for pair in effects], dtype=object)
    dr = odd-even

    def jets(p):
        return cert.enc.geometry_jet(even + dr*p)

    # Monotonicity on the whole physical parity interval, including either pure parity.
    derivatives = []
    for k in range(32):
        lo, hi = cert.ri(sp.Rational(k, 32)), cert.ri(sp.Rational(k+1, 32))
        _, volume = jets(cert.iv.mpf([lo.a, hi.b]))
        derivatives.append(volume.g @ dr)
    assert all(v.b < cert.ri('-.075') for v in derivatives)
    lo, hi = cert.ri('.10750779'), cert.ri('.10750781')
    assert jets(lo)[1].v.a > target.b and jets(hi)[1].v.b < target.a
    area, volume = jets(cert.iv.mpf([lo.a, hi.b]))
    # YZ and ZZ are measured Pauli products; both perturbations are supported in the even block.
    _, _, y, z = map(cert.prior.sym_matrix, station.phase.AXES)
    directions = [sp.kronecker_product(left, z)/2 for left in (y, z)]
    slopes = np.array([[cert.ri_band(cert.inner(pair[0], d), 2*eps) for d in directions]
                       for pair in effects], dtype=object)
    az, vz = area.g @ slopes, volume.g @ slopes
    assert vz[1].a > cert.ri('.035')
    adjust = -vz[0]/vz[1]
    response = gain*(az[0]+az[1]*adjust)
    assert response.b < cert.ri('-.004') and response.a > cert.ri('-.006')

    # The separate local-existence certificate encloses the true equilibrium within 2e-7 HS.
    from ..material_clocks.measure import density_center
    center, _, _, radius = density_center()
    zz = sp.kronecker_product(z, z)
    mean = sum(cert.inner(cert.prior.sym_matrix(block), zz) for block in center)
    # density_center converts the normalized rational matrix to binary64 once; each
    # real diagonal error is below 2^-50 since each absolute entry is below one.
    expect = cert.ri_band(mean, 3*sp.Rational(radius)+8*sp.Rational(1, 2**50))
    assert expect.a > cert.ri('.4608') and expect.b < cert.ri('.4610')
    assert sp.Rational('.14')*sp.Rational('.4608')**2 > sp.Rational('.0297')
    return dict(odd_weight_interval=['.10750779', '.10750781'], calibrated_volume=str(target),
                volume_strictly_decreasing=True, volume_control=str(vz[1]), matching_coefficient=str(adjust),
                fixed_volume_objective_slope=str(response), optimizer_joint_ZZ=str(expect),
                bath_joint_ZZ=0, feasible_gap_lower_nats='0.0297')


def run():
    return dict(native_collision=generator_checks(), matched_measurement=matched_readout())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    output = json.dumps(run(), indent=2, allow_nan=False)+'\n'
    if args.write:
        Path(__file__).with_suffix('.json').write_text(output)
    else:
        print(output, end='')
