"""Stationary source of an equally scheduled, head-reusing version of the four-clock survey."""

import argparse
import json
from pathlib import Path

import numpy as np
import sympy as sp

from . import certificate as c
from . import global_survey, station


def sector_center(channels, error, operators):
    scale = 2**35
    assert abs(channels).max() < 2
    # Round 24 small binary64 matrices separately on an exactly represented binary grid.
    re = np.rint(channels.real*scale).astype(np.int64).sum(axis=0)
    im = np.rint(channels.imag*scale).astype(np.int64).sum(axis=0)
    channel = (sp.Matrix(re) + sp.I*sp.Matrix(im))/(len(channels)*scale)
    eps = sp.Rational(float(error)) + sp.Rational(64, scale)
    basis = sp.Matrix.hstack(*(op.reshape(16, 1) for op in operators))
    uniform = sp.eye(4)/4
    jac = (basis.H * channel * basis).applyfunc(sp.re)
    shift = (basis.H * channel * uniform.reshape(16,1)).applyfunc(sp.re)
    upper = sp.Rational('.918')
    _, pivots = (upper**2*sp.eye(15)-jac.T*jac).LDLdecomposition(hermitian=False)
    assert min(pivots.diagonal()) > 0 and upper+2*eps < sp.Rational('.92')
    proposal = np.linalg.solve(np.eye(15)-np.array(jac, dtype=float), np.array(shift, dtype=float).ravel())
    coords = sp.Matrix([sp.Rational(round(v*scale), scale) for v in proposal])
    rho = uniform + (basis*coords).reshape(4,4)
    assert rho == rho.H and rho.trace() == 1 and global_survey.positive_definite(rho-sp.Rational('.19')*sp.eye(4))
    residual = (sp.eye(15)-jac)*coords-shift
    assert sum(v*v for v in residual) < sp.Rational(1, 10**16)
    radius = (eps + sp.Rational(1, 10**8))/(1-sp.Rational('.92'))
    assert radius < sp.Rational('.000025')
    return rho, radius


def verify(bank, channel_error):
    c.iv.dps = 30
    ops = [c.prior.sym_matrix(h) for h in station.iso.basis(False)]
    items = [sector_center(bank[:,b].reshape(24,16,16), channel_error, ops) for b in (0,1)]
    nums, unit, native_error = c.exact_effects()
    error = sp.Rational(native_error)
    effects = c.rational_effects(nums, unit)
    global_survey.certify_ranges(effects, error)
    gain, target = c.calibration(effects, error)
    fractions = [(sp.Rational(hi-lo, 10000)) for lo, hi in global_survey.BOUNDS]
    range_centers = [np.array([c.ri(c.inner(pair[b],rho)) for pair in effects], dtype=object)
                     for b, (rho, _) in enumerate(items)]
    half_widths = np.array([c.ri(error+width*max(x[1] for x in items)) for width in fractions],
                           dtype=object)
    centered_drift = range_centers[1]-range_centers[0]
    unit = c.iv.mpf([-1,1])
    drift = centered_drift + half_widths*(2*unit)

    def jets(p):
        # First interpolate the nominal means, then use the convex mixture's uniform error.
        center = range_centers[0]+centered_drift*p
        return c.enc.geometry_jet(center+half_widths*unit)

    def endpoint_volume(p):
        center = range_centers[0]+centered_drift*p
        _, nominal = c.enc.geometry_jet(center)
        _, enclosed = jets(p)
        return nominal.v + unit*sum(abs(j)*w for j,w in zip(enclosed.g,half_widths))

    for i in range(32):
        p = c.iv.mpf([i, i+1])/32
        assert (jets(p)[1].g @ drift).b < c.ri('-.075')
    low, high = c.ri('.134'), c.ri('.1346')
    endpoints = endpoint_volume(low), endpoint_volume(high)
    assert endpoints[0].a > target.b and endpoints[1].b < target.a, str(endpoints)
    p = c.iv.mpf([low.a, high.b])
    ar, vr = jets(p)
    chosen = [ops[i] for i in (2, 10, 6, 14)]  # IZ/2, YZ/2, XZ/2, ZZ/2 in the even block
    directions = np.array([[c.ri_band(c.inner(pair[0],op), 2*error) for op in chosen]
                           for pair in effects], dtype=object)
    az, vz = ar.g @ directions, vr.g @ directions
    assert vz[-1].a > c.ri('.03')
    even, radius = items[0]
    grad, err = c.prior.entropy_coordinates(even,radius,sp.Rational('.18'),chosen)
    ds = [c.ri_band(v,2*err) for v in grad]
    response, entropy_slopes, area_slopes = [], [], []
    for i in (0,1,2):
        adj = -vz[i]/vz[-1]
        entropy_slopes.append(ds[i]+adj*ds[-1])
        area_slopes.append(az[i]+adj*az[-1])
        response.append(entropy_slopes[-1]+gain*area_slopes[-1])
    assert entropy_slopes[0].b < 0 and entropy_slopes[1].a > 0
    assert area_slopes[0].b < 0 and area_slopes[1].b < 0 and area_slopes[2].a > 0
    gains = [-s/a for s,a in zip(entropy_slopes,area_slopes)]
    assert gains[1].b < gains[2].a, [str(z) for z in gains]
    # The area-preserving XZ transfer includes -(A_XZ/A_YZ) times the YZ transfer.
    invisible_entropy = area_slopes[2]*(gains[1]-gains[2])
    assert invisible_entropy.b < c.ri('-.00033')
    assert response[0].b < c.ri('-.06'), str(response[0])
    assert response[1].b < 0, str(response[1])

    # A Pauli expansion with sum of absolute nonidentity coefficients below
    # the identity weight is a strictly interior mixture of native product atoms.
    moments = [np.array([c.ri_band(2*c.inner(rho,op),2*rad) for op in ops],dtype=object)
               for rho,rad in items]
    slack = [(1-2*p)-sum(abs(v) for v in (moments[0]*(1-p)-moments[1]*p)),
             2*p-sum(abs(v) for v in (moments[1]*(2*p)))]
    assert all(value.a > c.ri('.20') for value in slack)
    zz = (1-p)*moments[0][14]+p*moments[1][14]
    return dict(conditional_HS_contraction_upper='.92', conditional_state_HS_error=str(c.ri(radius)),
                matched_odd_weight=['.134','.1346'], target_volume=str(target), joint_ZZ=str(zz),
                fixed_volume_marginal_coin_slope=str(response[0]), fixed_volume_correlation_only_slope=str(response[1]),
                entropy_slopes=[str(v) for v in entropy_slopes], area_slopes=[str(v) for v in area_slopes],
                would_require_gains=[str(g) for g in gains], fixed_area_and_volume_entropy_slope=str(invisible_entropy),
                native_product_frame_slack=[str(s) for s in slack])


def run():
    from .reuse.measure import channels
    return verify(*channels())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    output = json.dumps(run(), indent=2, allow_nan=False)+'\n'
    if args.write:
        Path(__file__).with_suffix('.json').write_text(output)
    else:
        print(output, end='')
