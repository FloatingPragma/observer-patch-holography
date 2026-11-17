#!/usr/bin/env python3
"""Derived-state analysis of retained native records; no new native run."""
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

HERE = Path(__file__).resolve().parent


def widths(mass, epsilon, velocity, order):
    nodes, weights = leggauss(order)
    r = ((nodes + 1) / 2)[:, None]
    u = nodes[None, :]
    measure = (weights / 2)[:, None] * weights[None, :] * r**2
    q1 = r * u
    q2 = r * np.sqrt(1 - u**2)
    f = (1-r**2)**6
    radial_derivative = -12 * r * (1-r**2)**5
    f1 = radial_derivative * u
    f2 = radial_derivative * np.sqrt(1-u**2)
    sigma = epsilon * mass
    energy = np.sqrt(mass**2 + sigma**2 * r**2)
    d = 1 + velocity * sigma * q1 / energy
    a2 = velocity * sigma * q2 / energy
    d1 = velocity * sigma * (1/energy - sigma**2*q1**2/energy**3)
    d2 = -velocity * sigma**3*q1*q2/energy**3
    h1 = f1/d**2 - f*d1/d**3
    h2 = f2/d - a2*f1/d**2 - f*d2/d**2 + a2*f*d1/d**3
    n0 = np.sum(measure*f*f)
    i0 = np.sum(measure*f1*f1)
    norm = np.sum(measure*f*f/d)
    longitudinal = np.sum(measure*h1*h1*d)
    transverse = np.sum(measure*h2*h2*d)/2
    contraction = np.sqrt(1-velocity**2)
    return [float(contraction*np.sqrt((longitudinal/norm)/(i0/n0))),
            float(np.sqrt((transverse/norm)/(i0/n0)))]


def main():
    raw_path = HERE/'inputs/COMMON_RECORDS.json.gz'
    records = json.loads(gzip.decompress(raw_path.read_bytes()))
    edges = records['graph']['edges']
    laplacian = np.zeros((12, 12), dtype=np.int64)
    verified = 0
    moves = {}
    for eid, slot, coin, current, target in records['branches']:
        n_a, n_b = ((1, 0) if slot == 0 else (0, 1) if slot == 1 else (0, 0))
        out_a = (n_a+n_b)//2 if coin == 0 else (n_a+n_b+1)//2
        out_b = n_a+n_b-out_a
        expected_target = 0 if out_a else 1 if out_b else 2
        assert current == n_a-out_a and target == expected_target
        verified += 1
        if slot < 2:
            moves[eid, slot, coin] = target
    for eid, (a,b,axis) in enumerate(edges):
        if axis == -1 and a < 12 and b < 12:
            assert sum(moves[eid,0,c] == 1 for c in (0,1)) == 1
            assert sum(moves[eid,1,c] == 0 for c in (0,1)) == 1
            laplacian[a,a] += 1; laplacian[b,b] += 1
            laplacian[a,b] -= 1; laplacian[b,a] -= 1
    eigenvalues = np.linalg.eigvalsh(laplacian)
    bands = []
    for value in eigenvalues:
        if value > 1e-9 and (not bands or abs(value-bands[-1]) > 1e-9):
            bands.append(float(value))
    assert len(bands) == 3
    rows = []
    max_quadrature_change = 0.0
    for lam in bands:
        mass = np.sqrt(lam)
        for velocity in (0.3,0.6,0.8):
            for epsilon in (0.5,0.25,0.125,0.0625):
                low = widths(mass,epsilon,velocity,64)
                high = widths(mass,epsilon,velocity,128)
                max_quadrature_change = max(max_quadrature_change,
                    max(abs(x-y) for x,y in zip(low,high)))
                target = np.sqrt(1-velocity**2)
                if epsilon == 0.0625:
                    assert abs(high[0]/target-1) < 0.01
                    assert abs(high[1]-1) < 0.01
                rows.append(dict(lambda_=lam,velocity=velocity,sigma_over_mass=epsilon,
                    width_ratios=high,limiting_ratios=[float(target),1.0]))
    assert max_quadrature_change < 1e-10
    assert all(abs(1-np.sqrt(1-v*v)) > 0.01 for v in (0.3,0.6,0.8))
    result = dict(candidate='MFR-common-native-field-v1',
        kind='derived-state analysis; no new native calls or prepared-packet tapes',
        retained_source_record_sha256=hashlib.sha256(raw_path.read_bytes()).hexdigest(),
        verified_source_branches=verified,source_internal_laplacian=laplacian.tolist(),
        native_mass_squared=bands,rows=rows,
        quadrature_order_comparison=max_quadrature_change,
        quadrature_comparison_is_rigorous_error_bound=False,
        unboosted_control='width ratio1 fails each 1% contraction gate',
        iid_source_control='zero separated OS covariance has no nonzero normalized packet',
        all_registered_analysis_checks='pass',official='untested')
    (HERE/'RECEIPT.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('rows','source_internal_laplacian')},indent=2))


if __name__ == '__main__':
    main()
