"""Empty-cell compatible extension of the complete seven-input native chain."""
from fractions import Fraction as F
from functools import lru_cache
import numpy as np
import scipy.sparse as sp
from scipy.linalg import eigh
from model.openwave_strong_step_scaling.source import chain as original_chain
import sys
from pathlib import Path
from model.openwave_hadron_refinement_flow.parameters import degree_coefficients


@lru_cache(None)
def chain(length,flow=False):
    d = original_chain(length)
    p = d['parameters']
    h = d['H'].copy()
    correction = sp.csr_matrix(h.shape)
    for x in range(length):
        excess = int(x > 0) + int(x + 1 < length) - 1
        if excess:
            correction += excess * (float(p['spin_gradient']) * d['angular']((x,), True)
                                     + float(p['eta']) * d['angular']((x,), False)
                                     + float(p['t']) * sp.diags(d['n'][:, x], dtype=np.int64))
    h -= correction
    if flow:
        dU,dm=degree_coefficients()
        for x in range(1,length-1):
            n=d['n'][:,x]
            h+=sp.diags(float(dm)*n+float(dU)*2*n*(3-n)/3)
    if length == 2:
        assert correction.nnz == 0
    return {**d, 'H': h, 'vacuum_self_correction': correction,'degree_flow':flow}


def probe(length):
    d = chain(length)
    n = np.flatnonzero((d['nu'] == 1) & (d['ns'] == 2))
    p = np.flatnonzero((d['nu'] == 2) & (d['ns'] == 2))
    hn = d['H'][n][:, n].toarray()
    e, v = eigh(hn)
    psi = v[:, 0]
    phi = d['Ip'][p][:, n] @ psi
    gv = np.linalg.norm(phi)
    assert abs(gv - 1) < 2e-8
    phi /= gv
    ga = float(phi @ (d['Jp'][p][:, n] @ psi))
    ep = e[0] - float(d['parameters']['delta'])
    assert np.linalg.norm(d['H'][p][:, p] @ phi - ep * phi) < 2e-8
    return {'cells': length, 'dimension': len(d['labels']), 'energy_source': float(e[0]),
            'energy_GeV': float(e[0]) * .00051099895069, 'gV': float(gv), 'gA': ga,
            'gap_source': float(e[1] - e[0]),
            'residual': float(np.linalg.norm(hn @ psi - e[0] * psi))}


if __name__ == '__main__':
    import json
    for length in (2, 3, 4):
        print(json.dumps(probe(length)), flush=True)
