"""Held fifth station at the independently selected four-station equilibrium."""

from itertools import permutations
import json
from pathlib import Path

import numpy as np
from scipy.linalg import helmert
from scipy.sparse.linalg import expm_multiply

from ..geometry_gain import pairing, ranging, simplex_area_volume


def leg_kernels(h, ports):
    indices = np.array([4*p+i for p in ports for i in range(4)])
    launch = np.zeros((h.shape[0], len(indices)), complex)
    launch[indices] = np.eye(len(indices))
    waves = expm_multiply(-1j*h, launch, start=0., stop=6., num=13, traceA=0)
    return waves[:, indices].reshape(13,len(ports),4,len(ports),4).transpose(0,1,3,2,4)


def fit(ranges):
    n = len(ranges)
    other = ~np.eye(n,dtype=bool)
    logs = np.zeros_like(ranges)
    logs[other] = np.log(ranges[other])
    gain = (logs-logs.T).mean(axis=1)
    adjusted = logs-gain[:,None]
    distances = np.exp((adjusted+adjusted.T)/2)
    np.fill_diagonal(distances,0)
    contrast = helmert(n)
    gram = -contrast @ (distances**2) @ contrast.T/2
    eigen = np.linalg.eigvalsh(gram)[::-1]
    return dict(gains=np.exp(gain).tolist(), distances=distances.tolist(),
                largest_directional_log_residual=float(np.max(abs(adjusted-adjusted.T))/2),
                centered_eigenvalues_descending=eigen.tolist(), last_to_first=float(eigen[-1]/eigen[0]),
                squared_distance_error_floor=float(2*eigen[-1]/(n-1)))


def run():
    receipt = json.loads((Path(__file__).resolve().parents[1]/"full_algebra/results.json").read_text())
    state = np.asarray(receipt["x_real"])+1j*np.asarray(receipt["x_imag"])
    board = pairing.formation.federation(1)
    source = pairing.formation.region(board)
    ports = tuple(source)+(5,)
    assert len(set(ports))==5
    ranges = np.zeros((5,5))
    for odd, block in enumerate(state):
        for charge in (-1,1):
            h = ranging.generator(board,source,charge,8,bool(odd))
            kernels = leg_kernels(h,ports)
            for p,q in permutations(range(5),2):
                out, ret = ranging.first_hits(kernels,p,q), ranging.first_hits(kernels,q,p)
                # A completed return advances the read below the cap; all other histories contribute 3.
                effect = 3*np.eye(4,dtype=complex)
                for n in range(2,12):
                    for t in range(1,n):
                        amp = ret[n-t] @ out[t]
                        effect += (n/4-3)*(amp.conj().T @ amp)
                ranges[p,q] += np.trace(block@effect).real/2
    geometry = simplex_area_volume(ranges[:4,:4])
    np.testing.assert_allclose(geometry,receipt["reports"]["6"]["area_volume"],rtol=0,atol=5e-12)
    return dict(ports=ports, cells=board.cells, tick=.5, horizon=6, eta=receipt["ordinary_eta"],
                range_matrix=ranges.tolist(), original_area_volume=np.asarray(geometry).tolist(),
                original_four=fit(ranges[:4,:4]), extended_five=fit(ranges))


if __name__=="__main__":
    print(json.dumps(run(),indent=2))
