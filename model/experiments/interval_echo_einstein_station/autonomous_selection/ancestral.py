"""Test whether compact stationary sector biases dilute on a fixed ancestral patch."""
import argparse
from collections import Counter
from itertools import combinations
import json
from math import comb
from pathlib import Path

import numpy as np

from model.executable_model import federation
from model.experiments.interval_echo_space.shell_source import laplacian
from ..shape_stationarity.geometry import selected_shapes
from .measure import SEED, local_modes
from .open_host import SHOTS


def run():
    parent = sorted(dict(selected_shapes(federation(0)))['cycle5'])
    rows = []
    for level in range(3):
        car = federation(level)
        cells = [c for c in range(car.cells) if c//(4**level) in parent]
        cell_pairs = list(combinations(cells, 2))
        per_side = len(cell_pairs)//2
        host_size = comb(car.n,2)
        max_residual = 0.
        participation = {}
        for site in cells:
            offset = 12*site
            gateways = tuple(np.flatnonzero(car.degree[offset:offset+12]==6))
            modes = local_modes(gateways)
            local = [(a-offset,b-offset) for a,b in car.seams if a//12==site==b//12]
            residual = (laplacian(12,local)/2) @ modes - 3*modes
            max_residual = max(max_residual,float(np.max(abs(modes[list(gateways)]))),
                               float(np.max(np.linalg.norm(residual,axis=0))))
            p = (abs(modes)**2).sum(axis=1)
            participation[site] = p
        bins = Counter()
        for sign,pairs in ((1,cell_pairs[:per_side]),(-1,cell_pairs[-per_side:])):
            for a,b in pairs:
                word_weights = np.outer(participation[a],participation[b]).ravel()
                bins.update((sign,float(w)) for w in word_weights.round(12) if w > 1e-12)
        classes = sorted(bins)
        signs = np.array([c[0] for c in classes])
        signed_w = np.array([s*w for s,w in classes])
        multiplicities = np.array([bins[c] for c in classes])
        tangent = multiplicities*signed_w/(4*per_side)
        np.testing.assert_allclose(tangent[signs>0].sum(),1,atol=1e-10)
        np.testing.assert_allclose(tangent[signs<0].sum(),-1,atol=1e-10)
        gradient = per_side*4/host_size
        prepared = []
        for s, alpha in enumerate((-.75,0.,.75)):
            prob = multiplicities*(1+alpha*signed_w)/host_size
            counts = np.random.default_rng([SEED,level,s,880]).multinomial(SHOTS,np.r_[prob,1-sum(prob)])
            read = counts[:-1]/SHOTS
            prepared.append(dict(alpha=alpha, selected_class_probabilities=[
                                     float(prob[signs>0].sum()),float(prob[signs<0].sum())],
                                 transferred_probability=abs(alpha)*gradient,
                                 signed_entropy_slope=float(-tangent @ np.log(prob/multiplicities)),
                                 counts=counts.tolist(), measured_slope=float(-tangent @ np.log(read/multiplicities)),
                                 standard_error=float(np.sqrt(np.sum(tangent**2/read)/SHOTS))))
        rows.append(dict(host_cells=car.cells, region_cells=len(cells), addressed_sites=cells,
                         pairs_per_side=per_side, rank_per_side=4*per_side, host_pair_dimension=host_size,
                         single_cell_gateway_and_eigenvector_residual=max_residual,
                         class_sign_and_participation=classes, multiplicities=multiplicities.tolist(),
                         preparations=prepared))
    return dict(shots_per_terminal_read=SHOTS, seed=SEED, ancestor=parent, rows=rows)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,type=Path)
    data = run()
    with parser.parse_args().output.open('x') as f:
        json.dump(data,f,separators=(',',':'),allow_nan=False)
        f.write('\n')
