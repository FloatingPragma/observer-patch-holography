"""All-seams-open two-load counterexample, with inclusive regional endpoint observations."""
import argparse
from collections import defaultdict
from math import comb
import json
from pathlib import Path

import numpy as np

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from ..equilibrium_refinement.measure import apparatus, seam
from ..shape_stationarity.geometry import selected_shapes
from .measure import SEED, stationary_pair_tangent

SHOTS = 2**34
HERE = Path(__file__).resolve().parent


def full_stencil(car, words, waves):
    h = seam()
    positive = np.flatnonzero(np.max(abs(waves),axis=1) > 1e-13)
    response = defaultdict(lambda: np.zeros(waves.shape[1], complex))
    for i in positive:
        a,b = words[i]
        active = int(car.degree[a] + car.degree[b]) - 2*(b in car.neighbours[a])
        diagonal = float(h[1,1]*active + 4*candidate.KAPPA*(12-int(car.degree[a])-int(car.degree[b])))
        response[a,b] += diagonal * waves[i]
        for src,other in ((a,b),(b,a)):
            for dest in car.neighbours[src]:
                if dest != other:
                    response[tuple(sorted((dest,other)))] += h[1,2]*waves[i]
    for i in positive:
        response[words[i]] -= 390*waves[i]
    return float(np.linalg.norm(np.stack(list(response.values())),axis=0).max())


def run():
    reports = []
    for level in range(3):
        car = federation(level)
        host_words = comb(car.n,2)
        for j, (name,cells) in enumerate(selected_shapes(car)):
            _, _, cut, words = apparatus(car,cells)
            direction,a,b = stationary_pair_tangent(car,cells,words)
            plus = np.flatnonzero(direction > 1e-13)
            minus = np.flatnonzero(direction < -1e-13)
            if len(plus) != 36 or len(minus) != 36:
                raise ArithmeticError("Two equal-size sharp interior sets are required")
            np.testing.assert_allclose(direction[plus],1/36,atol=5e-15)
            np.testing.assert_allclose(direction[minus],-1/36,atol=5e-15)
            assert all(set(words[i]).isdisjoint(cut) for i in np.r_[plus,minus])
            residual = full_stencil(car,words,np.column_stack([a,b]))
            if residual > 1e-10:
                raise ArithmeticError("Free full-host evolution leaves the proposed eigenspace")
            rows = []
            for s,alpha in enumerate((-.75,0.,.75)):
                bins = np.array([(len(plus)+4*alpha)/host_words,(len(minus)-4*alpha)/host_words])
                counts = np.random.default_rng([SEED,level,j,s,238]).multinomial(
                    SHOTS,np.r_[bins,1-sum(bins)])
                measured = float(-np.log(counts[0]/counts[1]))
                se = float(np.sqrt(1/counts[0]+1/counts[1]))
                rows.append(dict(alpha=alpha, inclusive_counts=counts.tolist(), class_probabilities=bins.tolist(),
                                 signed_slope_nats=float(-np.log(bins[0]/bins[1])),
                                 measured_slope_nats=measured, standard_error=se))
            reports.append(dict(board_cells=car.cells, host_pair_words=host_words, shape=name,
                                source_cells=sorted(cells)[:3], region_cells=sorted(cells),
                                affected_sharp_words=[[words[i] for i in ids] for ids in (plus,minus)],
                                full_generator_eigenvector_residual=residual,
                                chance_both_residents_in_region=len(words)/host_words, preparations=rows))
    return dict(shots_per_terminal_read=SHOTS, seed=SEED, rows=reports)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,type=Path)
    data=run()
    with parser.parse_args().output.open('x') as f:
        json.dump(data,f,separators=(',',':'),allow_nan=False)
        f.write('\n')
