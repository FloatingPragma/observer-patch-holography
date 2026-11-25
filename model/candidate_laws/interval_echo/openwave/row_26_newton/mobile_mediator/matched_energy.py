"""Orientational reversal with identical native static source energy on equal-degree pins."""

import argparse
from itertools import product
import json

import numpy as np

from . import measure


def regular_pair(board, gate):
    first = min(a for a in board.neighbours[gate] if board.degree[a] == 5)
    # Orient the adjacent pair from the near to the far side of the gateway.
    second = min(b for b in board.neighbours[first] if b//12 == gate//12
                 and board.degree[b] == 5 and b not in board.neighbours[gate])
    return first, second


def source_energy(board, pair, reversed_source):
    blank = np.ones(board.n, dtype=int)
    blank[list(pair)] += [1,-1] if reversed_source else [-1,1]
    seams = np.array(board.seams)
    contrast = blank[seams[:,0]] - blank[seams[:,1]]
    return float(8*(contrast @ contrast) + np.count_nonzero(contrast)/2
                 + 8*(6-board.degree) @ (blank**2-1))


def run():
    board = measure.federation(0)
    panels = []
    for (distance, path, original), q, reverse in product(
            tuple(measure.sources(board)), (1,-1), (False,True)):
        a, b = regular_pair(board,original[0]), regular_pair(board,original[2])
        pins = (*a,*b)
        graph,h,va,vb = measure.operators(board,pins,8,q,reverse)
        assert not np.any(va*vb)
        assert all(np.min(h.diagonal()+i*va+j*vb-graph.sum(axis=1)) >= 0
                   for i, j in measure.SETTINGS)
        heat = measure.thermal_trace(h,va,vb)
        remainder = 16*len(va)**2*1e-72
        panels.append(dict(distance=distance,path=path,pins=pins,mobile_charge=q,
                           reversed_receiver=reverse,
                           first_source_energy=source_energy(board,a,False),
                           second_source_energy=source_energy(board,b,reverse),
                           first_mixed=measure.first_mixed_moment(graph,va,vb),
                           heat_trace_beta_eighth=heat,
                           potential_per_activity_beta_eighth_z1e_minus72=[
                               float(np.nextafter(-8*heat[1]-remainder, -np.inf)),
                               float(np.nextafter(-8*heat[0]+remainder, np.inf))]))
    return dict(kappa=8, mobile_positions=236, panels=panels)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    output = json.dumps(run(),indent=2,allow_nan=False) + "\n"
    if args.write:
        (measure.HERE/"matched_energy.json").write_text(output)
    else:
        print(output,end="")
