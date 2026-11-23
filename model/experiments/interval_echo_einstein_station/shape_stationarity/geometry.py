"""Select equal-size regions by seam incidence and independently measure their apertures."""
from functools import lru_cache

import networkx as nx
import numpy as np
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from model.experiments.contact import radar
from model.experiments.interval_echo_einstein_station.regional_geometry.measure import hand_times

SHOTS = 2**30
SEED = 26112362


def selected_shapes(board):
    """Use the lexicographically first eligible five-cycle and neighbouring path at home."""
    graph = nx.Graph()
    graph.add_edges_from((a//12, b//12) for a, b in board.external)
    walks = [(0,)]
    for _ in range(4):
        walks = [walk+(n,) for walk in walks for n in graph[walk[-1]] if n not in walk]
    cycles = [tuple(sorted(w)) for w in walks
              if graph.has_edge(w[-1], 0) and graph.subgraph(w).number_of_edges()==5]
    cycle = set(min(cycles))
    paths = [tuple(sorted(w)) for w in walks
             if graph.subgraph(w).number_of_edges()==4 and len(set(w)&cycle)==4]
    path = set(min(paths))
    # The same new address attaches by one seam in both shapes. On the tree
    # it attaches to an endpoint, giving a held six-cell path.
    common_extension = min(n for n in graph if n not in cycle|path
                           and len(set(graph[n])&cycle)==1 and len(set(graph[n])&path)==1
                           and max(dict(graph.subgraph(path|{n}).degree).values())==2)
    return [("cycle5", cycle), ("path5", path), ("cycle5_spur", cycle|{common_extension}),
            ("path6", path|{common_extension})]


def meter(h, region, hand, rng):
    """Release every interior sharp port and measure endpoint membership, retaining quiet."""
    dimension = h.shape[0]
    outside = np.ones(dimension, dtype=bool)
    outside[region] = False
    preparations = np.zeros((dimension, len(region)))
    preparations[region, np.arange(len(region))] = 1
    states = expm_multiply(-1j*h, preparations, start=0., stop=2*radar.DT, num=3,
                           traceA=-1j*h.diagonal().sum())
    crossing = np.sum(abs(states[1:,outside])**2, axis=(1,2))/len(region)
    area = 4*len(region)*crossing/hand**2
    successes = rng.binomial(SHOTS, crossing)
    read = successes/SHOTS
    observed_areas = 4*len(region)*read/hand**2
    # A real symmetric flight has an even expansion in time. Extrapolate
    # the two finite appointments to remove the leading t² aperture bias.
    weights = np.array([hand[1]**2, -hand[0]**2])/(hand[1]**2-hand[0]**2)
    observed_se = np.sqrt(np.sum((weights*4*len(region)/hand**2)**2*read*(1-read)/SHOTS))
    return dict(hand=hand.tolist(), cross_probability=crossing.tolist(),
                finite_areas=area.tolist(), extrapolated_area=float(weights@area),
                crossing_counts=successes.tolist(), trials_per_appointment=SHOTS,
                sampled_extrapolated_area=float(weights@observed_areas),
                sampled_extrapolated_standard_error=float(observed_se))


@lru_cache(maxsize=1)
def geometry():
    """Blank probe calibrated on separate copies of the addressed screen for each shape."""
    hand = hand_times()[[1,2]]
    rng = np.random.default_rng(SEED)
    records = []
    for level in range(3):
        board = federation(level)
        h = candidate.generator(board).tocsr()
        for name, cells in selected_shapes(board):
            region = np.array([12*cell+offset for cell in sorted(cells) for offset in range(12)])
            ends = sorted((a,b) if a//12 in cells else (b,a) for a,b in board.external
                          if (a//12 in cells)!=(b//12 in cells))
            interior_ports = [a for a,_ in ends]
            if len(set(interior_ports))!=len(ends):
                raise ArithmeticError("multiple crossing modes share an interior gateway")
            flight = meter(h,region,hand,rng)
            area = float(len(ends))
            relative_error = abs(flight["extrapolated_area"]/area-1)
            if relative_error >= .05:
                raise ArithmeticError("finite-hand aperture failed 5% agreement with independent seam count")
            records.append(dict(board_level=level, cells=sorted(cells), shape=name,
                                accessible_volume=float(len(region)), aperture=area, cut_ports=interior_ports,
                                crossing_seams=[list(e) for e in ends], blank_meter=flight,
                                extrapolation_relative_error=relative_error))
    return records


if __name__=="__main__":
    for r in geometry():
        print(r["board_level"], r["shape"], r["accessible_volume"], r["aperture"],
              r["blank_meter"]["sampled_extrapolated_area"], r["extrapolation_relative_error"])
