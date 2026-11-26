"""Material quadrupole and absolute finite-time energy instrument."""
from itertools import combinations

import numpy as np

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.executable_model import carrier

ARMS = 8192
STEP = 1 / 128
NATIVE_LINE = 132.924
STARTS = 2**26


def time_port(energies, reference):
    z = (np.asarray(energies) - reference) * STEP / (2 * np.pi)
    return np.exp(-1j * np.pi * (ARMS - 1) * z) * np.sinc(ARMS * z) / np.sinc(z)


def calibrate(rng):
    echo = law.echo()
    reads = {}
    for load in range(3):
        column = law.SRC.index((load, 0))
        weights = abs(echo[9:, column].reshape(10, 9))**2
        by_difference = np.zeros(5)
        for i, (_, difference) in enumerate(law.LABELS):
            by_difference[difference + 2] += weights[i].sum()
        reads[str(load)] = rng.multinomial(STARTS, by_difference/by_difference.sum()).tolist()
    empty = [0., 1., NATIVE_LINE]
    return dict(difference_order=list(range(-2, 3)), signed_echo=reads,
                vacuum_bright={str(f): int(rng.binomial(STARTS, abs(time_port(0, f))**2))
                               for f in empty})


class Shape:
    def __init__(self, sector):
        graph = sector.car
        vertices, edges = carrier.icosahedron_from_coordinates()
        labels = carrier.isomorphism(graph.seams, edges)
        vertices = vertices[labels] / np.linalg.norm(vertices[0])
        antipodes = np.argmin(vertices @ vertices.T, axis=1)
        self.axes = [j for j, k in enumerate(antipodes) if j < k]
        self.antipodes = antipodes
        self.axis = vertices[self.axes]
        self.sites = np.array(sector.states)
        positions = vertices[self.sites]
        projected = positions @ self.axis.T
        self.cell_q = (3 * projected**2 - 1).sum(axis=1)
        centered = positions - positions.mean(axis=1, keepdims=True)
        radii = (centered**2).sum(axis=2)
        self.centroid_q = (3 * (centered @ self.axis.T)**2 - radii[:, :, None]).sum(axis=1)
        unions = {tuple(sorted(a+b)) for a, b in combinations(graph.faces(), 2)
                  if set(a).isdisjoint(b)}
        self.pairs = np.array([word in unions for word in sector.states])
        self.caps = sector.agreed == 10
        self.sources = {"cap": sector.index[tuple(sorted([0, *graph.neighbours[0]]))],
                        "two_faces": min(np.flatnonzero(self.pairs &
                                           (sector.agreed == max(sector.agreed[self.pairs]))))}

    def read(self, counts):
        hits = int(sum(counts))
        if not hits:
            return {"hits": 0}
        frac = np.array(counts) / hits
        return dict(hits=hits, two_faces=int(sum(counts[self.pairs])),
                    cap=int(sum(counts[self.caps])), cell_q=(frac @ self.cell_q).tolist(),
                    centroid_q=(frac @ self.centroid_q).tolist())


def sample(rng, probabilities):
    """Last category is every rejected time-port event; probabilities are unconditional."""
    p = np.maximum(np.asarray(probabilities, dtype=float).ravel(), 0)
    if p.sum() > 1 + 1e-10:
        raise ValueError("Time-port effects exceed identity")
    weights = np.append(p, max(0., 1 - p.sum()))
    # Use an explicit probability grid so tiny eigensolver differences cannot
    # create platform-dependent multinomial branches at mathematical zeros.
    weights = np.rint(weights * 2**40)
    return rng.multinomial(STARTS, weights / weights.sum())
