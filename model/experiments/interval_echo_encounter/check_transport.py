"""Independent bounded polynomial propagation for the fifteen complete terminal laws."""
import json
import math

import numpy as np
from scipy.sparse import eye

from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.executable_model.carrier import federation
from .transport import PROTOCOL, RECEIPT, census, release


def checked_steps(h, start, exposures):
    # Center the Gershgorin interval. All matrix entries and the step are dyadic.
    off = np.asarray(abs(h).sum(axis=1)).ravel() - abs(h.diagonal())
    low, high = min(h.diagonal() - off), max(h.diagonal() + off)
    center = (low + high) / 2
    step = 1 / 8
    matrix = ((h - center * eye(h.shape[0], format="csr")) * step).tocsr()
    radius = float(np.asarray(abs(matrix).sum(axis=1)).max())
    terms = max(np.diff(matrix.indptr))
    degree = 24
    u = np.finfo(float).eps / 2
    def gamma(k):
        return k*u / (1-k*u)
    # Each Taylor term incurs one sparse multiply and scalar division per order.
    # Every addition contributes to at most degree later additions.
    rounding = math.exp(radius) * math.expm1(math.log1p(gamma(degree)) +
                                            radius * gamma(terms + 4))
    tail = math.exp(radius) * radius**(degree+1) / math.factorial(degree+1)
    local_error = rounding + tail
    state = start.copy()
    for k in range(1, 8 * max(exposures) + 1):
        term = state
        total = state.copy()
        for j in range(1, degree + 1):
            term = (matrix @ term) * (-1j / j)
            total += term
        state = total
        if k % 8 == 0 and k//8 in exposures:
            delta = math.expm1(k * math.log1p(local_error))
            # L1 error for a complete orthogonal read and its numerical sums.
            probability_error = delta*(2+delta) + gamma(len(state)+8)*(1+delta)**2
            yield k//8, state, probability_error


def verify():
    spec = json.loads(PROTOCOL.read_text())["release"]
    stored = json.loads(RECEIPT.read_text())
    sector = PairSector(federation(spec["level"]))
    home, source = release(sector, spec["cell"])
    state = np.zeros(sector.P, complex)
    state[sector.index(*source)] = 1
    for strength in spec["strengths"]:
        h = sector.hamiltonian(strength)
        for t, evolved, budget in checked_steps(h, state, spec["exposures"]):
            bins = census(sector, spec["cell"], home, evolved)
            reference = stored["arms"][str(strength)][str(t)]
            if np.shape(reference) != (34,):
                raise ValueError("Incomplete stored census")
            discrepancy = sum(abs(x-y) for x,y in zip(bins, reference))
            bound = discrepancy + budget + len(bins)*0.5e-12
            if not bound < 1e-9:
                raise ArithmeticError(f"Terminal law exceeds budget at {(strength,t)}: {bound}")
            print(strength, t, f"L1 error <= {bound:.3g}")


if __name__ == "__main__":
    verify()
