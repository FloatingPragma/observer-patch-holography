"""Independent bounded polynomial propagation for the fifteen complete terminal laws."""
import json
import math

import numpy as np
from scipy.sparse import eye

from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.executable_model.carrier import federation
from .transport import PROTOCOL, RECEIPT, probability_vector, read_protocol, release, release_metadata


def terminal_labels(sector, home, internal):
    """Decode individual port pairs independently of the propagation replay's masks."""
    local = {tuple(seam): i for i, seam in enumerate(internal)}
    seams = set(sector.car.seams)
    labels = np.empty(sector.P, dtype=np.intp)
    for i, (a, b) in enumerate(zip(sector.aa, sector.bb, strict=True)):
        pair = (a, b)
        if pair in local:
            labels[i] = local[pair]
        elif a // 12 == b // 12 == home:
            labels[i] = len(local)
        elif pair in seams:
            labels[i] = len(local) + (1 if a // 12 != b // 12 else 2)
        else:
            labels[i] = len(local) + 3
    return labels


def checked_steps(h, start, exposures, *, degree=24):
    # Center the Gershgorin interval. All matrix entries and the step are dyadic.
    off = np.asarray(abs(h).sum(axis=1)).ravel() - abs(h.diagonal())
    low, high = min(h.diagonal() - off), max(h.diagonal() + off)
    center = (low + high) / 2
    step = 1 / 8
    matrix = ((h - center * eye(h.shape[0], format="csr")) * step).tocsr()
    radius = float(np.asarray(abs(matrix).sum(axis=1)).max())
    terms = max(np.diff(matrix.indptr))
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
    spec = read_protocol(PROTOCOL)["release"]
    stored = json.loads(RECEIPT.read_text())
    sector = PairSector(federation(spec["level"]))
    home, source = release(sector, spec["cell"])
    if {k: v for k, v in stored.items() if k != "arms"} != release_metadata(sector, home, source):
        raise ValueError("Stored release or terminal addresses differ")
    if set(stored["arms"]) != {str(k) for k in spec["strengths"]}:
        raise ValueError("Stored interactions differ")
    labels = terminal_labels(sector, spec["cell"], home)
    state = np.zeros(sector.P, complex)
    state[sector.index(*source)] = 1
    for strength in spec["strengths"]:
        h = sector.hamiltonian(strength)
        readings = stored["arms"][str(strength)]
        if set(readings) != {str(t) for t in spec["exposures"]}:
            raise ValueError("Stored exposures differ")
        for t, evolved, budget in checked_steps(h, state, spec["exposures"]):
            bins = np.bincount(labels, weights=evolved.real**2 + evolved.imag**2, minlength=34)
            reference = probability_vector(readings[str(t)], 34)
            discrepancy = np.abs(bins - reference).sum()
            bound = discrepancy + budget
            if not bound < 1.2e-10:
                raise ArithmeticError(f"Terminal law exceeds budget at {(strength,t)}: {bound}")
            print(strength, t, f"L1 error <= {bound:.3g}")


if __name__ == "__main__":
    verify()
