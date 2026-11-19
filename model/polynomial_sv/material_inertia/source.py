"""Exact sector construction from local seam terms and an optional three-port pair hop."""
from dataclasses import dataclass
from itertools import combinations

import numpy as np
from scipy import sparse

from .. import law


def require(condition, message):
    if not condition:
        raise ValueError(message)


def pair_hop(kappa):
    """Source-internal matching of the pair's curvature to its rest ratio."""
    delta = 1 + 2*kappa
    rest_ratio = (5 + 1/delta)/3
    return (1/rest_ratio - 1/(2*delta))/(2*(1 - delta**-2))


def parameters(member, kappa):
    require(member in ('contrast', 'assisted', 'exchange'), 'unknown local law')
    return (0 if member == 'exchange' else kappa,
            pair_hop(kappa) if member == 'assisted' else 0)


def rotate(mask, size):
    return ((mask << 1) & ((1 << size)-1)) | (mask >> (size-1))


@dataclass
class Ring:
    size: int
    number: int
    states: list
    index: dict
    orbits: list

    @classmethod
    def build(cls, size, number):
        require(1 <= number and 2*number < size, 'unsupported ring sector')
        states = [sum(1 << j for j in sites) for sites in combinations(range(size), number)]
        return cls.from_states(size,number,states,lambda mask: rotate(mask,size))

    @classmethod
    def from_states(cls, size, number, states, shift):
        index = {mask:i for i,mask in enumerate(states)}
        seen = set()
        orbits = []
        for mask in states:
            if mask in seen:
                continue
            orbit = []
            while mask not in seen:
                seen.add(mask)
                orbit.append(index[mask])
                mask = shift(mask)
            orbits.append(orbit)
        require(sum(map(len, orbits)) == len(index) == len(states), 'lost translation orbit')
        return cls(size, number, states, index, orbits)

    def hamiltonian(self, member, kappa):
        k, hop = parameters(member, kappa)
        term = law.seam_term(k)
        rows, cols, data = [], [], []
        for col, mask in enumerate(self.states):
            # Four inactive contrast references per occupied site keep the official degree-six cost.
            diagonal = 4*k*self.number
            for j in range(self.size):
                a, b, c = (1 << ((j+d) % self.size) for d in range(3))
                x, y = int(bool(mask & a)), int(bool(mask & b))
                loc, flipped = 3*x+y, 3*y+x
                diagonal += term[loc,loc]
                if x != y:
                    rows.append(self.index[mask ^ a ^ b])
                    cols.append(col)
                    data.append(term[flipped,loc])
                if hop and y and x != bool(mask & c):
                    diagonal += hop
                    rows.append(self.index[mask ^ a ^ c])
                    cols.append(col)
                    data.append(-hop)
            rows.append(col)
            cols.append(col)
            data.append(diagonal)
        return sparse.csr_matrix((data,(rows,cols)), shape=(len(self.states),)*2)

    def frame(self, ell):
        rows, cols, values = [], [], []
        included = [o for o in self.orbits if (ell*len(o)) % self.size == 0]
        for col, orbit in enumerate(included):
            rows.extend(orbit)
            cols.extend([col]*len(orbit))
            values.extend(np.exp(2j*np.pi*ell*np.arange(len(orbit))/self.size)/np.sqrt(len(orbit)))
        frame = sparse.csc_matrix((values,(rows,cols)),shape=(len(self.states),len(included)))
        # Enumeration starts with the compact cluster; it has a full translation orbit.
        require(included[0][0] == self.index[(1 << self.number)-1], 'wrong preparation orbit')
        return frame


def spectral_measure(H, frame):
    gram_defect = frame.conj().T@frame - sparse.eye(frame.shape[1])
    require(np.max(abs(gram_defect.data), initial=0) < 3e-12, 'sector frame is not isometric')
    image = H@frame
    block = (frame.conj().T@image).toarray()
    defect = image - frame@ sparse.csc_matrix(block)
    require(np.max(abs(defect.data), initial=0) < 3e-12, 'translation reduction leaks')
    # eigh reads one triangle and discards imaginary diagonal entries.
    require(np.max(abs(block-block.conj().T)) < 3e-12, 'sector generator is not Hermitian')
    values, vectors = np.linalg.eigh(block)
    return values, abs(vectors[0])**2
