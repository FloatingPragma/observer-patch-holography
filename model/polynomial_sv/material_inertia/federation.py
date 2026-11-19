"""Complete two-read histories for native single and doubled loads on the official federations."""
import numpy as np
from scipy.sparse.linalg import expm_multiply

from .. import carrier
from ..law import seam_term
from .instrument import counts
from .source import require


def propagate(H, vector, time):
    if time == 0:
        return vector
    return expm_multiply(-1j*time*H, vector, traceA=-1j*time*H.diagonal().sum())


def collect(spec):
    records = {}
    for level in range(3):
        car = carrier.federation(level)
        n = car.n
        mesh = car.adjacency()
        source_port = carrier.EMITTER_SEAM[level][0]
        start = np.zeros(n, complex)
        start[source_port] = 1
        towns = np.arange(n)//12
        regions = np.full(n, len(spec['flight_ports']))
        for label, town in enumerate(spec['flight_ports']):
            regions[towns == town] = label
        for q in (1,2):
            term = seam_term(spec['kappa'])
            # Read the off diagonal and both types of diagonal directly from the qutrit law.
            hop = term[q,3*q]
            local = car.degree*term[q,q] + (6-car.degree)*spec['kappa']*q*q
            from scipy.sparse import diags
            H = (hop*mesh + diags(local)).tocsr()
            samples = []
            for tick in spec['flight_ticks']:
                first = propagate(H, start, tick/spec['tick_divisor'])
                parts = first[:,None]*(regions[:,None] == np.arange(regions.max()+1))
                # Both source sectors evolve, including their offsets and each unnormalized meter branch.
                later = propagate(H, parts, 2)
                joint = np.stack([np.sum(abs(later[regions==r])**2,axis=0)
                                  for r in range(regions.max()+1)],axis=1)
                require(abs(joint.sum()-1)<1e-10,'lost a two-read outcome')
                edges = counts(joint.ravel().cumsum(),spec['shots'])[:,1]
                samples.append(np.diff(edges,prepend=np.int32(0)).reshape(joint.shape))
            records[f'flight_{level}_{q}'] = np.asarray(samples,dtype=np.int32)
    return records
