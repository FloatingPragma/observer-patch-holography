"""Local ring energy under uninterrupted finite or infinite randomized terminal waits."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

from .measure import HERE, SEED, SHOTS, ring_laplacian, search, shapes


def run():
    patch = search.Patch()
    words = np.array(patch.field_component())
    h = ring_laplacian(patch,words).toarray()
    energies, vectors = eigh(h)
    groups = np.r_[0, np.cumsum(np.diff(energies)>1e-10)]
    centers = np.array([energies[groups==g].mean() for g in np.unique(groups)])
    deviations = abs(energies-centers[groups])
    min_gap = min(np.diff(centers))
    if deviations.max()>1e-12 or min_gap<1e-8:
        raise ArithmeticError("Energy clusters poorly separated")
    face = patch.faces.index(tuple(shapes(patch)[0]))
    ids = {word:i for i,word in enumerate(words)}
    pairs = [(i,ids[end]) for i,word in enumerate(words) for f,end in patch.ring_moves(word)
             if f==face and i<ids[end]]
    a,b=np.array(pairs).T
    differences=vectors[a]-vectors[b]
    pointer=differences.T@differences/2
    equal=groups[:,None]==groups[None,:]
    dephased=pointer*equal
    local_limit=np.sum(vectors*(vectors@dephased),axis=1)
    native=h.diagonal()
    cases=[]
    for level in np.unique(native):
        indices=np.flatnonzero(native==level)
        lo=indices[np.argmin(local_limit[indices])]
        hi=indices[np.argmax(local_limit[indices])]
        cases.append((local_limit[hi]-local_limit[lo],lo,hi,level))
    spread,lo,hi,level=max(cases,key=lambda x:x[0])
    # Bound the summed off-block amplitude by ||P||_F using Cauchy-Schwarz
    # against the squared products of normalized spectral amplitudes.
    wait=1e5*np.sqrt(len(pairs))/min_gap
    delta=centers[groups,None]-centers[groups][None,:]
    # All states and read operators are real. Only the cosine term survives.
    average=pointer*np.sinc(delta*wait/np.pi)
    rng=np.random.default_rng([SEED,3])
    records=[]
    for i in (lo,hi):
        held=float(vectors[i]@average@vectors[i])
        if abs(held-local_limit[i])>2.01e-5:
            raise ArithmeticError("Random appointment violates the declared spectral bound")
        records.append(dict(word=int(words[i]),mean_energy=float(level),ground_probability=1/len(words),
                            limit_probability=float(local_limit[i]),finite_probability=held,
                            finite_clicks=int(rng.binomial(SHOTS,held)),
                            second_energy_moment=float(h[i]@h[i]),
                            energy_probabilities=np.bincount(groups,vectors[i]**2).tolist()))
    return dict(protocol_sha256=hashlib.sha256((HERE/'FREE_PROTOCOL.md').read_bytes()).hexdigest(),
                ring_face=patch.faces[face],ring_projector_rank=len(pairs),shots=SHOTS,records=records,
                ground_ring_probability=float(np.sum(abs(differences[:,0])**2)/2),
                mixed_stationary_ring_probability=len(pairs)/len(words),
                max_infinite_time_spread_at_equal_mean_energy=float(spread),
                ground_dimension=int(sum(groups==0)),energy_cluster_sizes=np.bincount(groups).tolist(),
                energy_centers=centers.tolist(),maximum_cluster_radius=float(deviations.max()),
                minimum_intercluster_gap=float(min_gap),terminal_wait_upper=float(wait),
                finite_average_error_bound=2e-5,
                max_relative_eigenpair_residual=float(np.linalg.norm(h@vectors-vectors*energies,axis=0).max()
                                                     /np.linalg.norm(h,ord=np.inf)))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    with parser.parse_args().output.open('x') as file:
        json.dump(run(),file,allow_nan=False,indent=2)
        file.write('\n')
