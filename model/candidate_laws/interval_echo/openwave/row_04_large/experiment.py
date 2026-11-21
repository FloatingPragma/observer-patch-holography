"""Certified radius-six two-load release on the registered 320-cell federation."""
import argparse
from collections import deque
import json
import math
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply
from scipy.stats import beta

from model.executable_model import federation
from model.executable_model.carrier import Carrier
from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.row_04.__main__ import evolve, verify_local_reader


HERE=Path(__file__).resolve().parent
HOME_CELL=13
RADIUS=6
SHOTS=16_384
NUMERIC=1e-9


def ball(car,radius):
    adjacent=[[] for _ in range(car.cells)]
    external_per_port=np.zeros(car.n,int)
    for a,b in car.external:
        adjacent[a//12].append(b//12)
        adjacent[b//12].append(a//12)
        external_per_port[[a,b]]+=1
    assert int(external_per_port.max())==1
    dist={HOME_CELL:0}
    queue=deque([HOME_CELL])
    while queue:
        c=queue.popleft()
        for n in adjacent[c]:
            if n not in dist:
                dist[n]=dist[c]+1
                queue.append(n)
    return sorted(c for c,v in dist.items() if v<=radius)


def subspace(car,cells):
    selected=np.array([12*c+j for c in cells for j in range(12)])
    lookup=np.full(car.n,-1,int)
    lookup[selected]=np.arange(len(selected))
    edges=[(int(lookup[a]),int(lookup[b])) for a,b in car.seams if lookup[a]>=0 and lookup[b]>=0]
    patch=Carrier('numerical projection',len(selected),edges)
    sector=PairSector(patch)
    missing=car.degree[selected]-patch.degree
    # One half per omitted outgoing hop for each occupied port. Adjacency between
    # the two occupied sites is identical in the full and projected graphs.
    boundary=(missing[sector.aa]+missing[sector.bb])/2
    exchange=sector.exchange()+sparse.diags(boundary,format='csr')
    return lookup,sector,exchange


def observe(level):
    car=federation(level)
    internal=[(a,b) for a,b in car.seams if a//12==b//12==HOME_CELL]
    assert len(internal)==30
    origin=next(e for e in internal if max(car.degree[list(e)])==6)
    cells=ball(car,RADIUS)
    lookup,sector,exchange=subspace(car,cells)
    addresses=np.array([sector.index(*lookup[list(e)]) for e in internal])
    initial=np.zeros(sector.P)
    initial[sector.index(*lookup[list(origin)])]=1
    tail=2*math.e/math.factorial(RADIUS+1) if len(cells)<car.cells else 0.
    error=2*tail+NUMERIC
    full=PairSector(car) if level==0 else None
    arms=[]
    for i,k in enumerate((candidate.KAPPA,0.)):
        H=exchange-sparse.diags(2*k*sector.adjacent.astype(float),format='csr')
        state=evolve(H,k,initial,[1.])[0]
        normalization=float(abs(np.vdot(state,state)-1))
        assert normalization<NUMERIC
        weights=abs(state[addresses])**2
        total=float(weights.sum())
        p=np.append(weights,1-total)
        counts=np.random.default_rng([202611214,level,i]).multinomial(SHOTS,p)
        hits=SHOTS-int(counts[-1])
        report=dict(kappa=k,seam_or_quiet_probabilities=p.tolist(),seam_or_quiet_counts=counts.tolist(),
                    local_probability=total,full_probability_lower_bound=max(0,total-error),
                    full_probability_upper_bound=min(1,total+error),norm_roundoff=normalization,
                    simulated_local_hits=hits,
                    sampling_95=[float(beta.ppf(.025,hits,SHOTS-hits+1)),
                                 float(beta.ppf(.975,hits+1,SHOTS-hits))])
        if full is not None:
            full_seed=np.zeros(full.P)
            full_seed[full.index(*origin)]=1
            exact=expm_multiply(-1j*full.hamiltonian(k),full_seed)
            full_address=[full.index(*e) for e in internal]
            direct=abs(exact[full_address])**2
            maxdiff=float(np.max(abs(weights-direct)))
            assert maxdiff < error
            report['full_exponential_max_address_error']=maxdiff
            report['full_exponential_local_probability']=float(direct.sum())
        arms.append(report)
    return dict(level=level,cells=car.cells,full_ports=car.n,full_physical_seams=len(car.seams),
                full_pair_configurations=car.n*(car.n-1)//2,source_cell=HOME_CELL,
                source_seam=list(origin),ordered_seams=[list(e) for e in internal],quiet_index=30,
                exposure=1.,radius=RADIUS,retained_cells=cells,retained_ports=sector.n,
                retained_pair_configurations=sector.P,state_cutoff_bound=tail,
                probability_error_bound=error,numerical_allowance=NUMERIC,shots=SHOTS,arms=arms,
                native_gate=arms[0]['full_probability_lower_bound']>.80,
                control_gate=arms[1]['full_probability_upper_bound']<.65)


def validate_projection():
    """Independently slice every matrix entry on a strict two-shell subset of the 20-cell sphere."""
    car=federation(0)
    full=PairSector(car)
    cells=ball(car,2)
    mapping,part,exchange=subspace(car,cells)
    selected=np.flatnonzero(mapping>=0)
    indices=np.asarray([full.index(int(a),int(b)) for a,b in zip(selected[part.aa],selected[part.bb])])
    arms=[]
    for k in (candidate.KAPPA,0.):
        compressed=exchange-sparse.diags(2*k*part.adjacent.astype(float),format='csr')
        difference=(compressed-full.hamiltonian(k)[indices][:,indices]).tocsr()
        difference.eliminate_zeros()
        assert difference.nnz==0
        arms.append(dict(kappa=k,disagreeing_entries=int(difference.nnz)))
    return dict(radius=2,retained_cells=len(cells),retained_configurations=part.P,arms=arms)


def portable(x):
    if isinstance(x,dict): return {k:portable(v) for k,v in x.items()}
    if isinstance(x,list): return [portable(v) for v in x]
    if isinstance(x,float): return 0.0 if abs(x)<1e-12 else round(x,12)
    return x


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    verify_local_reader()
    result=portable(dict(full_matrix_compression_check=validate_projection(),screens=[observe(level) for level in (0,2)]))
    for r in result['screens']:
        print(r['cells'],r['source_seam'],len(r['retained_cells']),r['retained_pair_configurations'],
              [(a['kappa'],a['local_probability'],a['simulated_local_hits']) for a in r['arms']],
              r['probability_error_bound'],r['native_gate'],r['control_gate'])
    output=HERE/'results.json'
    if args.write:
        with output.open('x') as f:
            json.dump(result,f,indent=2,allow_nan=False)
            f.write('\n')
    else:
        assert result==json.loads(output.read_text())
