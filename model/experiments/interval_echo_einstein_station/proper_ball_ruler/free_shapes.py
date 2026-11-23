"""One free, faithful, equal-energy family at fixed clock-reconstructed area and volume."""
import argparse
from itertools import combinations
import json
from pathlib import Path

import numpy as np

from model.executable_model import federation
from ..autonomous_selection.measure import local_modes, pair_hamiltonian
from ..equilibrium_refinement.measure import apparatus
from ..shape_stationarity.measure import partition

HERE = Path(__file__).resolve().parent
SHOTS = 2**38


def wave_blocks(car,cells,words,h,cut):
    weights=[]
    products=[]
    energy=None
    gates = {c:tuple(np.flatnonzero(car.degree[12*c:12*c+12]==6)) for c in cells}
    modes = {c:local_modes(gates[c]) for c in cells}
    pairs=np.asarray(words)
    for first,last in combinations(sorted(cells),2):
        use=(pairs[:,0]//12==first) & (pairs[:,1]//12==last)
        a,b=modes[first][pairs[use,0]%12],modes[last][pairs[use,1]%12]
        wave=np.zeros((len(words),4),dtype=complex)
        wave[use]=(a[:,:,None]*b[:,None,:]).reshape(-1,4)
        if energy is None:
            energy=float(np.vdot(wave[:,0],h@wave[:,0]).real)
        products.append(wave)
        weights.append(np.sum(abs(wave)**2,axis=1)/4)
    joint=np.column_stack(products)
    error=float(np.max(np.linalg.norm(h@joint-energy*joint,axis=0)))
    np.testing.assert_allclose(joint.conj().T@joint,np.eye(joint.shape[1]),atol=2e-12)
    np.testing.assert_allclose(joint[np.any(np.isin(pairs,list(cut)),axis=1)],0,atol=1e-14)
    if error>1e-11:
        raise ArithmeticError('Stationary cell pairs lack a common native energy')
    return np.array(weights), energy, error


def report(cohorts,scores,eps):
    scores=np.atleast_2d(scores)
    means = np.asarray(cohorts)@scores.T / SHOTS
    variances = [((scores-m[:,None])**2)@(n/SHOTS) for n,m in zip(cohorts,means,strict=True)]
    return dict(signed_slopes=((means[1]-means[0])/(2*eps)).tolist(),
                standard_errors=(np.sqrt(np.sum(variances,axis=0)/SHOTS)/(2*eps)).tolist())


def observe(row,seed):
    car=federation(row['level'])
    cells=row['cells']
    ports,contacts,cut,words=apparatus(car,cells)
    inside,edge,k,_,_=partition(words,cut)
    size=len(words)
    h=pair_hamiltonian(car,contacts,words)
    blocks, energy, error=wave_blocks(car,cells,words,h,cut)
    alpha=np.zeros(len(blocks))
    alpha[:2]=[.75,-.75]
    geometry=np.array(row['per_word_volume_area_rms'])[:,:2]
    references=blocks @ geometry
    features=np.column_stack([np.ones(len(blocks)),references])
    orthogonal,rank,_=np.linalg.svd(features,full_matrices=False)
    active = rank > 1e-12 * rank[0]
    orthogonal = orthogonal[:,active]
    quantum_scores=-np.log1p(alpha)
    coeff=quantum_scores-orthogonal@(orthogonal.T@quantum_scores)
    coeff/=np.maximum(coeff,0).sum()
    tangent=coeff@blocks
    eigenvalues=(1+alpha)/size
    eps=.5*np.min(4*eigenvalues[abs(coeff)>1e-10]/np.abs(coeff[abs(coeff)>1e-10]))
    base=(1+4*(alpha@blocks))/size
    sharp_ref=np.bincount(inside,weights=base)
    sharp_scores=-np.log(sharp_ref)
    analytic_sharp=float(np.bincount(inside,weights=tangent) @ sharp_scores)
    capacity=float(np.sum(tangent*np.log(np.bincount(inside)[inside])))
    analytic_q=float(coeff@quantum_scores)
    geo,bulk,coherent=[],[],[]
    for sign in (-1,1):
        p=base+sign*eps*tangent
        np.testing.assert_allclose(p.sum(),1,atol=2e-12)
        rng_geo,rng_bulk,rng_coh=[np.random.default_rng([621397,seed,sign+1,i]) for i in range(3)]
        geo.append(rng_geo.multinomial(SHOTS,p))
        bulk.append(rng_bulk.multinomial(SHOTS,np.bincount(inside,weights=p)))
        modal=4*eigenvalues+sign*eps*coeff
        coherent.append(rng_coh.multinomial(SHOTS,np.r_[modal,1-modal.sum()]))
    cut_gradient=np.bincount(edge,weights=tangent)
    if np.max(np.abs(cut_gradient))>3e-13:
        raise ArithmeticError('Stationary mixture changed a cut word')
    return dict(level=row['level'],shape=row['shape'],cells=cells,
                common_energy=energy,eigenvector_residual=error,orthogonal_blocks=len(blocks),
                signed_block_weights=coeff.tolist(),signed_displacement=eps,
                independent_constraints=int(active.sum()),
                exact_volume_area_slopes=(geometry.T @ tangent).tolist(),
                exact_sharp_interior_entropy_slope=analytic_sharp,
                exact_quantum_interior_entropy_slope=analytic_q,cut_capacity_slope=capacity,
                max_cut_word_probability_slope=float(np.max(np.abs(cut_gradient))),
                geometry_observations=report(geo,geometry.T,eps),
                sharp_interior_observations=report(bulk,sharp_scores,eps),
                coherent_interior_observations=report(coherent,np.r_[quantum_scores,0],eps),
                geometry_counts=[a.tolist() for a in geo],sharp_counts=[a.tolist() for a in bulk],
                interior_mode_counts=[a.tolist() for a in coherent])


def run():
    archive=json.loads((HERE/'results.json').read_text())
    return dict(seed=621397,shots_per_cohort=SHOTS,
                rows=[observe(row,i) for i,row in enumerate(archive['cases'])])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    path=p.parse_args().output
    result=json.dumps(run(),separators=(',',':'),allow_nan=False)
    with path.open('x') as out:
        out.write(result+'\n')
