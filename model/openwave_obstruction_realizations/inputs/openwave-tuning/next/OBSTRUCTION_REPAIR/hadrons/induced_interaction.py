"""Derive low-energy interactions from actual staggered colour hopping.

The effective degree of freedom is baryon occupancy, not a physical spin.
"""
from fractions import Fraction as F
from itertools import combinations,product
from pathlib import Path
import json
import numpy as np
from staggered_vacuum import creation,annihilation
from staggered_chain4 import colour6
HERE=Path(__file__).resolve().parent


def main():
    records=[]
    for L,k in [(2,1),(4,1),(4,2),(4,3)]:
        choices=list(combinations(range(L),k))
        words=sorted(sum(1<<(3*x+c) for c,sites in enumerate(sets) for x in sites)
                     for sets in product(choices,repeat=3))
        ix={w:i for i,w in enumerate(words)}
        n=len(words);K=np.zeros((n,n),dtype=np.int64)
        for j,w in enumerate(words):
            for x in range(L-1):
                for c in range(3):
                    for a,b in [(3*x+c,3*(x+1)+c),(3*(x+1)+c,3*x+c)]:
                        z,s=annihilation(b,w)
                        if z is None:continue
                        y,t=creation(a,z)
                        if y is not None:K[ix[y],j]-=s*t
        ids=[i for i,w in enumerate(words) if all(((w>>(3*x))&7) in [0,7] for x in range(L))]
        Z=np.eye(n,dtype=np.int64)[:,ids]
        E6=sum((colour6(words,list(range(c))) for c in range(1,L)),np.zeros((n,n),dtype=np.int64))
        C=colour6(words,list(range(L)))
        W=K@Z
        assert not np.any(E6@Z) and not np.any(C@Z) and not np.any(C@W)
        assert np.array_equal(E6@W,8*W)
        assert not np.any(Z.T@W)
        second=W.T@W
        third=W.T@K@W
        domainwalls=[];hops=np.zeros_like(second)
        pwords=[words[i] for i in ids];pix={w:i for i,w in enumerate(pwords)}
        for j,w in enumerate(pwords):
            occ=[int((w>>(3*x))&7==7) for x in range(L)]
            domainwalls.append(sum(occ[x]!=occ[x+1] for x in range(L-1)))
            for x in range(L-1):
                if occ[x]!=occ[x+1]:
                    y=w^(7<<(3*x))^(7<<(3*(x+1)))
                    hops[pix[y],j]=1
        assert np.array_equal(second,np.diag([3*d for d in domainwalls]))
        assert np.array_equal(third,-6*hops)
        records.append(dict(L=L,baryon_count=k,weight_dimension=n,low_manifold_dimension=len(ids),
            exact_one_hop_colour_energy='4 kappa/3',domainwalls=domainwalls,
            exact_PV2P_without_t=second.tolist(),exact_PV3P_without_t=third.tolist(),
            exact_global_colour_support=True))
    rows=[]
    t=F(1,8)
    for kappa in [F(1,2),F(1),F(2)]:
        rows.append(dict(kappa=str(kappa),domainwall_energy_order2=str(-9*t*t/(4*kappa)),
            baryon_neighbour_hopping_order3=str(-27*t**3/(8*kappa*kappa))))
    out=dict(schema='oph.obstruction_repair.staggered_induced_interaction.v1',
        exact_matrix_checks=records,derived_coefficients=rows,
        H_effective_order2='-(9 t^2/(4 kappa)) sum_edges 1_(occupied/empty or empty/occupied)',
        H_effective_order3='-(27 t^3/(8 kappa^2)) sum_edges (B_left† B_right + h.c.)',
        approximation_domain='perturbation expansion about zero staggered mass and t=0 in the same exact source. Finite nonzero-t exact spectra are in the vacuum result files. No truncated-series equality to full dynamics is asserted.',
        physical_spin_identification=False,coefficient_fitted_to_data=False,
        physical_nucleon_mass_solved=False,
        conclusion='The new source generates its baryon-density interaction and hopping from literal colour-current dynamics; they are not installed as independent phenomenological coefficients. This does not establish physical spin, a continuum, or the observed mass scale.')
    (HERE/'INDUCED_INTERACTION_RESULT.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'checks':len(records),'coefficients':rows},indent=2))


if __name__=='__main__':main()
