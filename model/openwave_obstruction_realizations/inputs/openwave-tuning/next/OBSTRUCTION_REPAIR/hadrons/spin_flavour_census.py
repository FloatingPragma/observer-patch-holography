"""Exact characters of 24 CAR orbitals; no spectral or empirical inputs.

Modes are (two cells, three colours, four physical spin/flavour labels).
Character product = product_{colour,alpha} (1+x_c y_alpha)^2.
Multiply by Weyl denominators to extract every SU3-singlet x SU4 irrep.
"""
from collections import defaultdict
from itertools import product,permutations
from math import prod, comb
from fractions import Fraction as F
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent

def sign(p):return (-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))
def shifts(n):
    rho=tuple(range(n-1,-1,-1))
    return [(tuple(rho[i]-rho[p[i]] for i in range(n)),sign(p)) for p in permutations(range(n))]
def dimension(lam):
    n=len(lam)
    return int(prod(F(lam[i]-lam[j]+j-i,j-i) for i in range(n) for j in range(i+1,n)))
def character():
    coeff=defaultdict(int)
    for occupation in product(range(3),repeat=12):
        rows=tuple(sum(occupation[4*i:4*i+4]) for i in range(3))
        cols=tuple(sum(occupation[4*i+j] for i in range(3)) for j in range(4))
        coeff[rows+cols]+=1<<occupation.count(1)
    assert sum(coeff.values())==1<<24
    return coeff

def census():
    weights=character();colourshifts=shifts(3);flavourshifts=shifts(4)
    out=[]
    for k in range(9):
        N=3*k;irreps=[]
        for lam in product(range(7),repeat=4):
            if sum(lam)!=N or list(lam)!=sorted(lam,reverse=True):continue
            multiplicity=0
            for ds,sg in colourshifts:
                cw=tuple(k+d for d in ds)
                for df,sf in flavourshifts:
                    fw=tuple(a+b for a,b in zip(lam,df))
                    multiplicity+=sg*sf*weights.get(cw+fw,0)
            assert multiplicity>=0
            if multiplicity:
                irreps.append(dict(U4_highest_weight=list(lam),SU4_Dynkin_labels=[lam[i]-lam[i+1] for i in range(3)],
                                   SU4_dimension=dimension(lam),multiplicity=multiplicity,
                                   joint_colour_and_flavour_weight_dimension=weights[(k,)*3+lam]))
        expected=dimension((3,)*k+(0,)*(8-k))
        assert sum(r['SU4_dimension']*r['multiplicity'] for r in irreps)==expected
        out.append(dict(N=N,colour_Cartan_zero_dimension=comb(8,k)**3,
                        colour_singlet_dimension=expected,SU4_irreps=irreps,
                        total_highest_weight_matrix_size=sum(r['multiplicity'] for r in irreps),
                        largest_highest_weight_matrix_size=max(r['multiplicity'] for r in irreps)))
    result=dict(schema='oph.obstruction_repair.spin_flavour_census.v1',full_CAR_dimension=1<<24,
        exact_character_total=1<<24,character_weight_count=len(weights),colour='SU3',
        spin_flavour='SU4 with alpha=(up-spin up-flavour,up-spin down-flavour,down-spin up-flavour,down-spin down-flavour)',
        architecture='two spatial cells with spin/flavour independent hopping and electric cut Casimir',
        sectors=out,full_symmetry_reconstruction_exact=True,
        spectrum_computed=False,ground_vacuum_irrep_established=False,
        nucleon_delta_split_forbidden_within_same_SU4_irrep=True,
        symmetric_three_index_SU4_branching={'spin3/2_isospin3/2_dimension':16,'spin1/2_isospin1/2_dimension':4})
    (HERE/'SPIN_FLAVOUR_CENSUS.json').write_text(json.dumps(result,indent=2)+'\n')
    for row in out:print(row,flush=True)
if __name__=='__main__':census()
