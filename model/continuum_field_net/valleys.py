"""Exact eight-valley phases and signed principal Clifford symbols."""
from itertools import product
import sympy as s


VALLEYS=tuple(product((0,1),repeat=3))


def phase(valley,site):
    if valley not in VALLEYS or len(site)!=3 or any(type(x) is not int for x in site):
        raise ValueError('binary valley and integer three-site required')
    return (-1)**sum(b*n for b,n in zip(valley,site,strict=True))


def certificate():
    pauli=(s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1))
    rows=[]
    for b in VALLEYS:
        signed=tuple((-1)**b[j]*pauli[j] for j in range(3))
        if any(signed[i]*signed[j]+signed[j]*signed[i]!=2*int(i==j)*s.eye(2)
               for i in range(3) for j in range(3)):
            raise ArithmeticError('signed Clifford relation changed')
        for n in product(range(-2,3),repeat=3):
            for j in range(3):
                neighbor=tuple(n[k]+int(k==j) for k in range(3))
                if phase(b,neighbor)!=(-1)**b[j]*phase(b,n):
                    raise ArithmeticError('valley hopping sign changed')
        rows.append({'valley':list(b),'signs':[(-1)**x for x in b],
                     'kinetic_orientation':(-1)**sum(b)})
    return {'valley_count':8,'valleys':rows,
            'signed_Clifford_identities_per_valley':9,
            'site_neighbor_phase_checks':8*125*3,
            'nearest_even_translation_phase_one':all(phase(b,n)==1 for b in VALLEYS
                                                      for n in product((-2,0,2),repeat=3)),
            'continuous_rotation_or_boost_symmetry_inferred':False}
