"""Exact bounded-matter SU(3) fusion improves the Gauss electric lower bound."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
from model.openwave_string_positive_cone.source import local_colour_census


def casimir(p, q):
    return F(p*p + q*q + p*q + 3*p + 3*q, 3)


@lru_cache(None)
def character(p, q):
    # Gelfand--Tsetlin patterns for GL3 highest partition(p+q,q,0).
    out = Counter()
    for a in range(q, p+q+1):
        for b in range(q+1):
            for c in range(b, a+1):
                n1, n2, n3 = c, a+b-c, p+2*q-a-b
                out[n1-n2, n2-n3] += 1
    assert sum(out.values()) == (p+1)*(q+1)*(p+q+2)//2
    return out


@lru_cache(None)
def fusion_certificate():
    local = {tuple(x['highest_weight']) for x in local_colour_census()['irreducible_rows']}
    maximum, gamma = max(casimir(*x) for x in local), F(11, 7)
    assert maximum == F(88, 3)
    cutoff = maximum / gamma
    assert cutoff == F(56, 3) and casimir(7, 0) > cutoff
    small = [x for x in product(range(7), repeat=2) if casimir(*x) < cutoff]
    roots = ((2,-1), (-1,2), (1,1))
    denominator = [(sum(x for x,y in z), sum(y for x,y in z), (-1)**n)
                   for n in range(4) for z in combinations(roots,n)]
    rows, highest = [], F(0)
    for left, right in product(small, repeat=2):
        total = casimir(*left) + casimir(*right)
        if total >= cutoff:
            continue
        tensor = Counter()
        # Opposite edge orientations require the conjugate right irrep.
        for x,u in character(*left).items():
            for y,v in character(right[1],right[0]).items():
                tensor[x[0]+y[0],x[1]+y[1]] += u*v
        admitted = []
        for c in sorted(local):
            mult = sum(sign*tensor[c[0]+x,c[1]+y] for x,y,sign in denominator)
            assert mult >= 0
            if mult:
                assert casimir(*c) <= gamma*total
                admitted.append([*c, mult])
                if total:
                    highest = max(highest, casimir(*c)/total)
        rows.append({'left_irrep':list(left), 'right_irrep':list(right),
                     'edge_Casimir_sum':str(total), 'admitted_local_channels':admitted})
    # The cutoff equality channel is controlled by the analytic large-pair branch.
    assert (8,0) in local
    sharp = Counter()
    for x,u in character(4,0).items():
        for y,v in character(4,0).items():
            sharp[x[0]+y[0],x[1]+y[1]] += u*v
    assert sum(sign*sharp[8+x,y] for x,y,sign in denominator) == 1
    assert casimir(8,0) == gamma*(casimir(4,0)+casimir(0,4))
    return {'complete_local_Fock_dimension':2**24, 'maximum_local_Casimir':str(maximum),
            'all_local_irreps':list(map(list,sorted(local))), 'edge_sum_cutoff':str(cutoff),
            'complete_below_cutoff_edge_pairs':rows,
            'large_edge_pairs_control':'C_local^2<=88/3<=(11/7)(C_left+C_right) for every edge sum>=56/3',
            'opposite_edge_orientation_uses_conjugate_right_representation':True,
            'sharp_local_factor':'11/7', 'uniform_open_chain_factor':'22/7',
            'sharp_channel':{'left':[4,0], 'right':[0,4], 'local':[8,0]},
            'scope':'Every physical local Gauss coupling admitted by the complete one-generation24-mode matter Fock space, with unrestricted SU3 electric-edge irreps and zero boundary flux.'}
