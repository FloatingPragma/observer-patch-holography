"""Exact spatial probe matrices and collision-retaining comparison kernels."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations
from math import comb

from model.native_correlation_sectors import sectors as s
from model.native_local_response import bounds


def probe(family, edges):
    s.fair_graph(family)
    if type(edges) is not tuple or len(edges)!=2:
        raise ValueError('an ordered pair of original disjoint seams required')
    h=s.contrast(*edges,len(family.weights))
    if any(edge not in family.seams for edge in edges):
        raise ValueError('each probe edge must be an actual source seam')
    return h


def menu(family, probes):
    if type(probes) is not tuple or not probes or len(probes)>12:
        raise ValueError('one through twelve finite probe settings required')
    return tuple(probe(family,edges) for edges in probes)


def killed_counts(family, values):
    """2M times the independent pair transfer killed upon coincidence.

    Missing column mass is an observed absorbing collision in the comparison,
    never a conditional surviving probability or a native repair outcome.
    """
    values,_=s.coefficients(values,len(family.weights),2)
    result=s.counts_step(family,values,2)
    edges=set(family.seams)
    for pair,value in values.items():
        if pair in edges:result[pair]=result.get(pair,0)-2*value
    return {pair:value for pair,value in result.items() if value}


def response_matrix(family, probes, attempts, *, killed=False):
    hs=menu(family,probes);s.integer(attempts)
    if attempts>12 or type(killed) is not bool:
        raise ValueError('at most twelve exact sparse attempts and a literal comparison flag required')
    columns=list(hs);denominator=1;scale=2*len(family.seams);rows=[]
    for t in range(attempts+1):
        rows.append(tuple(tuple(F(s.dot(left,right),4*denominator) for right in columns) for left in hs))
        if t<attempts:
            columns=[killed_counts(family,h) if killed else s.counts_step(family,h,2) for h in columns]
            denominator*=scale
    return tuple(rows)


def collision_law(family, initial, attempts):
    """Complete endpoint distribution plus accumulated absorbing collision."""
    s.fair_graph(family);initial,_=s.coefficients(initial,len(family.weights),2);s.integer(attempts)
    if attempts>64 or any(v<0 for v in initial.values()) or sum(initial.values())!=1:
        raise ValueError('a normalized nonnegative exact pair law and at most64 attempts required')
    law={k:F(v) for k,v in initial.items()};failure=F();scale=2*len(family.seams)
    for _ in range(attempts):
        nxt={pair:F(value,scale) for pair,value in killed_counts(family,law).items()}
        if any(v<0 for v in nxt.values()):raise ArithmeticError('killed pair kernel is not substochastic')
        failure+=sum(law.values())-sum(nxt.values());law=nxt
    if sum(law.values())+failure!=1:raise ArithmeticError('absorbing comparison lost complete mass')
    return {'surviving_pairs':law,'collision':failure}


def pair_distance(family,pair):
    s.coefficients({pair:1},len(family.weights),2)
    start,end=pair;adj=[[] for _ in family.weights]
    for a,b in family.seams:adj[a].append(b);adj[b].append(a)
    distance={start:0};queue=[start]
    for a in queue:
        if a==end:return distance[a]
        for b in adj[a]:
            if b not in distance:distance[b]=distance[a]+1;queue.append(b)
    raise ArithmeticError('declared source graph disconnected')


def comparison_error(family, past, future, attempts):
    """Absolute normalized cross-response error; no collision postselection."""
    left=probe(family,future);right=probe(family,past);s.integer(attempts)
    size=sum(abs(v) for v in right.values());initial={p:F(abs(v),size) for p,v in right.items()}
    exact=collision_law(family,initial,attempts)['collision']
    degrees=[0]*len(family.weights)
    for a,b in family.seams:degrees[a]+=1;degrees[b]+=1
    distance=min(pair_distance(family,pair) for pair in right)
    path=(F(comb(attempts,distance))*F(max(degrees),len(family.seams))**distance
          if distance<=attempts else F())
    upper=min(F(1),F(attempts,len(family.seams)),path)
    if exact>upper:raise ArithmeticError('complete collision exceeds the finite path budget')
    factor=F(max(abs(v) for v in left.values())*size,4)
    return {'exact_collision_probability':exact,'collision_upper':upper,
            'normalized_cross_error_upper':factor*exact,'minimum_initial_pair_distance':distance}


def local_error(family, occupied, past, future, attempts, radius):
    probe(family,past);probe(family,future);s.integer(attempts);s.integer(radius)
    variance=4*s.norm_factor(len(family.weights),2,occupied)
    roots=tuple(sorted({x for edge in past+future for x in edge}))
    escape=bounds.radius_bound(family,roots,radius,attempts)
    return {'roots':roots,'window':bounds.ball(family,roots,radius),
            'complete_local_halfdiamond_upper':escape,
            'unnormalized_score_error_upper':2*escape,
            'normalized_response_error_upper':2*escape/variance,
            'normalization_variance':variance}


def cell_probe(family,cell):
    s.integer(cell)
    if 12*(cell+1)>len(family.weights):raise ValueError('an existing twelve-port source cell required')
    edges=tuple(e for e in family.seams if e[0]//12==e[1]//12==cell)
    for first,second in combinations(edges,2):
        if not set(first)&set(second):return first,second
    raise ValueError('two disjoint original internal seams required')


def federation_menu(level):
    family=s.source.federation_family(level,(1,)*12)
    # Keep the two local supports and a separated support in one source.
    first=cell_probe(family,0)
    alternate=next((a,b) for a,b in combinations(tuple(e for e in family.seams if e[1]<12),2)
                   if not set(a)&set(b) and (a,b)!=first)
    inter=next(e for e in family.seams if e[0]<12<=e[1])
    return family,(first,alternate,cell_probe(family,inter[1]//12))
