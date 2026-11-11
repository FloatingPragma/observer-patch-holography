"""Exact attempted-event law on a closed triangular two-face cell.

The cell has oriented seams (0,1),(1,2),(2,0), opposite faces k and -k mod 6,
and three bounded integer loads. Repair and acceptance are the native Law C functions.
"""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import comb

from model.defect_integer_lift.model import COST, KAPPA, accepts, attempt

SEAMS = ((0, 1), (1, 2), (2, 0))
STATES = tuple((*loads, k) for loads in product(range(6), repeat=3) for k in range(6))


def valid_state(state):
    return isinstance(state,tuple) and len(state)==4 and all(type(x) is int for x in state) and state in STATE_SET


def event(state, seam, coin):
    if not valid_state(state) or type(seam) is not int or seam not in range(3) or type(coin) is not int or coin not in (0, 1):
        raise ValueError('state, seam or coin outside the declared finite domain')
    a, b = SEAMS[seam]
    loads, k = list(state[:3]), state[3]
    na, nb, j = attempt(loads[a], loads[b], coin)
    accepted = bool(j and accepts(k, (-k) % 6, j, loads[a] - loads[b]))
    loads[a], loads[b] = na, nb
    nk = (k + j) % 6 if accepted else k
    return tuple(loads) + (nk,), dict(seam=seam, coin=coin, current=j, applied=accepted,
                                     no_operation=(tuple(loads) + (nk,) == state))


STATE_SET = frozenset(STATES)


def step(state, seam, coin):
    return event(state, seam, coin)[0]


@lru_cache(maxsize=1)
def kernel():
    out = {}
    for state in STATES:
        row = defaultdict(F)
        for seam, coin in product(range(3), range(2)):
            row[step(state, seam, coin)] += F(1, 6)
        out[state] = dict(row)
    return out


def initial_distribution():
    return {(*loads, 0): F(1, 216) for loads in product(range(6), repeat=3)}


def evolve(distribution, steps=1):
    if (type(steps) is not int or steps < 0 or any(not valid_state(s) for s in distribution)
            or any(type(v) not in (int,F) or v < 0 for v in distribution.values())
            or sum(distribution.values()) != 1):
        raise ValueError('normalized nonnegative probability and integer attempt count required')
    p = dict(distribution)
    for _ in range(steps):
        q = defaultdict(F)
        for state, mass in p.items():
            for target, weight in kernel()[state].items():
                q[target] += mass * weight
        p = dict(q)
    return p


def energy(state):
    # Twice the accepted register cost change is bounded by KAPPA times the load square decrease.
    return KAPPA * sum(x*x for x in state[:3]) + 2 * (COST[state[3]] + COST[(-state[3]) % 6])


def balanced(state):
    return max(state[:3]) - min(state[:3]) <= 1


def extinct_or_frozen(state):
    return balanced(state) and (state[3] == 0 or len(set(state[:3])) == 1)


def balance_tail_bound(n):
    # At most eight strict load-square descents; while unbalanced their chance is at least 1/3.
    return sum(F(comb(n, k) * 2**(n-k), 3**n) for k in range(min(7, n) + 1))


def observations(p):
    classes = [F(0) for _ in range(6)]
    for state, mass in p.items():
        classes[state[3]] += mass / 2
        classes[(-state[3]) % 6] += mass / 2
    return dict(class_shares=classes, occupied_face_share=1 - classes[0],
                mean_energy=sum(mass*energy(state) for state, mass in p.items()),
                unbalanced_probability=sum(mass for state, mass in p.items() if not balanced(state)),
                unresolved_probability=sum(mass for state, mass in p.items() if not extinct_or_frozen(state)),
                frozen_defect_probability=sum(mass for state, mass in p.items()
                                              if state[3] and len(set(state[:3])) == 1),
                total_probability=sum(p.values()))


def _solve(a, b):
    """Exact elimination for the finite transient blocks of the native kernel."""
    rows=[list(row)+[value] for row,value in zip(a,b)]
    for col in range(len(rows)):
        pivot=next((i for i in range(col,len(rows)) if rows[i][col]),None)
        if pivot is None:raise ValueError('unclassified recurrent block')
        rows[col],rows[pivot]=rows[pivot],rows[col]
        scale=rows[col][col];rows[col]=[v/scale for v in rows[col]]
        for i in range(len(rows)):
            if i==col:continue
            scale=rows[i][col]
            rows[i]=[v-scale*w for v,w in zip(rows[i],rows[col])]
    return [row[-1] for row in rows]


@lru_cache(maxsize=1)
def persistent_probability():
    """Exact eventual occupied-face probability from every state, by descending load-square blocks.

    An equal load-square transition keeps k or lowers its register cost. States in a transient block
    share total load, load square and k. Closed terminal classes have known occupied-face value.
    """
    groups=defaultdict(list)
    for state in STATES:
        groups[(sum(x*x for x in state[:3]),2*COST[state[3]],sum(state[:3]),state[3])].append(state)
    values={}
    for _,group in sorted(groups.items()):
        if all(extinct_or_frozen(s) for s in group):
            values.update({s:F(int(s[3]!=0)) for s in group});continue
        indices={s:i for i,s in enumerate(group)};n=len(group)
        a=[[F(int(i==j)) for j in range(n)] for i in range(n)];b=[F(0) for _ in group]
        for i,state in enumerate(group):
            for target,weight in kernel()[state].items():
                if target in indices:a[i][indices[target]]-=weight
                else:b[i]+=weight*values[target]
        values.update(zip(group,_solve(a,b)))
    return values


@lru_cache(maxsize=None)
def survival_probabilities(steps):
    """Probability of staying outside the closed terminal set through a fixed number of attempts."""
    if type(steps) is not int or steps<0:
        raise ValueError('nonnegative integer attempt count required')
    p={state:F(not extinct_or_frozen(state)) for state in STATES}
    for _ in range(steps):
        p={state:sum(weight*p[target] for target,weight in row.items()) for state,row in kernel().items()}
    return p
