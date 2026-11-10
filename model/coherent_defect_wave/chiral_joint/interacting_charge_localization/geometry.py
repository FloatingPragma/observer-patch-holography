"""Actual multigraph, cut charges and all minimal flux supports."""
from collections import Counter,deque
from itertools import combinations
from functools import lru_cache
import sympy as s
from ..periodic_family.model import Configuration,kinetic_terms,cells
from ..peter_weyl_links.clebsch import irreps


@lru_cache(maxsize=1)
def graph():
    terms=kinetic_terms(Configuration())
    assert len(terms)==24 and all(len(t.path)==1 for t in terms)
    edges=tuple((t.path[0],t.left_cell,t.right_cell) for t in terms)
    assert sorted(e[0] for e in edges)==list(range(24))
    distance={0:0};queue=deque([0])
    while queue:
        u=queue.popleft()
        for _,a,b in edges:
            if u not in (a,b):continue
            v=b if a==u else a
            if v not in distance:distance[v]=distance[u]+1;queue.append(v)
    return edges,distance


def component(selected,origin=0):
    visited={origin}
    while True:
        expanded=visited|{v for _,a,b in selected for u,v in ((a,b),(b,a)) if u in visited}
        if expanded==visited:return visited
        visited=expanded


def certificate():
    edges,distance=graph();rows=[]
    for y,ell in sorted(distance.items()):
        if y==0:continue
        count=0
        for chosen in combinations(edges,ell):
            if y not in component(chosen):continue
            degrees=Counter(v for _,a,b in chosen for v in (a,b))
            assert degrees[0]==degrees[y]==1
            assert all(deg==2 for v,deg in degrees.items() if v not in (0,y))
            assert len(degrees)==ell+1 and component(chosen)==set(degrees)
            count+=1
        assert count>0
        rows.append({'cell':y,'coordinate':list(cells(2)[y]),'distance':ell,'minimal_support_count':count})
    # Every link carries equal and opposite abelian endpoint charge, including
    # the complete nonblank carrier. Therefore summing Gauss over a component
    # removes every interior link charge without a state truncation.
    reps=irreps()
    for charge,row in reps.items():
        q=row['currents'][-1]
        assert q==s.Rational(charge,6)*s.eye(row['dimension'])
        assert s.kronecker_product(q,s.eye(row['dimension']))-s.kronecker_product(s.eye(row['dimension']),q.T)==s.zeros(row['dimension']**2)
    return {'actual_stored_edges':[list(edge) for edge in edges],
        'all_nonorigin_cells':rows,'number_of_link_states':sum(row['dimension']**2 for row in reps.values()),
        'complete_link_U1_endpoint_charge_cancellation':True,
        'distance_shell_multiplicities':{str(k):v for k,v in sorted(Counter(distance[y] for y in distance if y).items())},
        'minimal_connected_supports_are_simple_paths':True}
