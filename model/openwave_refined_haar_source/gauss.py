"""Exact full-Haar SU(2) spin networks on the 2-by-2-by-2 spatial cube.

Twice-spins and left-associated fusion channels label an orthonormal Gauss basis.
A finite proof frame is only a compression; every Wilson successor is returned.
"""
from fractions import Fraction
from functools import lru_cache
from itertools import product
import sympy as sy
from sympy.physics.wigner import clebsch_gordan


@lru_cache(None)
def graph():
    vertices=tuple(product(range(3),repeat=3)); edges=[]
    for p in vertices:
        for axis in range(3):
            if p[axis]<2:
                q=list(p);q[axis]+=1;q=tuple(q)
                edges.append((p,q)if sum(p)%2==0 else(q,p))
    faces=[]
    for normal in range(3):
        varying=[a for a in range(3)if a!=normal]
        for side in range(3):
            for u0,v0 in product(range(2),repeat=2):
                points=[]
                for u,v in ((0,0),(1,0),(1,1),(0,1)):
                    p=[0]*3;p[normal]=side;p[varying[0]]=u0+u;p[varying[1]]=v0+v;points.append(tuple(p))
                faces.append(tuple(points))
    return vertices,tuple(edges),tuple(faces)


@lru_cache(None)
def incidence():
    vertices,edges,_=graph()
    return tuple(tuple(i for i,e in enumerate(edges)if v in e)for v in vertices)


@lru_cache(None)
def fusion_paths(spins):
    """All left-associated irreducible channels ending at the singlet."""
    if not spins:return ((),)
    partial=[(spins[0],())]
    for k in spins[1:]:
        partial=[(n,path+(n,))for old,path in partial for n in range(abs(old-k),old+k+1,2)]
    return tuple(path for n,path in partial if n==0)


@lru_cache(None)
def invariant(spins,path,magnetic):
    """Normalized CG-chain tensor; magnetic labels are twice-magnetic numbers."""
    if path not in fusion_paths(spins)or len(magnetic)!=len(spins):return sy.S.Zero
    if any(abs(m)>k or (k-m)%2 for k,m in zip(spins,magnetic)):return sy.S.Zero
    if not spins:return sy.S.One
    angular=spins[0];m=magnetic[0];value=sy.S.One
    for k,n,q in zip(spins[1:],path,magnetic[1:]):
        value*=clebsch_gordan(sy.Rational(angular,2),sy.Rational(k,2),sy.Rational(n,2),sy.Rational(m,2),sy.Rational(q,2),sy.Rational(m+q,2))
        if not value:return sy.S.Zero
        angular=n;m+=q
    return sy.simplify(value)if angular==m==0 else sy.S.Zero


def canonical_cycle(loop):
    loop=tuple(loop);choices=[]
    for word in (loop,tuple(reversed(loop))):
        choices.extend(word[i:]+word[:i]for i in range(len(word)))
    return min(choices)


@lru_cache(None)
def short_cycles():
    vertices,edges,_=graph();adj={v:[]for v in vertices}
    for a,b in edges:adj[a].append(b);adj[b].append(a)
    found=set()
    for start in vertices:
        def visit(path):
            for q in adj[path[-1]]:
                if q==start and len(path)in(4,6):found.add(canonical_cycle(path))
                elif len(path)<6 and q>start and q not in path:visit(path+(q,))
        visit((start,))
    return tuple(sorted(found,key=lambda p:(len(p),p)))


@lru_cache(None)
def cycle_edges(loop):
    index={frozenset(e):i for i,e in enumerate(graph()[1])}
    return tuple(index[frozenset((p,loop[(j+1)%len(loop)]))]for j,p in enumerate(loop))


def state_from_cycle(loop):
    spins=[0]*len(graph()[1])
    for i in cycle_edges(loop):spins[i]=1
    paths=[]
    for incident in incidence():
        choices=fusion_paths(tuple(spins[i]for i in incident));assert len(choices)==1
        paths.append(choices[0])
    return tuple(spins),tuple(paths)


@lru_cache(None)
def basis():
    zero=tuple(0 for _ in graph()[1]);paths=tuple(fusion_paths(tuple(0 for _ in row))[0]for row in incidence())
    states=((zero,paths),)+tuple(state_from_cycle(loop)for loop in short_cycles())
    assert len(states)==len(set(states))==201
    assert sum(len(p)==4 for p in short_cycles())==36
    assert sum(len(p)==6 for p in short_cycles())==164
    return states


def electric(state):return Fraction(sum(k*(k+2)for k in state[0]),4)


@lru_cache(None)
def vertex_contraction(left,left_path,right,right_path,positive,negative):
    total=sy.S.Zero
    for old in product(*(range(-k,k+1,2)for k in right)):
        original=invariant(right,right_path,old)
        if not original:continue
        for s in (-1,1):
            new=list(old);new[positive]+=s;new[negative]-=s;new=tuple(new)
            transformed=invariant(left,left_path,new)
            if not transformed:continue
            first=clebsch_gordan(sy.Rational(right[positive],2),sy.Rational(1,2),sy.Rational(left[positive],2),sy.Rational(old[positive],2),sy.Rational(s,2),sy.Rational(new[positive],2))
            second=clebsch_gordan(sy.Rational(right[negative],2),sy.Rational(1,2),sy.Rational(left[negative],2),sy.Rational(old[negative],2),sy.Rational(-s,2),sy.Rational(new[negative],2))
            total+=original*transformed*first*second*s
    return sy.simplify(total)


@lru_cache(None)
def successors(state,loop):
    """Every spin and fusion channel of normalized fundamental Wilson action."""
    vertices,edges,_=graph();vertex_index={v:i for i,v in enumerate(vertices)}
    spins,old_paths=state;chosen=cycle_edges(loop);answers=[]
    for delta in product(*( (-1,1)if spins[e]else(1,)for e in chosen)):
        new=list(spins)
        for e,d in zip(chosen,delta):new[e]+=d
        new=tuple(new);local=[];valid=True
        factor=sy.Rational(1,2)
        for e in chosen:factor*=sy.sqrt(sy.Rational(spins[e]+1,new[e]+1))
        for j,vertex in enumerate(loop):
            vi=vertex_index[vertex];incident=incidence()[vi]
            outgoing=chosen[j];incoming=chosen[(j-1)%len(loop)]
            positive=outgoing if edges[outgoing][0]==vertex else incoming
            negative=incoming if positive==outgoing else outgoing
            left=tuple(new[e]for e in incident);right=tuple(spins[e]for e in incident)
            choices=[]
            for path in fusion_paths(left):
                value=vertex_contraction(left,path,right,old_paths[vi],incident.index(positive),incident.index(negative))
                if value:choices.append((path,value))
            if not choices:valid=False;break
            local.append((vi,tuple(choices)))
        if not valid:continue
        for channels in product(*(choices for _,choices in local)):
            paths=list(old_paths);value=factor
            for (vi,_),(path,coefficient)in zip(local,channels):paths[vi]=path;value*=coefficient
            value=sy.simplify(value)
            if value:answers.append(((new,tuple(paths)),value))
    assert len(answers)==len(set(state for state,_ in answers))
    return tuple(answers)


@lru_cache(None)
def compressed_wilson(loop):
    states=basis();index={s:i for i,s in enumerate(states)};entries={}
    for j,state in enumerate(states):
        for target,value in successors(state,loop):
            if target in index:entries[index[target],j]=value
    matrix=sy.SparseMatrix(len(states),len(states),entries)
    assert matrix==matrix.T
    return matrix


@lru_cache(None)
def loop_phase(loop):
    """Evaluate the normalized spin network at identity; trace(I)=2 oracle."""
    spins,paths=state_from_cycle(loop);active=cycle_edges(loop);total=sy.S.Zero
    for signs in product((-1,1),repeat=len(active)):
        m=[0]*len(spins)
        for i,s in zip(active,signs):m[i]=s
        value=sy.S.One
        for incident,path in zip(incidence(),paths):
            value*=invariant(tuple(spins[i]for i in incident),path,tuple(m[i]for i in incident))
        total+=value
    value=sy.simplify(total*sy.sqrt(2)**len(active)/2)
    assert value in(-1,1)
    return value


def haar_trace_product(loops):
    """Independent exact fundamental Haar oracle for at most three traces.

    Every link appears at most three times. Central parity kills odd powers;
    the remaining two occurrences integrate by delta or epsilon contractions.
    """
    _,edges,_=graph();index={frozenset(e):i for i,e in enumerate(edges)}
    occurrences={};count=0
    for loop in loops:
        ids=list(range(count,count+len(loop)));count+=len(loop)
        for j,p in enumerate(loop):
            q=loop[(j+1)%len(loop)];edge=index[frozenset((p,q))]
            orientation=1 if edges[edge]==(p,q)else -1
            occurrences.setdefault(edge,[]).append((orientation,ids[j],ids[(j+1)%len(loop)]))
    assert all(len(row)<=3 for row in occurrences.values())
    if any(len(row)%2 for row in occurrences.values()):return sy.S.Zero
    parent=list(range(count));parity=[0]*count
    def find(i):
        if parent[i]!=i:
            p,d=find(parent[i]);parity[i]^=d;parent[i]=p
        return parent[i],parity[i]
    def join(a,b,d):
        ra,pa=find(a);rb,pb=find(b)
        if ra==rb:return(pa^pb)==d
        parent[rb]=ra;parity[rb]=pa^pb^d;return True
    signs=[]
    for row in occurrences.values():
        (s,a,b),(t,c,d)=row
        if s!=t:
            if not join(a,d,0)or not join(b,c,0):return sy.S.Zero
        else:
            if not join(a,c,1)or not join(b,d,1):return sy.S.Zero
            signs.extend((a,b))
    roots={};constant=0
    for i in signs:
        root,p=find(i);roots[root]=roots.get(root,0)^1;constant^=p
    if any(roots.values()):return sy.S.Zero
    components=len({find(i)[0]for i in range(count)})
    return sy.Rational((-1)**constant*2**components,2**len(occurrences))


@lru_cache(None)
def source_actions():
    return tuple(tuple(successors(state,face)for state in basis())for face in graph()[2])


@lru_cache(None)
def total_action():
    out=[]
    for j in range(len(basis())):
        vector={}
        for action in source_actions():
            for target,value in action[j]:vector[target]=vector.get(target,0)+value
        out.append(tuple((s,sy.simplify(v))for s,v in vector.items()if v))
    return tuple(out)


@lru_cache(None)
def compressed_sum():
    index={s:i for i,s in enumerate(basis())};entries={}
    for j,column in enumerate(total_action()):
        for target,value in column:
            if target in index:entries[index[target],j]=value
    matrix=sy.SparseMatrix(len(index),len(index),entries)
    assert matrix==matrix.T
    return matrix


def coarse_boundary_loop(face=0):
    """One original coarse face pulled back through two links per edge."""
    from model.openwave_cubic_haar_source.cube import graph as coarse_graph
    points=coarse_graph()[2][face];fine=[]
    for j,p in enumerate(points):
        q=points[(j+1)%4];fine.append(tuple(2*x for x in p));fine.append(tuple(a+b for a,b in zip(p,q)))
    return tuple(fine)
