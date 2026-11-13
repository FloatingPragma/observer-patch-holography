from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import sympy as sy
from sympy.physics.wigner import wigner_3j,clebsch_gordan


@lru_cache(None)
def graph():
    vertices=tuple(product(range(2),repeat=3));edges=[]
    for p in vertices:
        for axis in range(3):
            if p[axis]==0:
                q=list(p);q[axis]=1;q=tuple(q);edges.append((p,q)if sum(p)%2==0 else(q,p))
    faces=[]
    for normal in range(3):
        varying=[a for a in range(3)if a!=normal]
        for side in range(2):
            points=[]
            for u,v in ((0,0),(1,0),(1,1),(0,1)):
                p=[0]*3;p[normal]=side;p[varying[0]]=u;p[varying[1]]=v;points.append(tuple(p))
            faces.append(tuple(points))
    return vertices,tuple(edges),tuple(faces)


def triangle(a,b,c):
    return abs(a-b)<=c<=a+b and (a+b+c)%2==0


@lru_cache(None)
def spin_census(max_quarter_energy=26):
    vertices,edges,_=graph();incident={v:[i for i,e in enumerate(edges)if v in e]for v in vertices}
    maxk=max(k for k in range(max_quarter_energy+1)if k*(k+2)<=max_quarter_energy)
    assignments=[];labels=[None]*len(edges)
    def recurse(energy):
        remaining=[v for v in vertices if any(labels[i] is None for i in incident[v])]
        if not remaining:
            assignments.append((energy,tuple(labels)));return
        # Finish the most constrained vertex first. Each physical intertwiner
        # is unique because every cube vertex is trivalent.
        vertex=min(remaining,key=lambda v:sum(labels[i]is None for i in incident[v]))
        unknown=[i for i in incident[vertex]if labels[i]is None]
        for choices in product(range(maxk+1),repeat=len(unknown)):
            cost=sum(k*(k+2)for k in choices)
            if energy+cost>max_quarter_energy:continue
            for i,k in zip(unknown,choices):labels[i]=k
            valid=True
            for v in vertices:
                row=[labels[i]for i in incident[v]]
                if all(k is not None for k in row)and not triangle(*row):valid=False;break
            if valid:recurse(energy+cost)
            for i in unknown:labels[i]=None
    recurse(0)
    assert len(assignments)==len(set(x[1]for x in assignments))
    return tuple(sorted(assignments))


def cycles(spins):
    vertices,edges,_=graph();active={i for i,k in enumerate(spins)if k}
    assert all(k in(0,1)for k in spins)
    incident={v:[i for i in active if v in edges[i]]for v in vertices}
    assert all(len(row)in(0,2)for row in incident.values())
    out=[]
    while active:
        edge=min(active);start=edges[edge][0];point=start;path=[]
        while True:
            path.append(point)
            available=[i for i in incident[point]if i in active]
            if not available:break
            take=min(available);active.remove(take)
            point=edges[take][1]if edges[take][0]==point else edges[take][0]
            if point==start:break
        assert point==start and len(path)>=4
        out.append(tuple(path))
    return tuple(out)


@lru_cache(None)
def basis():
    census=spin_census()
    selected=[(energy,spins)for energy,spins in census if energy<=24]
    assert all(all(k in(0,1)for k in spins)for _,spins in selected)
    assert min(energy for energy,_ in census if energy>24)==26
    rows=tuple((F(energy,4),spins,cycles(spins))for energy,spins in selected)
    assert len(rows)==32
    assert sorted(len(loop)for _,_,loops in rows for loop in loops).count(8)==6
    assert sum(len(loops)==2 for _,_,loops in rows)==3
    return rows


def haar_trace_product(loops):
    """Exact pairwise fundamental Haar integral using signed index constraints.

    This oracle handles the required cube matrix elements. Every edge occurs
    at most three times. Odd multiplicity vanishes under the central sign.
    Two equal orientations use epsilon pairings; opposite ones use deltas.
    """
    _,edges,_=graph();edge_index={frozenset(e):i for i,e in enumerate(edges)}
    occurrences={};count=0
    for loop in loops:
        ids=list(range(count,count+len(loop)));count+=len(loop)
        for j,p in enumerate(loop):
            q=loop[(j+1)%len(loop)];e=edge_index[frozenset((p,q))]
            orientation=1 if edges[e]==(p,q)else -1
            occurrences.setdefault(e,[]).append((orientation,ids[j],ids[(j+1)%len(loop)]))
    assert all(len(row)<=3 for row in occurrences.values())
    if any(len(row)%2 for row in occurrences.values()):return sy.S.Zero
    parent=list(range(count));parity=[0]*count
    def find(i):
        if parent[i]!=i:
            p,d=find(parent[i]);parity[i]^=d;parent[i]=p
        return parent[i],parity[i]
    def join(a,b,d):
        ra,pa=find(a);rb,pb=find(b)
        if ra==rb:return (pa^pb)==d
        parent[rb]=ra;parity[rb]=pa^pb^d;return True
    sign_indices=[]
    for row in occurrences.values():
        (s,a,b),(t,c,d)=row
        if s!=t:
            if not join(a,d,0)or not join(b,c,0):return sy.S.Zero
        else:
            if not join(a,c,1)or not join(b,d,1):return sy.S.Zero
            sign_indices.extend((a,b))
    signs={};constant=0
    for i in sign_indices:
        root,p=find(i);signs[root]=signs.get(root,0)^1;constant^=p
    if any(signs.values()):return sy.S.Zero
    components=len({find(i)[0]for i in range(count)})
    return sy.Rational((-1)**constant*2**components,2**len(occurrences))


@lru_cache(None)
def matrices():
    rows=basis();faces=graph()[2];n=len(rows)
    gram=sy.Matrix(n,n,lambda i,j:haar_trace_product(rows[i][2]+rows[j][2]))
    assert gram==sy.eye(n)
    wilson=tuple(sy.Matrix(n,n,lambda i,j:haar_trace_product(rows[i][2]+(face,)+rows[j][2])/2)for face in faces)
    assert all(w==w.T for w in wilson)
    electric=sy.diag(*[sy.Rational(row[0])for row in rows])
    return electric,wilson


def hamiltonian(g2,spacing=F(1)):
    g2,spacing=F(g2),F(spacing)
    if g2<=0 or spacing<=0:raise ValueError('positive declared coupling and spacing required')
    E,faces=matrices();alpha=g2/(2*spacing);beta=2/(g2*spacing)
    return sy.Rational(alpha)*E+sy.Rational(beta)*(6*sy.eye(E.rows)-sum(faces,sy.zeros(E.rows)))


@lru_cache(None)
def face_edges(face):
    _,edges,faces=graph();vertices=faces[face]
    return tuple(i for i,e in enumerate(edges)if e[0]in vertices and e[1]in vertices)


def successors(spins,face):
    selected=face_edges(face)
    for delta in product((-1,1),repeat=4):
        target=list(spins)
        for e,d in zip(selected,delta):target[e]+=d
        if min(target)<0:continue
        target=tuple(target);value=cg_entry(target,spins,face)
        if value:yield target,value


@lru_cache(None)
def vertex_contraction(left,right,positive,negative):
    half=sy.Rational(1,2);total=sy.S.Zero
    magnetic=lambda k:tuple(half*m for m in range(-k,k+1,2))
    for old in product(*(magnetic(k)for k in right)):
        invariant=wigner_3j(*(half*k for k in right),*old)
        if not invariant:continue
        for s in (-half,half):
            new=list(old);new[positive]+=s;new[negative]-=s
            if any(new[i]not in magnetic(left[i])for i in range(3)):continue
            transformed=wigner_3j(*(half*k for k in left),*new)
            first=clebsch_gordan(right[positive]*half,half,left[positive]*half,old[positive],s,new[positive])
            second=clebsch_gordan(right[negative]*half,half,left[negative]*half,old[negative],-s,new[negative])
            total+=invariant*transformed*first*second*(-1 if s==-half else 1)
    return sy.simplify(total)


@lru_cache(None)
def cg_entry(left,right,face):
    vertices,edges,faces=graph();selected=face_edges(face)
    if any(left[e]!=right[e]for e in range(12)if e not in selected):return sy.S.Zero
    if any(abs(left[e]-right[e])!=1 for e in selected):return sy.S.Zero
    value=sy.Rational(1,2)
    for e in selected:value*=sy.sqrt(sy.Rational(right[e]+1,left[e]+1))
    walk=faces[face]
    for vertex in walk:
        incident=sorted(e for e in range(12)if vertex in edges[e]);two=[e for e in selected if vertex in edges[e]]
        j=walk.index(vertex);following=walk[(j+1)%4];previous=walk[(j-1)%4]
        outgoing=next(e for e in two if following in edges[e]);incoming=next(e for e in two if previous in edges[e])
        # Each vertex need not be a source/sink for the coordinate orientation.
        # Convert inverse fundamental matrices to their epsilon-dual form.
        pos=outgoing if edges[outgoing][0]==vertex else incoming
        neg=incoming if pos==outgoing else outgoing
        # This compact formula assumes a bipartite edge orientation.
        assert all((sum(edges[e][0])%2)==0 for e in incident)
        l=tuple(left[e]for e in incident);r=tuple(right[e]for e in incident)
        if not triangle(*l)or not triangle(*r):return sy.S.Zero
        value*=vertex_contraction(l,r,incident.index(pos),incident.index(neg))
    return sy.simplify(value)


@lru_cache(None)
def recoupled_matrices(max_quarter_energy):
    census=spin_census(max_quarter_energy);states=[spins for _,spins in census];index={s:i for i,s in enumerate(states)}
    electric,faces=matrices_on_states(tuple(states))
    return states,electric,faces


@lru_cache(None)
def matrices_on_states(states):
    index={s:i for i,s in enumerate(states)};assert len(index)==len(states)
    faces=[]
    for face in range(6):
        entries={}
        for j,spins in enumerate(states):
            for target,value in successors(spins,face):
                if target in index:entries[index[target],j]=value
        W=sy.SparseMatrix(len(states),len(states),entries);assert W==W.T;faces.append(W)
    electric=sy.diag(*[sy.Rational(sum(k*(k+2)for k in spins),4)for spins in states])
    return electric,tuple(faces)


def trial_state_set():
    states={s for e,s in spin_census(56)}
    for face in range(6):
        spins=[0]*12
        for edge in face_edges(face):spins[edge]=3
        states.add(tuple(spins))
    assert len(states)==785
    return tuple(sorted(states,key=lambda s:(sum(k*(k+2)for k in s),s)))
