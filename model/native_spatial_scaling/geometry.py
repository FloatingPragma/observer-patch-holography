"""Source topology and auxiliary affine coordinates, without routing fits."""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from model.federated_observer import source
from model.physical_space_time import uniform_green as affine


def level(value):
    if type(value) is not int or value<0:
        raise ValueError('a nonnegative integer refinement level required')
    return value


def census(value):
    s=4**level(value)
    return {'level':value,'cells':20*s,'ports':240*s,'seams':630*s,
            'inter_seams':30*s,'degree_five_ports':180*s,'degree_six_ports':60*s}


@dataclass(frozen=True)
class Graph:
    level:int
    seams:tuple
    neighbours:tuple
    inter:tuple

    @property
    def ports(self):return len(self.neighbours)

    @property
    def cells(self):return self.ports//12


def audit(graph):
    expected=census(graph.level)
    if graph.ports!=expected['ports'] or len(graph.seams)!=expected['seams']:
        raise ValueError('source census mismatch')
    if graph.seams!=tuple(sorted(set(graph.seams))):
        raise ValueError('canonical unique source seams required')
    adjacency=[[] for _ in range(graph.ports)];within=[set() for _ in range(graph.cells)]
    inter=[];used=set()
    for a,b in graph.seams:
        if not 0<=a<b<graph.ports:raise ValueError('a bounded source edge required')
        adjacency[a].append(b);adjacency[b].append(a)
        if a//12==b//12:within[a//12].add((a%12,b%12))
        else:
            if a in used or b in used:raise ValueError('routing must use distinct port slots')
            used.update((a,b));inter.append((a,b))
    local=set(source.transfer.native_counts().seams)
    if any(edges!=local for edges in within):raise ValueError('original carrier seams changed')
    if tuple(tuple(sorted(row)) for row in adjacency)!=graph.neighbours or tuple(inter)!=graph.inter:
        raise ValueError('source adjacency and inter-seam declarations disagree')
    _,faces=affine.affine_prefix(graph.level)
    dual=tuple(sorted((a//12,b//12) for a,b in inter))
    if dual!=affine.dual_edges(faces):raise ValueError('routing must realize each source dual edge once')
    if sorted(map(len,adjacency))!=[5]*expected['degree_five_ports']+[6]*expected['degree_six_ports']:
        raise ValueError('source degree census differs')
    return expected


@lru_cache(maxsize=4)
def graph(value):
    value=level(value)
    if value>3:raise ValueError('the literal graph control evaluates only levels zero through three')
    raw=source._source().build_federation(value,'port_pair')
    seams=tuple(sorted(tuple(sorted((int(a),int(b)))) for a,b in zip(raw.seam_a,raw.seam_b,strict=True)))
    neighbours=tuple(tuple(sorted(map(int,row))) for row in raw.graph.neighbour_lists)
    result=Graph(value,seams,neighbours,tuple((a,b) for a,b in seams if a//12!=b//12))
    audit(result)
    return result


def halo(g,cell):
    if type(cell) is not int or not 0<=cell<g.cells:raise ValueError('an actual cell required')
    root=set(range(12*cell,12*(cell+1)))
    outer={v for u in root for v in g.neighbours[u]}-root
    if len(outer)!=3 or len({v//12 for v in outer})!=3:
        raise ArithmeticError('three distinct neighbouring cells required')
    window=root|outer
    for v in outer:
        if len(set(g.neighbours[v])&root)!=1 or len(set(g.neighbours[v])-window)!=5:
            raise ArithmeticError('the native three-port halo attachment changed')
        if set(g.neighbours[v])&outer:raise ArithmeticError('halo ports must not be joined')
    return tuple(sorted(root)),tuple(sorted(window))


@lru_cache(maxsize=4)
def coordinates(value):
    """Six integer coefficients; divide by3*2^(level+1) and circumradius.

    An integer pair(a,b) represents a+b sqrt(5). These are comparison
    coordinates, not physical positions or the floating routing vectors.
    """
    g=graph(value);vertices,faces=affine.affine_prefix(value)
    rows=tuple(tuple(sum(vertices[v][i] for v in face) for i in range(6)) for face in faces)
    return tuple(rows[p//12] for p in range(g.ports)),3*2**(value+1)


def add(x,y):return x[0]+y[0],x[1]+y[1]
def multiply(x,y):return x[0]*y[0]+5*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def divide(x,y):
    d=y[0]*y[0]-5*y[1]*y[1]
    if not d:raise ValueError('nonzero quadratic-field divisor required')
    z=multiply(x,(y[0],-y[1]));return F(z[0],d),F(z[1],d)
def sign(x):return affine._qsign(x[0],x[1])
def at_most(x,bound):return sign((x[0]-bound,x[1]))<=0
def squared(row):return affine._squared(row)
def display(x):return [str(F(x[0])),str(F(x[1]))]


def normalized_square(row,denominator):
    # The raw regular icosahedron has circumradius squared(5+sqrt5)/2.
    return divide(squared(row),(F(5*denominator*denominator,2),F(denominator*denominator,2)))


def edge_energy(value):
    g=graph(value);rows,denominator=coordinates(value);energy=(F(),F())
    for a,b in g.inter:
        difference=tuple(x-y for x,y in zip(rows[a],rows[b],strict=True))
        norm=normalized_square(difference,denominator)
        if not at_most(norm,F(4,3*4**value)):
            raise ArithmeticError('the existing affine hinge bound failed')
        energy=add(energy,norm)
    if energy!=energy_formula(value) or not at_most(energy,F(12)):
        raise ArithmeticError('the exact all-level crossing energy changed')
    return energy


def energy_formula(value):
    n=2**level(value)
    return F(20)-F(40,3*n),F(-4)+F(16,3*n)


def geometry_certificate():
    rows=[]
    for ell in range(4):
        g=graph(ell)
        for cell in range(g.cells):halo(g,cell)
        rows.append({**audit(g),'halos_checked':g.cells,'affine_crossing_energy':display(edge_energy(ell)),
                     'native_expanded_geometry_supported':ell<=2})
    return {'finite_source_controls':rows,'formal_level200_census':census(200),
            'all_level_crossing_energy':'20-4sqrt(5)+(-40+16sqrt(5))/(3*2^level)',
            'all_level_argument':'shared midpoint identifiers, simple closed triangulation, distinct routed ports',
            'all_level_floating_constructor_return_claimed':False,
            'new_physical_coordinates_or_spatial_dimension_claimed':False}
