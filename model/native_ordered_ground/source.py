"""Exact finite staggered-field antiferromagnet, in the Marshall frame."""
from fractions import Fraction as F
from itertools import product
from math import lcm
import numpy as np


def parameters(sites,edges,coupling,field):
    if type(sites) is not int or sites<1 or type(edges) is not tuple:
        raise ValueError('positive site count and literal edge tuple required')
    if type(coupling) is not F or type(field) is not F or coupling<=0 or field<=0:
        raise ValueError('positive rational coupling and selection field required')
    for edge in edges:
        if type(edge) is not tuple or len(edge)!=2 or any(type(i) is not int for i in edge) or not 0<=edge[0]<edge[1]<sites:
            raise ValueError('ordered valid endpoints required; parallel bonds are counted')
    colors={}
    adjacency=[[] for _ in range(sites)]
    for x,y in edges:adjacency[x].append(y);adjacency[y].append(x)
    for start in range(sites):
        if start in colors:continue
        colors[start]=0;queue=[start]
        while queue:
            x=queue.pop()
            for y in adjacency[x]:
                if y in colors:
                    if colors[y]==colors[x]:raise ValueError('bipartite finite graph required')
                else:colors[y]=1-colors[x];queue.append(y)


def cube():
    sites=tuple(product(range(2),repeat=3));index={x:i for i,x in enumerate(sites)}
    edges=[]
    for i,x in enumerate(sites):
        for a in range(3):
            if x[a]==0:
                y=list(x);y[a]=1;edges.append((i,index[tuple(y)]))
    return tuple(edges),tuple(i for i,x in enumerate(sites) if sum(x)%2)


def matrix(sites=8,edges=None,coupling=F(1),field=F(1,4)):
    if edges is None:
        if sites!=8:raise ValueError('explicit edges required outside the demonstration cube')
        edges=cube()[0]
    parameters(sites,edges,coupling,field)
    if sites>8:raise ValueError('explicit host matrix is limited to eight sites; theorem is unrestricted')
    scale=lcm((coupling/2).denominator,(field/2).denominator)
    if scale>10**6:raise ValueError('explicit coefficient denominator exceeds host bound')
    j=int(coupling*scale/2);b=int(field*scale/2);size=1<<sites
    if max(j,b)*max(len(edges),sites)>10**8:raise ValueError('explicit integer matrix exceeds host bound')
    h=np.zeros((size,size),dtype=np.int64)
    for mask in range(size):
        for x,y in edges:
            if ((mask>>x)^(mask>>y))&1:
                h[mask,mask]-=j;h[mask^(1<<x)^(1<<y),mask]-=j
        for x in range(sites):h[mask^(1<<x),mask]-=b
    h.flags.writeable=False
    return h,scale


def perron_bound(sites,edges,coupling,field):
    parameters(sites,edges,coupling,field)
    s=coupling*len(edges)+field*sites
    delta=field/(2*s);minor=delta**sites
    kappa=minor**2/F(2**sites)
    return {'shift_scale':str(s),'positive_N_step_entry_floor':str(minor),
      'Doeblin_minorization':str(kappa),'Hamiltonian_gap_lower':str(s*kappa/sites),
      'block_power_contraction_upper':str(1-kappa),
      'normalized_vector_error_after_m_blocks':'2*sqrt(2^N)*(1-kappa)^m',
      'field_is_strictly_positive':True,'host_matrix_allocation_used':False}


def rational_power(sites,edges,coupling,field,blocks):
    """Exact all-positive power iterates; small-host implementation only."""
    from math import gcd
    parameters(sites,edges,coupling,field)
    if type(blocks) is not int or blocks<0:raise ValueError('nonnegative block count required')
    h,scale=matrix(sites,edges,coupling,field);s=coupling*len(edges)+field*sites
    multiplier=(s*scale).numerator;denominator=(s*scale).denominator
    a=multiplier*np.eye(1<<sites,dtype=object)-denominator*h.astype(object)
    v=np.ones(1<<sites,dtype=object)
    for _ in range(blocks*sites):
        v=a@v;divisor=0
        for x in v:divisor=gcd(divisor,int(x))
        v//=divisor
    return tuple(int(x) for x in v)
