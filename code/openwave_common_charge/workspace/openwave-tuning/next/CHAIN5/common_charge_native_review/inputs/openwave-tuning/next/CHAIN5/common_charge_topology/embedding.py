"""Exact rational PL embedding of the actual periodic source graph in T³."""
from bridge import *
from itertools import product
import math

DEN=1<<14

def choose_points(seed):
    graph=native.spatial.graph(2);raw=native.spatial.frame();points=[]
    state=seed;motif=[]
    for row in raw:
        p=[]
        for a,b in row:
            state=(1664525*state+1013904223)%(1<<32)
            p.append(round(DEN*(.5+.05*(float(a)+math.sqrt(5)*float(b))))+(state%15)-7)
        motif.append(p)
    for cell in graph.cell_coordinates:
        for p in motif:points.append([DEN*cell[j]+p[j] for j in range(3)])
    return np.array(points,np.int64)

def edge_lifts(points):
    graph=native.spatial.graph(2);cell_index={tuple(c):j for j,c in enumerate(graph.cell_coordinates)}
    shifts={}
    for cell,c in enumerate(graph.cell_coordinates):
        for axis,row in enumerate(native.spatial.links()):
            other=list(c);other[axis]=(other[axis]+1)%2;destination=cell_index[tuple(other)]
            for a,b in row:
                u,v=12*cell+a,12*destination+b;shift=np.zeros(3,np.int64)
                if c[axis]==1:shift[axis]=2*DEN
                edge=tuple(sorted((u,v)));assert edge not in shifts
                shifts[edge]=shift if u<v else -shift
    starts=[];ends=[]
    for u,v in graph.seams:
        starts.append(points[u]);ends.append(points[v]+shifts.get((u,v),np.zeros(3,np.int64)))
    return np.array(starts),np.array(ends)

def intersection(P,Q,R,S):
    # Python integers keep this certificate independent of machine overflow.
    P,Q,R,S=(tuple(map(int,x)) for x in (P,Q,R,S))
    sub=lambda a,b:tuple(x-y for x,y in zip(a,b))
    cross_product=lambda a,b:(a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
    dot=lambda a,b:sum(x*y for x,y in zip(a,b))
    u=sub(Q,P);v=sub(S,R);w=sub(R,P);cross=cross_product(u,v)
    determinant=dot(cross,w)
    if determinant:return False,'skew'
    norm=dot(cross,cross)
    if norm:
        t=dot(cross_product(w,v),cross);z=dot(cross_product(w,u),cross)
        if 0<=t<=norm and 0<=z<=norm:
            # Intersections are allowed only at the same lifted graph endpoint.
            endpoint=t in (0,norm) and z in (0,norm)
            return not endpoint,'shared_endpoint' if endpoint else 'crossing'
        return False,'coplanar_disjoint'
    if any(cross_product(w,u)):return False,'parallel_disjoint'
    axis=int(np.flatnonzero(u)[0]);left=max(min(P[axis],Q[axis]),min(R[axis],S[axis]));right=min(max(P[axis],Q[axis]),max(R[axis],S[axis]))
    if left>right:return False,'collinear_disjoint'
    if left==right:return False,'shared_endpoint'
    return True,'collinear_overlap'

def certify(points):
    P,Q=edge_lifts(points);lo=np.minimum(P,Q);hi=np.maximum(P,Q);counts={};tested=0
    assert len({tuple(p) for p in points})==96
    # Image offsets beyond one period cannot meet these short edges.
    assert np.max(abs(Q-P))<2*DEN
    assert np.max(np.max(hi,axis=0)-np.min(lo,axis=0))<4*DEN
    for shift in product((-1,0,1),repeat=3):
        displacement=np.array(shift,np.int64)*(2*DEN);R=P+displacement;S=Q+displacement
        low=lo+displacement;high=hi+displacement
        for i in range(len(P)):
            candidates=np.flatnonzero(np.all(hi[i]>=low,axis=1)&np.all(high>=lo[i],axis=1))
            for j in candidates:
                if shift==(0,0,0) and j<=i:continue
                bad,kind=intersection(P[i],Q[i],R[j],S[j]);tested+=1;counts[kind]=counts.get(kind,0)+1
                if bad:return None
    # Every incident tangent is distinct. This gives disjoint tube caps on
    # the boundary of a sufficiently small vertex ball.
    graph=native.spatial.graph(2)
    for vertex in range(96):
        tangents=[]
        for i,(a,b) in enumerate(graph.seams):
            if vertex==a:tangents.append(Q[i]-P[i])
            if vertex==b:tangents.append(P[i]-Q[i])
        assert len(tangents)==6
        for i,u in enumerate(tangents):
            for v in tangents[i+1:]:
                if not np.any(np.cross(u,v)):return None
    maximum=max(1,int(np.max(abs(Q-P))))
    radius=F(1,100*maximum**2*DEN)
    return {'tested_overlapping_bounding_boxes':tested,'exact_pair_classes':counts,
        'rational_vertex_ball_radius':str(radius),'rational_tube_radius':str(radius/(100*maximum**2)),
        'all_actual_edges_embedded':288,'all_vertices':96,'source_cell_size':2,
        'all_incident_tangent_pairs_distinct':True,'native_seams_changed':False,
        'selection_geometry':'Regular-icosahedron port coordinates rounded and perturbed rationally inside each original cubic cell; same motif in every cell. This is a spatial dictionary, not a change to source seam/energy data.',
        'metric_boundary':'No force, speed or energy is inferred from the auxiliary tube radii.'}

def main():
    pins=verify()
    for seed in range(1,65):
        points=choose_points(seed);result=certify(points)
        if result is not None:break
    else:raise AssertionError('no embedded motif selected')
    P,Q=edge_lifts(points)
    result.update({'seed':seed,'coordinate_denominator':DEN,'vertex_integer_coordinates':points.tolist(),
        'oriented_edge_start_lifts':P.tolist(),'oriented_edge_end_lifts':Q.tolist(),
        'topology':'The actual graph embedded in R^3/(2Z)^3, with vertex balls and oriented tube cross-sections. Each degree6 vertex dual surface is a sphere with six labelled tube disks.',
        'image_domain_checked':[-1,0,1]})
    assert verify()==pins;save('EMBEDDING.json',result)
if __name__=='__main__':main()
