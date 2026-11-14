"""Independent integer/rational geometry and charge audit; no producer imports."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations, product
from collections import deque
import hashlib, json

HERE=Path(__file__).resolve().parent
INPUT=json.loads((HERE/'INPUTS.json').read_text())

def read(name):
    h=INPUT['snapshots'][name];b=(HERE/'inputs'/h).read_bytes()
    assert hashlib.sha256(b).hexdigest()==h
    return json.loads(b)

BASE='openwave-tuning/next/CHAIN5/common_charge_topology/'
FIELD='openwave-tuning/next/CHAIN3/native_field_connection/'

def graph():
    # Exact coordinates in Q(sqrt5), with the original source ordering.
    z=(F(0),F(0));o=(F(1),F(0));m=(F(-1),F(0));p=(F(1,2),F(1,2));q=(-p[0],-p[1])
    frame=[(m,p,z),(o,p,z),(m,q,z),(o,q,z),(z,m,p),(z,o,p),(z,m,q),(z,o,q),(p,z,m),(p,z,o),(q,z,m),(q,z,o)]
    internal=[]
    for a,b in combinations(range(12),2):
        ds=[(x[0]-y[0],x[1]-y[1]) for x,y in zip(frame[a],frame[b])]
        norm=(sum(x*x+5*y*y for x,y in ds),sum(2*x*y for x,y in ds))
        if norm==(4,0):internal.append((a,b))
    assert len(internal)==30
    links=[]
    for axis in range(3):
        row=[]
        for i,v in enumerate(frame):
            if v[axis]!=z:continue
            a,b=v[(axis+1)%3],v[(axis+2)%3]
            # Each nonzero component has coefficients of one common sign.
            if (a[0]+a[1])*(b[0]+b[1])>0:
                w=list(v);j=(axis+2)%3;w[j]=(-w[j][0],-w[j][1]);row.append((i,frame.index(tuple(w))))
        assert len(row)==2;links.append(row)
    cells=list(product(range(2),repeat=3));edges=[]
    for c,coord in enumerate(cells):
        edges += [(12*c+a,12*c+b) for a,b in internal]
        for axis,row in enumerate(links):
            end=list(coord);end[axis]^=1;d=cells.index(tuple(end))
            edges += [tuple(sorted((12*c+a,12*d+b))) for a,b in row]
    assert len(set(edges))==len(edges)==288
    return sorted(edges),internal,links,cells

def divergence(edges,values):
    out=[F(0)]*96
    for (a,b),v in zip(edges,values):out[a]+=v;out[b]-=v
    return out

def subtract(a,b):return tuple(x-y for x,y in zip(a,b))
def dot(a,b):return sum(x*y for x,y in zip(a,b))

def segment_intersection(a,b,c,d):
    """Solve in an invertible coordinate plane, then check the third coordinate."""
    u,v,w=subtract(b,a),subtract(d,c),subtract(c,a)
    for i,j in combinations(range(3),2):
        det=u[j]*v[i]-u[i]*v[j]
        if not det:continue
        t=F(w[j]*v[i]-w[i]*v[j],det)
        s=F(u[i]*w[j]-u[j]*w[i],det)
        if not (0<=t<=1 and 0<=s<=1):return 'none'
        if any(t*u[k]-s*v[k]!=w[k] for k in range(3)):return 'none'
        return 'endpoint' if t in (0,1) and s in (0,1) else 'crossing'
    # Parallel lines: reject overlap, permit only a common endpoint.
    axis=next(i for i,x in enumerate(u) if x)
    if any(w[k]*u[axis]!=w[axis]*u[k] for k in range(3)):return 'none'
    lo=max(min(a[axis],b[axis]),min(c[axis],d[axis]))
    hi=min(max(a[axis],b[axis]),max(c[axis],d[axis]))
    return 'none' if lo>hi else 'endpoint' if lo==hi else 'crossing'

def main():
    edges,internal,links,cells=graph();embed=read(BASE+'EMBEDDING.json')
    topology=read(BASE+'TOPOLOGY.json');data=read(FIELD+'COORDINATES.json')
    points=[tuple(x) for x in embed['vertex_integer_coordinates']];den=embed['coordinate_denominator'];period=2*den
    lookup={tuple(v%period for v in p):i for i,p in enumerate(points)};assert len(lookup)==96
    starts=list(map(tuple,embed['oriented_edge_start_lifts']));ends=list(map(tuple,embed['oriented_edge_end_lifts']))
    rebuilt=[(lookup[tuple(x%period for x in a)],lookup[tuple(x%period for x in b)]) for a,b in zip(starts,ends)]
    assert rebuilt==edges
    lows=[tuple(min(a[k],b[k]) for k in range(3)) for a,b in zip(starts,ends)]
    highs=[tuple(max(a[k],b[k]) for k in range(3)) for a,b in zip(starts,ends)]
    globlo=[min(x[k] for x in lows) for k in range(3)];globhi=[max(x[k] for x in highs) for k in range(3)]
    checked=skipped_images=0
    for shift in product(range(-2,3),repeat=3):
        delta=tuple(period*x for x in shift)
        if any(globhi[k]<globlo[k]+delta[k] or globhi[k]+delta[k]<globlo[k] for k in range(3)):
            skipped_images+=1;continue
        for i,(a,b) in enumerate(zip(starts,ends)):
            for j in range(288):
                if shift==(0,0,0) and j<=i:continue
                if any(highs[i][k]<lows[j][k]+delta[k] or highs[j][k]+delta[k]<lows[i][k] for k in range(3)):continue
                c=tuple(x+y for x,y in zip(starts[j],delta));d=tuple(x+y for x,y in zip(ends[j],delta))
                assert segment_intersection(a,b,c,d)!='crossing',(i,j,shift)
                checked+=1
    assert segment_intersection((0,0,0),(2,2,0),(0,2,0),(2,0,0))=='crossing'
    assert segment_intersection((0,0,0),(2,2,0),(0,2,1),(2,0,1))=='none'
    assert segment_intersection((0,0,0),(2,0,0),(1,0,0),(3,0,0))=='crossing'
    # Integral H2: each edge relates adjacent face coefficients by equality.
    sphere=topology['sphere_complex'];faces=sphere['oriented_faces'];se=sphere['edges'];face_adj=[set() for _ in faces]
    boundary=[[0]*len(faces) for _ in se]
    for f,vertices in enumerate(faces):
        for a,b in zip(vertices,vertices[1:]+vertices[:1]):boundary[se.index(sorted((a,b)))][f]+=1 if a<b else -1
    assert boundary==sphere['boundary'];assert sphere['vertices']-len(se)+len(faces)==2
    for row in boundary:
        nz=[(i,x) for i,x in enumerate(row) if x];assert len(nz)==2 and sum(x for i,x in nz)==0
        a,b=nz[0][0],nz[1][0];face_adj[a].add(b);face_adj[b].add(a)
    seen={0};todo=[0]
    while todo:
        for v in face_adj[todo.pop()]:
            if v not in seen:seen.add(v);todo.append(v)
    assert len(seen)==len(faces) # kernel over Z is exactly Z*(1,...,1).
    parts=[list(map(F,x['particular'])) for x in data['charge_sectors']]
    for cell,p in enumerate(parts):
        expect=[F(0)]*96;expect[0]+=1;expect[12*cell]-=1
        assert divergence(edges,p)==expect
        # Independently integrate the cochain to a vertex potential.
        adj=[[] for _ in range(96)]
        for e,(a,b) in enumerate(edges):adj[a].append((b,-p[e]));adj[b].append((a,p[e]))
        potential={0:F(0)};todo=deque([0])
        while todo:
            a=todo.popleft()
            for b,difference in adj[a]:
                proposed=potential[a]+difference
                if b not in potential:potential[b]=proposed;todo.append(b)
                else:assert potential[b]==proposed
        assert len(potential)==96
    # Fresh divergence-free perturbations from every internal triangular face.
    triangles=[t for t in combinations(range(12),3) if all(tuple(x) in internal for x in combinations(t,2))]
    perturb=[F(0)]*288
    for j,(a,b,c) in enumerate(triangles):
        coeff=F((-1)**j,j+19)
        for x,y in ((a,b),(b,c),(c,a)):
            perturb[edges.index(tuple(sorted((x,y))))]+=coeff*(1 if x<y else -1)
    assert not any(divergence(edges,perturb))
    local_checks=0
    for p in parts:
        expected=divergence(edges,p);actual=divergence(edges,[a+b for a,b in zip(p,perturb)])
        assert actual==expected
        for value in actual:assert value.denominator==1;local_checks+=1
    half=[(a+b)/2 for a,b in zip(parts[1],parts[2])];halfdiv=divergence(edges,half)
    assert halfdiv[0]==1 and halfdiv[12]==halfdiv[24]==F(-1,2)
    # The clutching function closes iff its total phase/(2pi) is integral.
    assert halfdiv[12].denominator!=1
    routes=0
    for item in data['currents']:
        path=item['path'];current=[0]*288
        for a,b in zip(path,path[1:]):current[edges.index(tuple(sorted((a,b))))]+=1 if a<b else -1
        assert current==item['integer_current']
        c,d=item['cell'],item['destination'];assert d==c^(1<<(2-item['axis']))
        delta=[b-a for a,b in zip(parts[c],parts[d])]
        assert divergence(edges,current)==divergence(edges,delta)
        assert not any(divergence(edges,[F(j)-v for j,v in zip(current,delta)]))
        routes+=1
    regions=0
    for size in (1,3,7,19,47,95,96):
        vertices={(17*j+5)%96 for j in range(size)};cut=[int(a in vertices)-int(b in vertices) for a,b in edges]
        for c,p in enumerate(parts):assert dot(cut,p)==int(0 in vertices)-int(12*c in vertices);regions+=1
    result={'schema':'oph.common_charge.independent_topology.v1','producer_pins_authenticated_at_snapshot':INPUT['producer_files_verified'],
        'independent_actual_graph_edges':len(edges),'exact_segment_pairs_checked':checked,'larger_periodic_image_offsets':125,'disjoint_image_boxes':skipped_images,
        'primitive_integral_sphere_kernel':'Z generated by all-one face coefficients; connected dual graph, no rank-only inference',
        'local_fiber_flux_checks':local_checks,'all_particulars_are_exact_gradients':True,'native_signed_current_checks':routes,'fresh_regional_flux_checks':regions,
        'half_daughter_clutching_rejected':True,'global_torus_bundle_claimed':False,
        'source_observable_reduction':'Every size2 signed route has the same axis endpoint, so U^dagger P_c U=P_(c xor7) per full macro. Applies to terminal cell-diagonal flux and clock effects with arbitrary field/reference; not intermediate probes or full native failure histories.',
        'full_criterion_decision':'pending separate complete native review and reduced-observable annex audit','official_rows_changed':[]}
    (HERE/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))

if __name__=='__main__':main()
