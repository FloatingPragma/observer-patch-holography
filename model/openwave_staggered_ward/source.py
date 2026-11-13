"""Exact cochain placement removes the centered-difference checkerboard."""
from fractions import Fraction as F
from itertools import product
import sympy as s
from model.openwave_source_ward import source as parent

add,scale,shift=parent.add,parent.scale,parent.shift
ETA=(-1,1,1,1)

def D(field,axis):
    # Coordinates are integers in half-spacing units. The two samples are
    # a physical spacing apart, so the denominator at a=1 is one.
    return add(shift(field,axis,-1),scale(-1,shift(field,axis,1)))

def parity(index,dimension=4):
    return tuple(sum(i==j for i in index)%2 for j in range(dimension))

def validate_tensor(h):
    assert set(h)==set(product(range(4),repeat=2))
    for (a,b),field in h.items():
        assert field==h[b,a]
        assert all(tuple(x%2 for x in p)==parity((a,b)) for p in field)

def box(field):return add(*(scale(ETA[i],D(D(field,i),i)) for i in range(4)))

def einstein(h):
    validate_tensor(h)
    trace=add(*(scale(ETA[i],h[i,i]) for i in range(4)))
    ricci={}
    for a,b in product(range(4),repeat=2):
        first=add(*(scale(ETA[c],add(D(D(h[c,b],a),c),D(D(h[c,a],b),c))) for c in range(4)))
        ricci[a,b]=scale(F(1,2),add(first,scale(-1,box(h[a,b])),scale(-1,D(D(trace,a),b))))
    scalar=add(*(scale(ETA[i],ricci[i,i]) for i in range(4)))
    G={(a,b):add(ricci[a,b],scale(-F(1,2)*ETA[a]*int(a==b),scalar))for a,b in ricci}
    validate_tensor(G)
    return G

def divergence(G,b):return add(*(scale(ETA[a],D(G[a,b],a)) for a in range(4)))

def gauge(xi):
    assert all(all(tuple(x%2 for x in p)==parity((b,)) for p in xi[b])for b in range(4))
    return {(a,b):add(D(xi[b],a),D(xi[a],b))for a,b in product(range(4),repeat=2)}

def fixture(seed):
    h={}
    for a in range(4):
        for b in range(a,4):
            base=parity((a,b));field={}
            for j in range(3):
                p=tuple(base[k]+2*((j*(k+1)+seed+a-b)%3-1)for k in range(4))
                field[p]=field.get(p,F(0))+F((a+1)*(b+2)*(j+1)-seed,seed+j+3)
            h[a,b]=h[b,a]={p:v for p,v in field.items()if v}
    return h

def pairing(h,k):
    return sum(ETA[a]*ETA[b]*sum(v*k[a,b].get(p,F(0)) for p,v in h[a,b].items())for a,b in h)

def periodic_gradient(side,collocated=False):
    # Scalar vertices x=0,...,n-1. Edge-centred derivative is f(x+1)-f(x).
    return s.Matrix(side,side,lambda i,j:(int(j==(i+1)%side)-int(j==(i-1)%side))/s.Integer(2)
                    if collocated else int(j==(i+1)%side)-int(j==i))

def static(phi,gamma):
    lap=add(*(D(D(phi,i),i) for i in range(3)))
    G={(0,0):scale(2*gamma,lap)}
    for i,j in product(range(3),repeat=2):
        G[i+1,j+1]=scale(gamma-1,add(D(D(phi,i),j),scale(-int(i==j),lap)))
    for i in range(1,4):G[0,i]=G[i,0]={}
    return G

def run():
    assert all(parent.symbolic_identities().values())
    h,k=fixture(1),fixture(2);G=einstein(h)
    assert all(not divergence(G,b) for b in range(4))
    xi={b:{parity((b,)):F(b+1,7)}for b in range(4)}
    pure=gauge(xi);assert all(not x for x in einstein(pure).values())
    assert pairing(h,einstein(k))==pairing(k,G)
    assert pairing(pure,G)==0
    scalar,center,stiffness=parent.profile()
    fine={tuple(2*x for x in p):v for p,v in scalar.items()}
    lap=add(*(D(D(fine,i),i)for i in range(3)))
    expected={tuple(2*x for x in p):v for p,v in parent.nearest_laplacian(scalar).items()}
    assert lap==expected
    phi=scale(-1,fine);actual=static(phi,F(1))
    assert all(not actual[i,j] for i,j in product(range(1,4),repeat=2))
    # On a periodic connected scalar graph, the kernel is exactly constants;
    # the old collocated central stencil has an extra checkerboard at even n.
    controls=[]
    for n in (4,5,6,8):
        B=periodic_gradient(n);old=periodic_gradient(n,True)
        assert B.rank()==n-1
        checker=s.Matrix([(-1)**j for j in range(n)])
        if n%2==0:assert old*checker==s.zeros(n,1) and B*checker!=s.zeros(n,1)
        L=B.T*B
        assert L==s.Matrix(n,n,lambda i,j:2*int(i==j)-int(j==(i+1)%n)-int(j==(i-1)%n))
        controls.append({'periodic_vertices':n,'staggered_gradient_rank':B.rank(),
            'collocated_central_gradient_rank':old.rank(),'constant_kernel_only':True,
            'complete_scalar_Laplacian_matches_nearest_neighbour':True})
    return {'schema':'oph.staggered_linearized_curvature.v1','verdict':'PASS_EXACT_STAGGERED_WARD_AND_CONNECTED_SCALAR_KERNEL',
        'construction':{'spacing':'1','stored_coordinate_unit':'one half of the declared record spacing',
            'scalar_and_diagonal_tensor_locations':'even coordinate lattice',
            'vector_component_locations':'parity e_b','tensor_component_locations':'parity e_a+e_b',
            'derivative':'D_a f(y)=f(y+e_a)-f(y-e_a); target parity changes by e_a',
            'full_zero_extension_and_all_halos_retained':True},
        'exact_identities':{'all_four_Bianchi_divergences_zero':True,'pure_gauge_Einstein_tensor_zero':True,
            'linearized_action_self_adjoint':True,'gauge_direction_action_variation_zero':True,
            'native_profile_Laplacian_equals_original_nearest_neighbour':True},
        'native_profile':{'sites':len(scalar),'stiffness_exact':str(stiffness),
            'center':list(center),'all_halo_laplacian_max':str(parent.maximum(lap)),
            'gamma_one_static_spatial_stress_exact_zero':True},
        'complete_periodic_controls':controls,
        'error_transport':{'edge_bound':'|lambda(x+e_j)-lambda(x)| <= a*w + 8*delta',
            'premises':'same staggered tensor components; every remainder entry bounded by delta; staggered divergence bounded by w',
            'path_bound':'sum edge budgets on the actual connected nearest-neighbour scalar graph',
            'nine_null_decoder_transfer':'requires an explicit co-location or transport map and its error; not inferred for generally staggered components'},
        'scope':{'exact_finite_differential_identity':True,'checkerboard_ambiguity_removed_on_declared_scalar_graph':True,
            'physical_metric_or_stress_derived':False,'nonlinear_Einstein_equation':False,
            'co_located_nine_null_observer_instrument_derived':False,'four_dimensional_continuum':False,
            'new_empirical_fit':False}}
