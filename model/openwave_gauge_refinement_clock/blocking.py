from fractions import Fraction as F
from itertools import product,combinations
from functools import lru_cache
from math import factorial,comb
import sympy as sy
from model.openwave_su2_spectral.spectrum import exp_bounds


def graph(side):
    vertices=tuple(product(range(side),repeat=3));edges=[];faces=[]
    for v in vertices:
        for a in range(3):
            if v[a]+1<side:
                w=list(v);w[a]+=1;edges.append((v,tuple(w)))
        for a,b in combinations(range(3),2):
            if v[a]+1<side and v[b]+1<side:
                p=list(v);p[a]+=1;q=p.copy();q[b]+=1;r=list(v);r[b]+=1
                faces.append((v,tuple(p),tuple(q),tuple(r)))
    return vertices,tuple(edges),tuple(faces)


def mask(loop,edges):
    index={frozenset(e):i for i,e in enumerate(edges)};out=0
    for v,w in zip(loop,loop[1:]+loop[:1]):out^=1<<index[frozenset((v,w))]
    return out


def subdivide(loop):
    out=[]
    for v,w in zip(loop,loop[1:]+loop[:1]):
        out.append(tuple(2*x for x in v));out.append(tuple(x+y for x,y in zip(v,w)))
    return tuple(out)


def pairwise_haar(loops,edges):
    """Independent exact SU2 fundamental two-matrix integral."""
    index={frozenset(e):i for i,e in enumerate(edges)};occurrences={};n=0
    for loop in loops:
        ids=list(range(n,n+len(loop)));n+=len(loop)
        for j,v in enumerate(loop):
            w=loop[(j+1)%len(loop)];e=index[frozenset((v,w))]
            occurrences.setdefault(e,[]).append((edges[e]==(v,w),ids[j],ids[(j+1)%len(loop)]))
    if any(len(row)%2 for row in occurrences.values()):return F(0)
    assert all(len(row)==2 for row in occurrences.values())
    parent=list(range(n));parity=[0]*n
    def find(i):
        if parent[i]!=i:
            p,d=find(parent[i]);parity[i]^=d;parent[i]=p
        return parent[i],parity[i]
    def join(a,b,d):
        ra,pa=find(a);rb,pb=find(b)
        if ra==rb:return (pa^pb)==d
        parent[rb]=ra;parity[rb]=pa^pb^d;return True
    sign=[]
    for row in occurrences.values():
        (s,a,b),(t,c,d)=row
        if s==t:
            if not join(a,c,1)or not join(b,d,1):return F(0)
            sign.extend((a,b))
        elif not join(a,d,0)or not join(b,c,0):return F(0)
    sign_counts={};constant=0
    for i in sign:
        root,p=find(i);sign_counts[root]=sign_counts.get(root,0)^1;constant^=p
    if any(sign_counts.values()):return F(0)
    components=len({find(i)[0]for i in range(n)})
    return F((-1)**constant*2**components,2**len(occurrences))


def outward(interval,bits=150):
    lo,hi=interval;den=2**bits
    return [str(F((lo*den).__floor__(),den)),str(F((hi*den).__ceil__(),den))]


@lru_cache(None)
def geometry():
    vertices,edges,faces=graph(3);coarse=graph(2);loops=tuple(subdivide(p)for p in coarse[2])
    fm=tuple(mask(p,edges)for p in faces);boundary=tuple(mask(p,edges)for p in loops)
    matches={b:[]for b in boundary};lower=[0,0,0,0];closed_four=0
    for size in range(1,5):
        for chosen in combinations(range(len(faces)),size):
            xor=0
            for p in chosen:xor^=fm[p]
            if size==4 and xor==0:closed_four+=1
            if xor in matches:
                if size<4:lower[size]+=1
                else:matches[xor].append(chosen)
    assert lower==[0,0,0,0]and all(len(x)==1 for x in matches.values())and closed_four==0
    integrals=[]
    for loop,b in zip(loops,boundary):
        tiled=matches[b][0]
        # Four W factors and one normalized character probe; Haar
        # character norm is one. This is independent of convolution.
        value=pairwise_haar((loop,)+tuple(faces[i]for i in tiled),edges)/16
        assert value==F(1,128)
        integrals.append({'coarse_boundary':list(map(list,loop)),
            'unique_four_face_tiling':list(tiled),
            'coarse_character_four_fine_Wilson_Haar_integral_exact':str(value)})
    # A central link sign flips every elementary face and fixes every
    # subdivided coarse link product. Thus the conditional map is even.
    signs=[]
    for v,w in edges:
        axis=next(i for i in range(3)if v[i]!=w[i]);signs.append((-1)**sum(v[:axis]))
    assert all(product_sign(fm[p],signs)==-1 for p in range(len(faces)))
    for v,w in coarse[1]:
        path=(tuple(2*x for x in v),tuple(x+y for x,y in zip(v,w)),tuple(2*x for x in w))
        ids={frozenset(e):i for i,e in enumerate(edges)}
        assert signs[ids[frozenset(path[:2])]]*signs[ids[frozenset(path[1:])]]==1
    return {'coarse_vertices':8,'coarse_links':12,'coarse_faces':6,
        'fine_vertices':len(vertices),'fine_links':len(edges),'fine_faces':len(faces),
        'all_six_tiling_integrals':integrals,'closed_distinct_four_face_sets':closed_four,'all_lower_degree_coarse_character_coefficients_zero':True,
        'central_link_flip_fixes_coarse_links_and_reverses_all_fine_faces':True}


def product_sign(bits,signs):
    result=1
    for i,s in enumerate(signs):
        if bits>>i&1:result*=s
    return result


def coefficient(kappa):
    kappa=F(kappa)
    if not 0<kappa<=F(1,1000):raise ValueError('a declared positive small magnetic step is required')
    scale=36*kappa
    remainder=exp_bounds(scale)[1]*scale**6/factorial(6)
    leading=kappa**4/64
    # Coefficient of W (whose squared Haar norm is 1/4) has error
    # at most 2||R||: equivalently use chi=2W and ||chi||_1<=1.
    interval=(leading-2*remainder,leading+2*remainder)
    return {'fine_magnetic_step_exact':str(kappa),
        'coarse_Wilson_coefficient_interval_exact':outward(interval),
        'leading_coarse_Wilson_coefficient_exact':str(leading),
        'conditional_operator_uniform_remainder_upper_exact':outward((0,remainder))[1],
        'positive_coarse_Wilson_coefficient_certified':interval[0]>0,
        'relative_coefficient_error_upper_exact':outward((0,2*remainder/leading))[1]}


def run():
    geo=geometry();second=F(36,4);fourth=F(36,8)+6*comb(36,2)*F(1,16)
    assert second/2==F(9,2)and fourth/24==F(321,32)
    rows=[coefficient(F(1,10**k))for k in(4,5,6)]
    assert not rows[0]['positive_coarse_Wilson_coefficient_certified']
    assert all(r['positive_coarse_Wilson_coefficient_certified']for r in rows[1:])
    return {'schema':'oph.exact_cubic_Haar_spatial_blocking.v1',
        'verdict':'PASS_NONZERO_COARSE_WILSON_FROM_ACTUAL_SPATIAL_HAAR_BLOCKING',
        'geometry':geo,'exact_fine_magnetic_moments':{'second':str(second),'fourth':str(fourth)},'full_Haar_spaces_and_all_Gauss_transformations_retained':True,
        'embedding':'Each coarse holonomy is the ordered product of its two fine links; all remaining fine Haar variables are integrated.',
        'electric_compression':'J* sum_fine C² J = 2 sum_coarse C²',
        'first_magnetic_compression':'J* sum_36_fine W J = 0',
        'exact_conditional_magnetic_expansion_through_order_four':'E[exp(kappa sum_36 Wfine)|coarse]=1+(9/2)kappa²+(321/32)kappa⁴+(kappa⁴/64)sum_6 Wcoarse+R',
        'uniform_remainder':'||R|| <= exp(36|kappa|)*(36|kappa|)^6/720',
        'declared_magnetic_steps':rows,
        'physical_clock':{'spacing_ratio_fine_over_coarse_exact':'1/2',
            'coarse_electric_coefficient':'x_coarse/(2a_coarse)',
            'compressed_fine_electric_coefficient':'2x_fine/a_coarse',
            'electric_matching_requires':'x_fine=x_coarse/4',
            'matching_electric_term_alone_matches_full_source':False,
            'same_source_Hamiltonian_time_rescaled':False},
        'scope':{'spatial_link_and_face_degrees_actually_added':True,
            'positive_conditional_magnetic_weight':True,'autonomous_full_Hamiltonian_RG_derived':False,
            'higher_coarse_operator_remainder_retained':True,'continuum_noncollapse_proved':False,
            'four_dimensional_Yang_Mills_or_GeV_gap_proved':False}}
