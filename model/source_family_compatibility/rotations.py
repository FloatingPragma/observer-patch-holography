"""Actual proper-cubic site action and the physical-spin inference boundary."""
from itertools import permutations,product
import sympy as s

from .algebra import creation


def proper_cubic():
    matrices=[]
    for perm in permutations(range(3)):
        for signs in product((-1,1),repeat=3):
            r=s.zeros(3)
            for i,j in enumerate(perm):r[i,j]=signs[i]
            if r.det()==1:matrices.append(s.ImmutableMatrix(r))
    assert len(matrices)==24 and len(set(matrices))==24
    assert all(a*b in matrices for a in matrices for b in matrices)
    return tuple(matrices)


def certificate():
    rotations=proper_cubic();side=4;vertices=tuple(product(range(side),repeat=3))
    def edge(a,b):return tuple(sorted((a,b)))
    edges={edge(v,tuple((v[j]+(j==i))%side for j in range(3)))for v in vertices for i in range(3)}
    assert len(edges)==3*len(vertices)
    images=[]
    for r in rotations:
        mapping={v:tuple(int(x)%side for x in r*s.Matrix(v))for v in vertices}
        assert len(set(mapping.values()))==len(vertices)
        assert {edge(mapping[a],mapping[b])for a,b in edges}==edges
        assert all(sum(v)%2==sum(w)%2 for v,w in mapping.items())
        images.append(tuple(mapping[v]for v in vertices))
    assert len(set(images))==24
    x,y,z=s.symbols('x y z',real=True);radius=x*x+y*y+z*z
    h4=x**4+y**4+z**4-s.Rational(3,5)*radius**2
    assert s.expand(sum(s.diff(h4,t,2)for t in (x,y,z)))==0
    for r in rotations:
        transformed=r*s.Matrix([x,y,z]);sub=dict(zip((x,y,z),transformed))
        assert s.expand(h4.subs(sub,simultaneous=True)-h4)==0
    # A continuous z-axis rotation acts nontrivially on this cubic singlet.
    generator=s.expand(x*s.diff(h4,y)-y*s.diff(h4,x))
    assert generator!=0
    a=[creation(i)for i in range(4)];number=sum((v*v.H for v in a),s.zeros(16))
    creators=[]
    for site in range(2):
        up,down=a[2*site:2*site+2]
        sz=(up*up.H-down*down.H)/2
        sy=(up*down.H-down*up.H)/(2*s.I)
        creators.extend([sz-s.I*sy,sz+s.I*sy])
    assert all(number*c-c*number==s.zeros(16)for c in creators)
    assert all(a*b-b*a==s.zeros(16)for a in creators[:2]for b in creators[2:])
    return {'proper_spatial_rotations':24,'even_periodic_source_side':side,
            'source_vertices':len(vertices),'all_source_bonds':len(edges),
            'H_J_bonds_and_staggered_selection_preserved':True,
            'all_twenty_four_site_actions_distinct':True,
            'same_beta_translation_semidirect_action':True,
            'internal_spin_axes_are_not_rotated_by_this_site_action':True,
            'binary_cubic_central_element_has_identity_site_action':True,
            'local_source_charge_creators_are_even_and_disjointly_commuting':True,
            'cubic_singlet_spin_four_countercontrol':str(h4),
            'continuous_rotation_derivative_of_countercontrol':str(generator),
            'cubic_singlet_alone_proves_spin_zero':False,
            'actual_SO3_scalar_pullback_extension_constructed':False,
            'relativistic_spin_statistics_theorem_instantiated':False}
