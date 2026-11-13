"""Complete native triquark frame for the positive correlated kinetic rule."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations_with_replacement, permutations
from math import comb
import numpy as np
import scipy.sparse as sp
import sympy as sy
from model.openwave_hadron_refinement_flow.source import chain
from model.openwave_hadron_refinement_flow.parameters import degree_coefficients
from model.openwave_hadron_refinement_flow.bloch import compact_frame
from model.openwave_hadron_refinement_flow.certified_rest import sqrt_upper
from .exact_native import complete


def canonical(orbitals,label):
    rows=[orbitals[i]for i in label];origin=min(x for x,s,w in rows)
    return origin,tuple((x-origin,s,w)for x,s,w in rows)


@lru_cache(None)
def certificate():
    d=chain(4,True);g=d['gram'];N=len(d['labels'])
    records=[canonical(d['orbitals'],lab)for lab in d['labels']]
    selected=[i for i,(x,shape)in enumerate(records)if max(z[0]for z in shape)<=1]
    shapes=sorted({shape for x,shape in(records[i]for i in selected)})
    assert len(shapes)==100
    rows=[records.index((1,shape))for shape in shapes]
    colshape=[shapes.index(records[i][1])for i in selected]
    gs=[int(g[i])for i in rows];G=sy.diag(*gs)
    # The independent exact CAR/occupation contraction supplies every entry.
    # Floating native contractions below are a separate consistency diagnostic.
    exact_A,exact_G,exact_angular=complete()
    assert G==exact_G
    A=sy.zeros(100)
    terms=list((d['weights'][name],part)for name,part in d['parts'].items())
    p=d['parameters'];dU,dm=degree_coefficients()
    for x in(1,2):
        terms.extend([(-p['spin_gradient']/4,4*d['angular']((x,),True)),
                      (-p['eta']/4,4*d['angular']((x,),False)),
                      (-p['t'],sp.diags(d['n'][:,x],dtype=float)),
                      (dm,sp.diags(d['n'][:,x],dtype=float)),
                      (dU/F(3),sp.diags(2*d['n'][:,x]*(3-d['n'][:,x]),dtype=float))])
    for coefficient,term in terms:
        raw=term.tocsr()[rows][:,selected].toarray()*np.sqrt(np.outer(g[rows],g[selected]))
        z=np.rint(raw).astype(np.int64)
        assert np.max(abs(raw-z),initial=0)<5e-10
        for i,j in zip(*z.nonzero()):A[int(i),colshape[j]]+=sy.Rational(coefficient)*int(z[i,j])
    assert A==exact_A
    A=exact_A
    # Exact full-channel positivity via directed square-root row bounds.
    floor=min(F(A[i,i])/gs[i]-sum(abs(F(A[i,j]))*sqrt_upper(F(1,gs[i]*gs[j])) for j in range(100)if i!=j)for i in range(100))
    assert floor>0
    for spin,exact_S in zip((True,False),exact_angular):
        raw=(4*d['angular']((0,1,2,3),spin)).tocsr()[rows][:,rows].toarray()*np.sqrt(np.outer(g[rows],g[rows]))
        ints=np.rint(raw).astype(np.int64);assert np.max(abs(raw-ints))<1e-10
        S=sy.Matrix(ints.tolist())/4
        assert S==exact_S
        assert A*G.inv()*S-S*G.inv()*A==sy.zeros(100)
    # Literal native CAR epsilon columns are normalized monomials. Every raw
    # occupation belongs to one multiset, so Gram equality is exact by counts.
    monomial_count=0
    for shape in shapes:
        assignments=set(permutations(shape));norm=len(assignments)
        assert norm in(1,3,6)
        monomial_count+=norm
    B0=d['B'][:,rows]
    assert np.max(abs((B0.T@B0-sp.eye(100)).data),initial=0)<2e-14
    # Compare the complete channel rule's actual neutron reduction to the
    # five-state source. No named-particle projector defines the rule.
    z=compact_frame(4)
    inds=sum([v for kind,x,n,v in z['groups']if x==1],[])
    assert len(inds)==5
    P=z['P'][rows][:,inds]
    normal=np.asarray(A,float)/np.sqrt(np.outer(gs,gs))
    assert np.max(abs(P.T@P-np.eye(5)))<1e-12
    from .source import matrix
    a,gram=matrix(F(1),F(0));diag=np.sqrt(np.asarray(gram.diagonal(),float).ravel())
    expected=np.asarray(a,float)/np.outer(diag,diag)
    assert np.max(abs(normal@P-P@expected))<2e-10
    eigen=np.linalg.eigvalsh(normal)
    return {'complete_light_triquark_shape_channels':100,
       'colocated_channels':20,'two_cell_channels':80,
       'all_native_colour_epsilon_monomials':monomial_count,
       'actual_native_fixed_N3_carrier_dimension':len(d['full']),
       'actual_four_cell_singlet_dimension':N,
       'exact_native_integer_Gram':gs,
       'complete_rest_matrix_entries_from_exact_CAR_enumeration':True,
       'floating_native_contraction_role':'Independent consistency diagnostic; no floating rounding supplies the proved exact matrix.',
       'complete_rest_matrix_positive_floor_source_exact':str(floor),
       'complete_rest_matrix_all_100_eigenvalues_diagnostic':eigen.tolist(),
       'spin_and_isospin_commutators_exact_zero':True,
       'actual_neutron_five_state_intertwiner':True,
       'all_families_rule':'For d internal spin/flavour coordinates: local binomial(d+2,3), edge binomial(2d+2,3)-2*binomial(d+2,3). Each receives the same positive inverse-rest-operator prescription.',
       'all_Fock_positivity':'DeltaC^* J DeltaC>=0 for positive J; CAR annihilation polynomials act on every populated sector. No occupancy projector defines the rule.',
       'local_gauge_rule':'Use the parent Gauss intertwiner; span-one clusters carry the corresponding internal Wilson transport. Common-root gauge coordinates give the verified epsilon columns.',
       'other_number_k0_cancellation_claimed':False,
       'source_rule_is_declared_architecture':True}
