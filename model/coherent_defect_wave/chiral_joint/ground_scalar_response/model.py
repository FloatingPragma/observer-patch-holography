"""Actual scalar star identity, neutral currents and complete source matrices."""
from functools import lru_cache
import sympy as s

from ..ground_composite_response.model import source,polynomial_map,validate_annihilator
from ..higgs_composites import model as car
from ..local_reference import algebra,characters
from ...chiral_higgs.model import higgs_memory


def validate_scalar_star(h,n):
    y=h-n
    if n*n!=n:raise ValueError('actual scalar excitation operator must be a projector')
    if n*y+y*n!=y:raise ValueError('actual scalar Yukawa anticommutator fails')
    squared=h*h-h
    if y*y!=squared or n*squared!=squared*n:
        raise ValueError('scalar squared-transfer conservation fails')
    return y,squared


@lru_cache(maxsize=1)
def matrices():
    scalar=higgs_memory()[2];operator=car.monomial(matrix=scalar)
    blocks=[]
    for two_sz in (-1,1):
        block=source(two_sz)
        raw=polynomial_map(block['keys'],block['keys'],operator)
        image=raw*block['basis']
        n=block['left']*image.extract(block['coordinates'],range(block['basis'].cols))
        assert block['basis']*n==image
        blocks.append((block,n))
    h=s.ImmutableSparseMatrix(s.diag(*(block['H'] for block,n in blocks)))
    g=s.ImmutableSparseMatrix(s.diag(*(block['gram'] for block,n in blocks)))
    n=s.ImmutableSparseMatrix(s.diag(*(n for block,n in blocks)))
    y,squared=validate_scalar_star(h,n)
    assert all(g*op==op.H*g for op in (h,n,y,squared))
    current=s.ImmutableSparseMatrix(s.I*(h*n-n*h))
    assert g*current==current.H*g and h*y-y*h==-(h*n-n*h)
    assert current*current==squared
    x=s.Symbol('x');p=s.Poly(h.charpoly(x).as_expr(),x).sqf_part()
    assert p.degree()==9
    validate_annihilator(h,p)
    return {'H':h,'G':g,'N':n,'Y':s.ImmutableSparseMatrix(y),
            'S':s.ImmutableSparseMatrix(squared),'current':current,'minimal_polynomial':p}


@lru_cache(maxsize=1)
def source_identities():
    n=car.monomial(matrix=higgs_memory()[2])
    parts,couplings=car.quantum_yukawa()
    y=car.add(*parts.values())
    y=car.clean({key:value.subs({v:1 for v in couplings.values()}) for key,value in y.items()})
    h=car.add(n,y)
    assert car.multiply(n,n)==n and n==car.adjoint(n) and y==car.adjoint(y)
    assert car.anticommutator(n,y)==y
    assert car.multiply(y,y)==car.add(car.multiply(h,h),car.scale(h,-1))
    current=car.scale(car.commutator(h,n),s.I)
    assert car.multiply(current,current)==car.multiply(y,y)
    for index in range(12):
        assert all(not car.commutator(car.gauge_current(index),op) for op in (n,y))
    for axis in range(3):
        assert all(not car.commutator(car.spin_current(axis),op) for op in (n,y))
    raw_l=s.diag(*(-1 if i%15==9 else 1 if i%15 in (10,11) else 0 for i in range(30)))
    raw_b=s.diag(*(1 if i%15 in (2,3,5,6,7,8) else -1 if i%15 in (0,1,4,12,13,14) else 0 for i in range(30)))
    assert all(not car.commutator(car.one_body(charge),op) for charge in (raw_l,raw_b) for op in (n,y))
    assert all(len(word)%2==0 for op in (n,y) for word,row,col in op)
    # This bound covers every color charge and every quark occupation.
    a=algebra.certificate();c=characters.certificate()
    assert c['quark_Casimir_operator_upper_bound']=='54'
    assert a['lepton_squared_transfer_operator_upper_bound']==6
    assert s.Rational(165-114,36)**2>2
    return {'full_native_operator_N_squared_equals_N':True,
        'full_native_anticommutator_N_Y_equals_Y':True,
        'full_native_Y_squared_equals_H_squared_minus_H':True,
        'full_native_scalar_current_squared_equals_Y_squared':True,
        'gauge_commutator_checks':24,'spin_commutator_checks':6,
        'raw_L_and_B3_conserved':True,'operator_fermion_parity':'even','operator_physical_spin':'0',
        'Y_definition':'sum of all three actual native Yukawa channels, each coefficient one',
        'individual_Yukawa_channels_substituted':False,
        'full_unconstrained_Y_squared_upper_bound':'114+36 sqrt(2)',
        'strict_rational_Y_squared_upper_bound':'165',
        'full_unconstrained_N_norm_upper_bound':'1'}


def matrix_certificate():
    data=matrices()
    return {'complete_spin_half_source_dimension':32,
        'exact_source_minimal_polynomial':str(data['minimal_polynomial'].as_expr()),
        'exact_scalar_star_and_current_identities':True,
        'matrices':{name:[[str(v) for v in row] for row in data[name].tolist()]
                    for name in ('H','G','N','Y','current')}}
