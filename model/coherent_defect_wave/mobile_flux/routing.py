"""Route the actual paired-flux object using a declared two-lane spinor coin.

An absent lane is the flat Gauss orbit. An occupied lane contains one of
five odd Gauss orbits and a two-component spin register. The sparse one-object
key records the cell, lane/spin coin and literal canonical odd link word.
"""
from functools import lru_cache
import sympy as s

from ..spatial.dirac_route import clifford,route_step,symbol
from .model import observer,representative,charge,STAR


def expanded(value):return s.expand(value)


def odd_words():
    return tuple(sorted({representative(g) for g in observer().base.link_words if charge(g)}))


def validate(state,length):
    if type(length) is not int or length<3:raise ValueError('periodic cube side at least three required')
    out={}
    for (position,coin,word),value in state.items():
        if len(position)!=3 or any(type(x) is not int or not 0<=x<length for x in position):
            raise ValueError('position on the declared periodic cube required')
        if type(coin) is not int or coin not in range(4) or word not in odd_words():
            raise ValueError('two-lane spinor and literal odd orbit required')
        value=s.sympify(value)
        if value.has(s.Float) or value.is_finite is not True or value.is_number is not True:
            raise ValueError('finite exact amplitude required')
        if value!=0:out[tuple(position),coin,word]=value
    return out


def local_code(word):
    """Literal positive Gauss orbit; all other source lanes use local_code(0)."""
    if word not in (0,*odd_words()):raise ValueError('canonical flat/odd word required')
    return {word:1/s.sqrt(2),word^STAR:1/s.sqrt(2)}


def matrix_at_cells(state,matrix):
    out={}
    for (position,column,word),value in state.items():
        for row in range(4):
            key=position,row,word
            out[key]=expanded(out.get(key,0)+matrix[row,column]*value)
    return {key:value for key,value in out.items() if value!=0}


@lru_cache(maxsize=3)
def lane_basis(axis):
    alpha,_=clifford();columns=[]
    for sign in (1,-1):
        projector=(s.eye(4)+sign*alpha[axis])/2
        group=[]
        for index in range(4):
            v=projector[:,index]
            for old in group:v-=old*(old.H*v)[0]
            norm=s.simplify((v.H*v)[0])
            if norm!=0:group.append((v/s.sqrt(norm)).applyfunc(s.simplify))
        if len(group)!=2:raise ArithmeticError('each routing lane must have spin dimension two')
        columns.extend(group)
    result=s.Matrix.hstack(*columns)
    assert result.H*result==s.eye(4)
    assert result*s.diag(1,1,-1,-1)*result.H==alpha[axis]
    return result


def lane_exchange(state):
    return {(position,coin^2,word):value for (position,coin,word),value in state.items()}


def cross_lane_exchange(state,length,axis):
    """SWAP tensor lanes (x,1) and (x+e_axis,0), including their flux words."""
    out={}
    for (position,coin,word),value in state.items():
        lane,spin=divmod(coin,2);target=list(position)
        target[axis]=(target[axis]+(1 if lane else -1))%length
        out[tuple(target),2*(1-lane)+spin,word]=value
    return out


def object_step(state,length,spacing,inverse=False):
    """Execute local rotations and two SWAP layers, with mass angle a*D=2a."""
    state=validate(state,length);a=s.sympify(spacing)
    if a.has(s.Float) or a.is_number is not True or a.is_real is not True or a.is_positive is not True:
        raise ValueError('positive finite exact spacing/tick duration required')
    if type(inverse) is not bool:raise ValueError('boolean inverse flag required')
    _,beta=clifford()
    # The source count is two actual nontrivial face holonomies per object.
    cost={observer().base.law.syndrome(word).bit_count() for word in odd_words()}
    if cost!={2}:raise ArithmeticError('source face-count mass changed')
    coin=s.cos(2*a)*s.eye(4)-s.I*s.sin(2*a)*beta
    if inverse:state=matrix_at_cells(state,coin.H)
    for axis in ((2,1,0) if inverse else range(3)):
        basis=lane_basis(axis);state=matrix_at_cells(state,basis.H)
        if inverse:state=lane_exchange(cross_lane_exchange(state,length,axis))
        else:state=cross_lane_exchange(lane_exchange(state),length,axis)
        state=matrix_at_cells(state,basis)
    if not inverse:state=matrix_at_cells(state,coin)
    return {key:s.simplify(value) for key,value in state.items() if s.simplify(value)!=0}


def position_records(state,length):
    """Complete PVM of the cell containing the odd paired-flux orbit."""
    from itertools import product
    state=validate(state,length)
    return {position:s.simplify(sum((s.conjugate(v)*v for (x,_,_),v in state.items() if x==position),s.Integer(0)))
            for position in product(range(length),repeat=3)}


def certificate():
    from itertools import product
    alpha,beta=clifford();p=s.symbols('p0:3',real=True)
    h=2*beta+sum((p[i]*alpha[i] for i in range(3)),s.zeros(4))
    assert (h*h-(sum(v*v for v in p)+4)*s.eye(4)).applyfunc(s.expand)==s.zeros(4)
    checked=0
    for word,coin in product(odd_words(),range(4)):
        state={((1,1,1),coin,word):s.Integer(1)}
        out=object_step(state,4,s.pi/12)
        reference=route_step({((1,1,1),coin):1},4,s.pi/6)
        assert out=={(x,c,word):s.simplify(v) for (x,c),v in reference.items()}
        assert object_step(out,4,s.pi/12,True)==state
        assert sum(position_records(out,4).values())==1
        checked+=1
    # Literal link-code checks establish that position labels are source charge.
    for word in (0,*odd_words()):
        orbit=local_code(word)
        assert s.simplify(sum(s.conjugate(v)*v for v in orbit.values()))==1
        assert all(charge(g)==int(word!=0) for g in orbit)
        assert {g^STAR:v for g,v in orbit.items()}==orbit
    return {'source':'coherent-mobile-flux-spinor-route-v1',
            'clock':'spacing a=tick duration a, c=1; three lane routes then exp(-i a D beta)',
            'physical_lane_dimension':11,'full_two_lane_local_dimension':121,
            'two_lane_dimensions_refer_to_invariant_logical_code':True,
            'ambient_link_times_three_state_memory_dimension_per_lane':36,
            'ambient_two_lane_dimension_with_empty_native_waves':1296,
            'vacuum_plus_one_object_local_dimension':21,
            'one_object_code':'ell2(cells) tensor C4_lane_spin tensor C5_flux_face',
            'literal_link_code_orbits_checked':6,'exact_source_route_intertwiner_columns':checked,
            'mass_angle':'2*a','source_face_count':2,
            'continuum_generator':'alpha.p + 2*beta',
            'continuum_energy_branches':'+-sqrt(|p|^2+4)',
            'positive_branch_rest_parameter':2,
            'fixed_band_product_formula_error':'t*a*sum_(r<s)|h_r*h_s|; h=(px,py,pz,2)',
            'compact_band_error':'t*a*(3*P^2+6*P) when |p_i|<=P',
            'finite_causal_support':'at most one cell per coordinate per tick',
            'object_position_is_actual_regional_odd_charge_cell':True,
            'coin_mass_cost_and_routing_are_declared_architecture':True,
            'CAR_statistics_of_flux_objects_derived':False,
            'interacting_relativistic_limit_proved':False,
            'positive_energy_many_body_vacuum_derived':False,
            'SM_identity_or_GeV_mass_derived':False}
