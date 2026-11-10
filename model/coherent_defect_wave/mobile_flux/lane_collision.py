"""Complete sparse two-lane controls and an exact coincident-flux collision."""
from itertools import product
import sympy as s

from ..spatial.dirac_route import clifford


def countershift(configuration,inverse=False):
    """Permutation of every lane tensor factor on a periodic one-axis row."""
    values=tuple(configuration)
    if len(values)<6 or len(values)%2 or any(type(v) is not int or v not in range(11) for v in values):
        raise ValueError('at least three cells with two eleven-state lanes required')
    length=len(values)//2;direction=-1 if inverse else 1
    return tuple(values[2*((x-direction)%length)] if lane==0
                 else values[2*((x+direction)%length)+1]
                 for x in range(length) for lane in (0,1))


def local_coin(state,cell,matrix):
    """W on each one-object C20 block; identity on N=0 and N=2."""
    out={}
    for configuration,value in state.items():
        if type(cell) is not int or cell not in range(len(configuration)//2):
            raise ValueError('cell within the configuration required')
        occupied=[lane for lane in (0,1) if configuration[2*cell+lane]]
        if len(occupied)!=1:
            out[configuration]=s.expand(out.get(configuration,0)+value)
            continue
        lane=occupied[0];flavor,spin=divmod(configuration[2*cell+lane]-1,2)
        column=2*lane+spin
        for row in range(4):
            if matrix[row,column]!=0:
                target=list(configuration);target[2*cell]=target[2*cell+1]=0
                target[2*cell+row//2]=1+2*flavor+row%2
                key=tuple(target)
                out[key]=s.expand(out.get(key,0)+matrix[row,column]*value)
    return {key:s.simplify(value) for key,value in out.items() if s.simplify(value)!=0}


def certificate():
    _,beta=clifford();coin=s.sqrt(3)*s.eye(4)/2-s.I*beta/2
    assert (coin.H*coin).applyfunc(s.simplify)==s.eye(4)
    # Every local basis state, including all 100 coincident-object states.
    for labels in product(range(11),repeat=2):
        original={labels:s.Integer(1)};image=local_coin(original,0,coin)
        assert local_coin(image,0,coin.H)==original
        assert s.simplify(sum(s.conjugate(v)*v for v in image.values()))==1
    # All occupation patterns with two representative nonempty labels.
    for configuration in product((0,1,4),repeat=6):
        image=countershift(configuration)
        assert countershift(image,True)==configuration
        assert sorted(image)==sorted(configuration)
    initial=(1,0,0,0,0,4)
    collision=countershift(initial)
    assert collision==(0,0,1,4,0,0)
    state={collision:s.Integer(1)}
    assert local_coin(state,1,coin)==state
    assert countershift(collision,True)==initial
    return {'local_code_dimension':121,'complete_local_coin_basis_columns_checked':121,
            'two_object_local_basis_states_retained':100,
            'row_permutation_configurations_checked':729,
            'incoming_two_object_lane_labels':list(initial),
            'coincident_two_object_lane_labels':list(collision),
            'outcome_probability':'1','coincident_coin_is_declared_identity':True,
            'inverse_and_object_number_checked':True,
            'CAR_or_interacting_continuum_claim':False}
