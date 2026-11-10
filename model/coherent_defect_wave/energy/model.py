"""Source-defined positive energy and exact one-wave Hamiltonian action.

H=sum_x[1-Re(U_F,x)]+(pi/4)sum_edges dGamma(1-S_edge).
This is a specified continuous clock, not a logarithm of the discrete tick.
"""
from fractions import Fraction as F
from itertools import product
import sympy as s

from ..model import JointObserver, SparseState, combine
from ..species.fock import FockObserver, FockState, mode
from ..species.spatial_fock import PairFockState,MODES,frame_sign,onsite as pair_onsite
from ..cubic_readout import lattice


def local_action(observer,state):
    """A=I-(U+U*)/2 on the literal one-wave source."""
    average=combine(observer.step(state),observer.inverse_step(state),F(1,2),F(1,2))
    return combine(state,average,1,-1)


def fock_local_action(observer,state):
    """The same functional of the complete occupation unitary U_F."""
    result=state.as_dict()
    for branch in (observer.step(state),observer.inverse_step(state)):
        for key,value in branch.amplitudes:
            result[key]=result.get(key,F(0))-value/2
    return FockState(result,state.norm_denominator)


def pair_hopping_laplacian(state):
    """Full exterior dGamma(I-S), with links fixed and CAR signs retained."""
    if not isinstance(state,PairFockState):
        raise ValueError('two-cell Fock state required')
    result={key:len(key[1])*value for key,value in state.amplitudes}
    for (links,occupied),value in state.amplitudes:
        for ordinal,index in enumerate(occupied):
            cell,local=divmod(index,MODES);target=(1-cell)*MODES+local
            rest=occupied[:ordinal]+occupied[ordinal+1:]
            if target in rest:
                continue
            sign=(-1)**(ordinal+sum(i<target for i in rest))
            sign*=frame_sign(links[cell],local)*frame_sign(links[1-cell],local)
            key=links,tuple(sorted((*rest,target)))
            result[key]=s.expand(result.get(key,0)-sign*value)
    return PairFockState(result,state.norm_denominator)


def pair_energy_action(state):
    """Full two-cell occupation H, including fermionic hopping signs."""
    hopping=pair_hopping_laplacian(state)
    result={key:2*value for key,value in state.amplitudes}
    for cell in (0,1):
        for branch in (pair_onsite(state,cell),pair_onsite(state,cell,True)):
            for key,value in branch.amplitudes:
                result[key]=s.expand(result.get(key,0)-value/2)
    for key,value in hopping.amplitudes:
        result[key]=s.expand(result.get(key,0)+s.pi*value/4)
    return PairFockState(result,state.norm_denominator)


def linear_combination(states,weights):
    if not states or len(states)!=len(weights):
        raise ValueError('a nonempty matching list of states and weights is required')
    denominator=states[0].norm_denominator
    if any(state.norm_denominator!=denominator for state in states):
        raise ValueError('one common normalization required')
    result={}
    for state,weight in zip(states,weights,strict=True):
        if type(weight) not in (int,F):
            raise ValueError('exact rational coefficients required')
        for key,value in state.amplitudes:
            result[key]=lattice.plus(result.get(key,lattice.ZERO),lattice.scale(value,F(weight)))
    return lattice.State(result,denominator)


def onsite_energy(state):
    return linear_combination((state,lattice.onsite(state),lattice.onsite(state,True)),
                              (1,F(-1,2),F(-1,2)))


def hopping_laplacian(state,length):
    """One-wave sum dGamma(I-S_edge); every link stays at its cell."""
    lattice.validate(state,length)
    result={key:lattice.scale(value,6) for key,value in state.amplitudes}
    for (position,index,changes),value in state.amplitudes:
        for axis,sign in product(range(3),(-1,1)):
            target=list(position);target[axis]=(target[axis]+sign)%length
            key=tuple(target),index,changes
            result[key]=lattice.plus(result.get(key,lattice.ZERO),lattice.scale(value,-1))
    return lattice.State(result,state.norm_denominator)


def h_action(state,length):
    """Return exact Gaussian-rational coefficients of H psi=Apsi+pi*Bpsi."""
    lattice.validate(state,length)
    return onsite_energy(state),linear_combination((hopping_laplacian(state,length),),(F(1,4),))


def energy_moments(state,length):
    if state.norm_squared()!=1:
        raise ValueError('normalized input required')
    a,b=h_action(state,length)
    ar,ai=lattice.overlap(state,a);br,bi=lattice.overlap(state,b)
    if ai or bi:
        raise ArithmeticError('Hermitian energy acquired an imaginary expectation')
    mean=s.Rational(ar)+s.pi*s.Rational(br)
    cross=lattice.overlap(a,b)[0]
    second=s.Rational(a.norm_squared())+2*s.pi*s.Rational(cross)+s.pi**2*s.Rational(b.norm_squared())
    return {'mean':s.expand(mean),'second':s.expand(second),'variance':s.expand(second-mean**2),
            'A_image':a,'hopping_coefficient_image':b}


def dark_wave(observer,internal=0,word=None):
    """Exact U=-1 wave: port2, blank repair, opposite address amplitudes."""
    if type(internal) is not int or internal not in range(5):
        raise ValueError('one of five internal components required')
    word=observer.initial_word if word is None else word
    a=observer.addresses.index((0,1));b=observer.addresses.index((1,0))
    return SparseState({(word,2,internal,0,a):F(1),
                        (word,2,internal,0,b):F(-1)},2)


def dark_pair(observer):
    """Two distinct internal dark waves on the same shared link register."""
    a=observer.base.addresses.index((0,1));b=observer.base.addresses.index((1,0))
    result={}
    for first,weight0 in ((a,1),(b,-1)):
        for second,weight1 in ((a,1),(b,-1)):
            occupied=tuple(sorted((mode(2,0,0,first),mode(2,1,0,second))))
            result[observer.base.initial_word,occupied]=F(weight0*weight1)
    return FockState(result,4)


def dark_plane_wave(length,momentum,internal=0):
    """Exact side-four plane waves in the Gauss-positive source dark band."""
    if length!=4 or len(momentum)!=3 or any(type(x) is not int or x not in range(4) for x in momentum):
        raise ValueError('three Z4 labels on the side-four torus required')
    if type(internal) is not int or internal not in range(5):
        raise ValueError('one of five internal components required')
    units=((F(1),F(0)),(F(0),F(1)),(F(-1),F(0)),(F(0),F(-1)))
    result={}
    for position in product(range(length),repeat=3):
        phase=units[sum(k*x for k,x in zip(momentum,position,strict=True))%4]
        for address,sign in ((0,1),(1,-1)):
            result[position,mode(2,internal,0,address),()]=lattice.scale(phase,sign)
    return lattice.State(result,2*length**3)
