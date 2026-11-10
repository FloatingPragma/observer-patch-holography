"""Exact Gauss-frame native action and complete routing with fixed native registers.

The weak exponential is defined by a bounded self-adjoint local generator.
Its Taylor vectors are approximations with a stated error, never normalized
and substituted for the exact unitary output.
"""
from dataclasses import dataclass
from functools import lru_cache
from itertools import product
import sympy as s

from ..mobile_flux.model import observer,representative,charge,STAR
from ..mobile_flux.routing import odd_words,lane_basis
from ..mobile_flux.lane_collision import local_coin
from ..species.fock import FockState,MODES,coordinates,mode
from ..spatial.dirac_route import clifford

PROGRAM_ID='coherent-native-routed-paired-flux-v1'


def native_sign(occupation):
    return (-1)**sum(coordinates(i)[0]==5 and coordinates(i)[1]>=3 for i in occupation)


def label_word(label):
    if type(label) is not int or label not in range(11):
        raise ValueError('one of eleven complete lane labels required')
    return 0 if label==0 else odd_words()[(label-1)//2]


def word_label(word,memory):
    rep=representative(word)
    if charge(word)==0:
        if memory!=0:raise ArithmeticError('flat flux requires empty routing memory')
        return 0
    if memory not in (1,2):raise ArithmeticError('odd flux requires occupied routing memory')
    return 1+2*odd_words().index(rep)+memory-1


def label_memory(label):
    return 0 if label==0 else 1+(label-1)%2


def encode_lane(label,occupation):
    word=label_word(label);occupation=tuple(occupation)
    return FockState({(word,occupation):1,(word^STAR,occupation):native_sign(occupation)},2)


@lru_cache(maxsize=4096)
def native_column(label,occupation,inverse=False):
    """Both gauge partners execute U; decoding is checked by full re-encoding."""
    occupation=tuple(occupation);memory=label_memory(label)
    source=encode_lane(label,occupation)
    image=observer().inverse_step(source) if inverse else observer().step(source)
    output={}
    for (word,wedge),value in image.amplitudes:
        if word==representative(word):
            output[word_label(word,memory),wedge]=s.sympify(value)
    encoded={}
    for (target,wedge),value in output.items():
        assert int(target!=0)==int(label!=0)
        assert len(wedge)==len(occupation)
        for key,coefficient in encode_lane(target,wedge).amplitudes:
            encoded[key]=s.expand(encoded.get(key,0)+value*coefficient)
    assert FockState(encoded,2)==image
    assert observer().gamma(image)==image
    return tuple(sorted(output.items()))


@dataclass(frozen=True)
class State:
    amplitudes:tuple
    norm_denominator:int=1

    def __post_init__(self):
        if type(self.norm_denominator) is not int or self.norm_denominator<=0:
            raise ValueError('positive integer normalization required')
        data=self.amplitudes.items() if isinstance(self.amplitudes,dict) else self.amplitudes
        out={};size=None
        for (labels,natives),value in data:
            labels=tuple(labels);natives=tuple(tuple(n) for n in natives)
            if len(labels)!=len(natives) or not labels or len(labels)%2:
                raise ValueError('two routed lanes and two independent native registers per cell required')
            if size is None:size=len(labels)
            if len(labels)!=size:raise ValueError('one fixed finite domain required')
            for label in labels:label_word(label)
            for n in natives:
                if tuple(sorted(set(n)))!=n or any(type(i) is not int or i not in range(MODES) for i in n):
                    raise ValueError('ordered native CAR occupation required')
            value=s.sympify(value)
            if value.is_number is not True or value.is_finite is not True or value.has(s.Float):
                raise ValueError('exact finite amplitudes required')
            key=labels,natives;out[key]=s.expand(out.get(key,0)+value)
        object.__setattr__(self,'amplitudes',tuple(sorted((key,value) for key,value in out.items() if value!=0)))

    def as_dict(self):return dict(self.amplitudes)

    def norm_squared(self):
        return s.simplify(sum((s.conjugate(value)*value for _,value in self.amplitudes),s.Integer(0))/self.norm_denominator)


def combine(parts,coefficients):
    if not parts:raise ValueError('at least one state required')
    denominator=parts[0].norm_denominator
    if any(part.norm_denominator!=denominator for part in parts):raise ValueError('common normalization required')
    out={}
    for part,coefficient in zip(parts,coefficients,strict=True):
        for key,value in part.amplitudes:
            out[key]=s.expand(out.get(key,0)+coefficient*value)
    return State(out,denominator)


def local_native(state,lane,inverse=False):
    out={}
    for (labels,natives),value in state.amplitudes:
        if type(lane) is not int or lane not in range(len(labels)):raise ValueError('native lane index required')
        for (label,occupation),coefficient in native_column(labels[lane],natives[lane],inverse):
            newlabels=list(labels);newlabels[lane]=label
            newnatives=list(natives);newnatives[lane]=occupation
            key=tuple(newlabels),tuple(newnatives)
            out[key]=s.expand(out.get(key,0)+value*coefficient)
    return State(out,state.norm_denominator)


def native_generator(state):
    """K_total=sum_lanes (I-Re U_F), skipping the exact native vacuum."""
    lanes=sorted({lane for (labels,natives),_ in state.amplitudes
                  for lane,n in enumerate(natives) if n})
    parts=[];coefficients=[]
    for lane in lanes:
        parts.extend((state,local_native(state,lane),local_native(state,lane,True)))
        coefficients.extend((1,-s.Rational(1,2),-s.Rational(1,2)))
    return combine(tuple(parts),tuple(coefficients)) if parts else State({},state.norm_denominator)


def coin_at_cell(state,cell,matrix):
    out={}
    for (labels,natives),value in state.amplitudes:
        for target,coefficient in local_coin({labels:value},cell,matrix).items():
            key=target,natives;out[key]=s.expand(out.get(key,0)+coefficient)
    return State(out,state.norm_denominator)


def coin_layer(state,matrix):
    if not state.amplitudes:return state
    cells=len(state.amplitudes[0][0][0])//2
    for cell in range(cells):state=coin_at_cell(state,cell,matrix)
    return state


@lru_cache(maxsize=12)
def lane_permutation(length,axis,inverse=False):
    if type(length) is not int or length<3 or axis not in range(3):
        raise ValueError('periodic cube of side at least three and spatial axis required')
    cells=tuple(product(range(length),repeat=3));lookup={x:i for i,x in enumerate(cells)}
    direction=-1 if inverse else 1
    result=[]
    for x in cells:
        for lane in (0,1):
            source=list(x);source[axis]=(source[axis]+direction*(1 if lane else -1))%length
            result.append(2*lookup[tuple(source)]+lane)
    return tuple(result)


def countertranslation(state,length,axis,inverse=False):
    permutation=lane_permutation(length,axis,inverse)
    out={}
    for (labels,natives),value in state.amplitudes:
        if len(labels)!=len(permutation):raise ValueError('state on the specified cube required')
        target=tuple(labels[i] for i in permutation)
        out[target,natives]=value
    return State(out,state.norm_denominator)


def route(state,length,spacing,inverse=False):
    a=s.sympify(spacing)
    if a.has(s.Float) or a.is_number is not True or a.is_real is not True or a.is_positive is not True:
        raise ValueError('exact positive spacing required')
    _,beta=clifford();mass=s.cos(2*a)*s.eye(4)-s.I*s.sin(2*a)*beta
    if inverse:state=coin_layer(state,mass.H)
    for axis in ((2,1,0) if inverse else range(3)):
        basis=lane_basis(axis)
        state=coin_layer(state,basis.H)
        state=countertranslation(state,length,axis,inverse)
        state=coin_layer(state,basis)
    if not inverse:state=coin_layer(state,mass)
    return state


def weak_taylor(state,epsilon,order=1):
    """Approximate exp(-i epsilon K_total), with a separate certified error."""
    if type(order) is not int or order<0:raise ValueError('nonnegative Taylor order required')
    epsilon=s.sympify(epsilon)
    if epsilon.has(s.Float) or epsilon.is_real is not True or epsilon.is_number is not True:
        raise ValueError('exact real weak duration required')
    parts=[state];coefficients=[s.Integer(1)];power=state
    for j in range(1,order+1):
        power=native_generator(power);parts.append(power)
        coefficients.append((-s.I*epsilon)**j/s.factorial(j))
    return combine(tuple(parts),tuple(coefficients))


def fixed_native_number(state):
    numbers={sum(map(len,natives)) for (_,natives),_ in state.amplitudes}
    if len(numbers)!=1:raise ValueError('a nonzero state in one fixed native-number sector required')
    return numbers.pop()


def taylor_error(state,epsilon,order=1):
    """Operator-norm error on the input's fixed native-number sector."""
    n=fixed_native_number(state)
    return s.Abs(s.sympify(epsilon)*2*n)**(order+1)/s.factorial(order+1)


def joined_taylor(state,length,spacing,coupling,order=1):
    if state.norm_squared()!=1:raise ValueError('normalized input required for this vector error')
    routed=route(state,length,spacing)
    epsilon=s.sympify(spacing)*s.sympify(coupling)
    return weak_taylor(routed,epsilon,order),taylor_error(state,epsilon,order)


def complete_binary(state,predicate):
    """Complete projectors; a Taylor input gives approximate squared norms."""
    return tuple((yes,State({key:value for key,value in state.amplitudes
                            if bool(predicate(key))==yes},state.norm_denominator))
                 for yes in (False,True))


def seed(length=3,with_native=True):
    """A paired flux at (1,1,1); a receiving native wave at (2,2,2), lane0."""
    cells=tuple(product(range(length),repeat=3));lookup={x:i for i,x in enumerate(cells)}
    labels=[0]*(2*len(cells));natives=[()]*(2*len(cells))
    labels[2*lookup[1,1,1]]=1
    if with_native:natives[2*lookup[2,2,2]]=(mode(7,3,0,0),)
    return State({(tuple(labels),tuple(natives)):1})


def changed_face(key):
    labels,_=key
    return any(label>2 for label in labels)
