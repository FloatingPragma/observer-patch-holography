"""Complete CAR matching transport in the occupation-dependent Gauss code."""
from itertools import product
from functools import lru_cache
import sympy as s

from ..native_routing_join.model import (State,route as complete_route,weak_taylor,taylor_error,
    native_generator,seed,changed_face,label_word,word_label,native_sign,MODES,mode)
from ..mobile_flux.model import PairState as LiteralPair,STAR
from ..mobile_flux.routing import object_step
from ..species.fock import wedge_insert
from ..species.spatial_fock import frame_sign

PROGRAM_ID='coherent-propagating-native-flux-v1'
A=(1-s.I)/2
B=(1+s.I)/2


@lru_cache(maxsize=256)
def flux_label_column(labels,length,spacing,inverse=False):
    """Use the proved one-object isometry; retain complete routing otherwise."""
    a=s.sympify(spacing)
    if type(length) is not int or length<3 or a.has(s.Float) or a.is_number is not True or a.is_positive is not True:
        raise ValueError('cube side at least three and exact positive spacing required')
    occupied=[(i,label) for i,label in enumerate(labels) if label]
    if len(labels)!=2*length**3:raise ValueError('cube lane configuration required')
    if not occupied:return ((labels,s.Integer(1)),)
    if len(occupied)==1:
        lane,label=occupied[0]
        cells=tuple(product(range(length),repeat=3));lookup={x:i for i,x in enumerate(cells)}
        source={(cells[lane//2],2*(lane%2)+(label-1)%2,label_word(label)):1}
        image=object_step(source,length,spacing,inverse)
        result=[]
        for (position,spin,word),value in image.items():
            target=[0]*len(labels);target[2*lookup[position]+spin//2]=word_label(word,1+spin%2)
            result.append((tuple(target),value))
        return tuple(sorted(result))
    empty=((),)*len(labels)
    image=complete_route(State({(labels,empty):1}),length,spacing,inverse)
    return tuple((target,value) for (target,_),value in image.amplitudes)


def route(state,length,spacing,inverse=False):
    out={}
    for (labels,natives),value in state.amplitudes:
        for target,weight in flux_label_column(labels,length,spacing,inverse):
            key=target,natives;out[key]=s.expand(out.get(key,0)+value*weight)
    return State(out,state.norm_denominator)


def flattened(natives):
    return tuple(lane*MODES+i for lane,n in enumerate(natives) for i in n)


def separated(occupied,lanes):
    values=[[] for _ in range(lanes)]
    for i in occupied:
        lane,local=divmod(i,MODES);values[lane].append(local)
    return tuple(tuple(n) for n in values)


def edge(state,left,right,inverse=False):
    """Exterior root-SWAP of native modes; flux/routing labels stay at cells."""
    if type(left) is not int or type(right) is not int or left==right:
        raise ValueError('two distinct native register indices required')
    aa,bb=(s.conjugate(A),s.conjugate(B)) if inverse else (A,B)
    output={}
    for (labels,natives),value in state.amplitudes:
        if left not in range(len(labels)) or right not in range(len(labels)):
            raise ValueError('native register within source required')
        terms={():value}
        for index in flattened(natives):
            lane,local=divmod(index,MODES)
            column=((index,1),) if lane not in (left,right) else (
                (index,aa),((right if lane==left else left)*MODES+local,bb))
            following={}
            for wedge,amplitude in terms.items():
                for target,weight in column:
                    inserted,parity=wedge_insert(wedge,target)
                    if parity:following[inserted]=s.expand(following.get(inserted,0)+parity*weight*amplitude)
            terms={k:v for k,v in following.items() if v!=0}
        for wedge,amplitude in terms.items():
            key=labels,separated(wedge,len(labels))
            output[key]=s.expand(output.get(key,0)+amplitude)
    return State(output,state.norm_denominator)


@lru_cache(maxsize=24)
def matching_edges(length,phase):
    if type(length) is not int or length<4 or length%2:
        raise ValueError('even periodic cube side at least four required')
    if type(phase) is not int or phase not in range(6):raise ValueError('six-state matching phase required')
    axis,parity=divmod(phase,2)
    cells=tuple(product(range(length),repeat=3));lookup={x:i for i,x in enumerate(cells)}
    edges=[]
    for x in cells:
        if x[axis]%2!=parity:continue
        y=list(x);y[axis]=(y[axis]+1)%length
        for lane in (0,1):edges.append((2*lookup[x]+lane,2*lookup[tuple(y)]+lane))
    assert sorted(i for pair in edges for i in pair)==list(range(2*length**3))
    return tuple(edges)


def native_matching(state,length,phase,inverse=False):
    pairs=matching_edges(length,phase)
    for left,right in (reversed(pairs) if inverse else pairs):
        if any(natives[left] or natives[right] for (_,natives),_ in state.amplitudes):
            state=edge(state,left,right,inverse)
    return state


def tick_taylor(state,length,spacing,coupling,phase,order=1):
    if state.norm_squared()!=1:raise ValueError('normalized input required for this vector error')
    transported=native_matching(route(state,length,spacing),length,phase)
    epsilon=s.sympify(spacing)*s.sympify(coupling)
    return weak_taylor(transported,epsilon,order),taylor_error(state,epsilon,order),(phase+1)%6


def moving_seed():
    state=seed(4,True);out={}
    cells=tuple(product(range(4),repeat=3));lookup={x:i for i,x in enumerate(cells)}
    receiver=2*lookup[2,2,2];sender=2*lookup[3,2,2]
    for (labels,natives),value in state.amplitudes:
        moved=list(natives);moved[sender]=moved[receiver];moved[receiver]=()
        out[labels,tuple(moved)]=value
    return State(out,state.norm_denominator),receiver,sender


def complete_records(state,receiver):
    return tuple(((present,changed),State({key:value for key,value in state.amplitudes
                    if bool(key[1][receiver])==present and changed_face(key)==changed},state.norm_denominator))
                 for present in (False,True) for changed in (False,True))


def encode_pair(state):
    """Literal two-link expansion for independent verification of CAR dressing."""
    out={}
    if len({labels for (labels,_),_ in state.amplitudes})>1:
        raise ValueError('fixed routing-memory labels required by this literal comparison')
    for (labels,natives),value in state.amplitudes:
        if len(labels)!=2:raise ValueError('two native registers required')
        for bits in product((0,1),repeat=2):
            words=tuple(label_word(labels[i])^(STAR if bits[i] else 0) for i in (0,1))
            coefficient=value*product_sign(native_sign(natives[i]) if bits[i] else 1 for i in (0,1))
            out[words,flattened(natives)]=coefficient
    return LiteralPair(out,4*state.norm_denominator)


def product_sign(values):
    result=1
    for value in values:result*=value
    return result


def literal_pair_edge(state,inverse=False):
    aa,bb=(s.conjugate(A),s.conjugate(B)) if inverse else (A,B)
    output={}
    for (links,occupied),value in state.amplitudes:
        terms={():value}
        for index in occupied:
            lane,local=divmod(index,MODES);other=1-lane
            dressing=frame_sign(links[lane],local)*frame_sign(links[other],local)
            following={}
            for wedge,amplitude in terms.items():
                for target,weight in ((index,aa),(other*MODES+local,bb*dressing)):
                    inserted,parity=wedge_insert(wedge,target)
                    if parity:following[inserted]=s.expand(following.get(inserted,0)+parity*weight*amplitude)
            terms={k:v for k,v in following.items() if v!=0}
        for wedge,amplitude in terms.items():
            output[links,wedge]=s.expand(output.get((links,wedge),0)+amplitude)
    return LiteralPair(output,state.norm_denominator)
