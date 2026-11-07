"""Complete E-address proof of the reflected four-read correlation probe."""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations, product

from . import sectors as s


def field(occupied, first, second):
    return (int(first[0] in occupied)-int(first[1] in occupied))*(
        int(second[0] in occupied)-int(second[1] in occupied))


def raw_maps(edge, ports=4, masks=None):
    """All25 fixed-seam slots of an E command, including zero C/hold slots."""
    s.integer(ports,2)
    if (type(edge) is not tuple or len(edge)!=2 or any(type(i) is not int for i in edge)
            or not 0<=edge[0]<edge[1]<ports):
        raise ValueError('ordered native endpoints required')
    if masks is None:
        if ports>12:
            raise ValueError('full binary raw-map evaluation is bounded to12 ports')
        masks=tuple(range(1 << ports))
    elif (type(masks) is not tuple or len(masks)>10000 or len(set(masks))!=len(masks)
          or any(type(i) is not int or i<0 or i.bit_length()>ports for i in masks)):
        raise ValueError('bounded distinct binary source masks required')
    result=[(('hold',),F(0),{},0,0)]
    for law in ('C','E'):
        for coin in (0,1):
            decoding=s.source.interpreter.circuit.record_values(law,coin)
            for raw in range(1 << s.source.interpreter.circuit.RECORD_BITS[law]):
                label=decoding.get(raw); mapping={}
                if law=='E' and label is not None:
                    for mask in masks:
                        target,current,contrast=s.source.image(mask,edge,coin,ports)
                        if contrast==label:
                            mapping[mask]=target
                    current=(label-(2*coin-1)*(abs(label)%2))//2
                    result.append(((law,coin,raw),F(1,2),mapping,label,label-2*current))
                else:
                    result.append(((law,coin,raw),F(0),{},0,0))
    if len(result)!=25:
        raise ArithmeticError('the complete native E alphabet changed')
    return tuple(result)


def pair_certificate():
    first=(0,1); second=(2,3)
    maps1=raw_maps(first);maps2=raw_maps(second)
    size=16;identity=[F(0)]*size;unital=[F(0)]*size
    incoming=[F(0)]*size;outgoing=[F(0)]*size;live=0
    # All pair maps are partial injections. This proves the diagonal matrix
    # identities on the full16-dimensional operator space and all references.
    for a,b in product(maps1,maps2):
        weight=a[1]*b[1]
        image={i:b[2][j] for i,j in a[2].items() if j in b[2]}
        if len(set(image.values()))!=len(image):
            raise ArithmeticError('E pair map is not a partial injection')
        if weight and image:live+=1
        sin=a[3]*b[3];sout=a[4]*b[4]
        for i,j in image.items():
            fi=field({k for k in range(4) if i >> k & 1},first,second)
            fj=field({k for k in range(4) if j >> k & 1},first,second)
            if fi!=sin or fj!=sout:
                raise ArithmeticError('native contrast record does not match its time cut')
            identity[i]+=weight;unital[j]+=weight
            incoming[i]+=weight*sin;outgoing[j]+=weight*sout
    diagonal=[field({k for k in range(4) if i >> k & 1},first,second) for i in range(size)]
    if identity!=[1]*size or unital!=[1]*size or incoming!=diagonal or outgoing!=diagonal:
        raise ArithmeticError('complete reflected instrument identity failed')
    return {'one_E_raw_slots':25,'two_E_raw_histories':625,'live_two_E_histories':live,
            'zero_two_E_histories':625-live,'full_binary_basis':16,
            'K_dagger_K_sum':'I16','K_K_dagger_sum':'I16',
            'incoming_scored_effect':'diag(phi)','outgoing_scored_identity_image':'diag(phi)',
            'scalar_scores':{'past':'(d1-2*j1)*(d2-2*j2)','future':'d3*d4'},
            'all_untouched_references_and_exterior_identity_extensions':True}


def native_step(occupied, edge, coin):
    ports=max((*occupied,*edge),default=0)+1
    mask=sum(1 << i for i in occupied)
    target,current,contrast=s.source.image(mask,edge,coin,max(2,ports))
    return tuple(i for i in range(ports) if target >> i & 1),current,contrast


def direct_probe(family, occupied, steps, first, second, *, reflected=True):
    """Independent complete native coin enumeration on bounded control graphs.

    Zero raw addresses have exactly zero weight by pair_certificate; no
    outcome conditioning or renormalization is performed here.
    """
    s.fair_graph(family);ports=len(family.weights)
    s.degree_data(ports,occupied);s.integer(steps)
    s.contrast(first,second,ports)
    if first not in family.seams or second not in family.seams or ports>8 or steps>3:
        raise ValueError('bounded native control graph and existing probe seams required')
    if type(reflected) is not bool:
        raise ValueError('a literal reflection control required')
    states=tuple(combinations(range(ports),occupied));table=defaultdict(F)
    for state in states:
        for coin1,coin2 in product((0,1),repeat=2):
            y,j1,d1=native_step(state,first,coin1)
            z,j2,d2=native_step(y,second,coin2)
            score=(d1-2*j1)*(d2-2*j2) if reflected else d1*d2
            table[(z,score)]+=F(1,4*len(states))
    for _ in range(steps):
        nxt=defaultdict(F)
        for (state,score),weight in table.items():
            for edge,coin in product(family.seams,(0,1)):
                y,_,_=native_step(state,edge,coin)
                nxt[(y,score)]+=weight/F(2*len(family.seams))
        table=nxt
    value=F(0);total=F(0)
    for (state,score),weight in table.items():
        for coin1,coin2 in product((0,1),repeat=2):
            y,_,d1=native_step(state,first,coin1)
            _,_,d2=native_step(y,second,coin2)
            value+=weight*score*d1*d2/4;total+=weight/4
    if total!=1:
        raise ArithmeticError('complete control lost probability')
    return value


def probe_program(level,occupied,steps):
    """Admit the entire four-E experiment with exported complete records."""
    graph,first,second,_=s.probe_data(level)
    s.degree_data(len(graph.weights),occupied);s.integer(steps)
    if steps>s.source.interpreter.MAX_ATTEMPTS-4:
        raise ValueError('complete probe exceeds the native attempt budget')
    pair=[{'op':'repair','selector':{'kind':'fixed','seam':list(e)},'law':'E','repeat':1}
          for e in (first,second)]
    body=([{'op':'repair','selector':{'kind':'uniform'},'law':'mixture','repeat':steps}]
          if steps else [])
    return s.source.interpreter.program(federation={'level':level,'total':occupied},
                                        commands=pair+body+pair,overflow='export_oldest')
