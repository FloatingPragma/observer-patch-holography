"""Native cross-probe programs and final-tree access to occupation snapshots."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations,product
from model.native_correlation_sectors import records as parent
from model.native_sector_readout import readout as tree
from model.native_local_response import cone
from . import response as r

source=tree.source


def program(level,occupied,attempts,past,future):
    family=r.s.source.federation_family(level,(1,)*12)
    r.probe(family,past);r.probe(family,future)
    r.s.norm_factor(len(family.weights),2,occupied);r.s.integer(attempts)
    if attempts>source.MAX_ATTEMPTS-4:raise ValueError('four probes and body exceed native admission')
    fixed=lambda edges:[{'op':'repair','selector':{'kind':'fixed','seam':list(e)},'law':'E','repeat':1} for e in edges]
    body=[{'op':'repair','selector':{'kind':'uniform'},'law':'mixture','repeat':attempts}] if attempts else []
    return source.program(federation={'level':level,'total':occupied},commands=fixed(past)+body+fixed(future),
                          capacity=1,overflow='export_oldest')


def direct_cross(family,occupied,attempts,past,future,*,reflected=True):
    """Independent literal endpoint enumeration on bounded control graphs."""
    r.probe(family,past);r.probe(family,future);ports=len(family.weights)
    r.s.degree_data(ports,occupied);r.s.integer(attempts)
    if ports>8 or attempts>3 or type(reflected) is not bool:
        raise ValueError('bounded complete native control and literal reflection flag required')
    states=tuple(combinations(range(ports),occupied));law=defaultdict(F)
    for state in states:
        for a,b in product((0,1),repeat=2):
            mid,j1,d1=parent.native_step(state,past[0],a)
            after,j2,d2=parent.native_step(mid,past[1],b)
            score=(d1-2*j1)*(d2-2*j2) if reflected else d1*d2
            law[(after,score)]+=F(1,4*len(states))
    for _ in range(attempts):
        nxt=defaultdict(F)
        for (state,score),weight in law.items():
            for edge,coin in product(family.seams,(0,1)):
                after,_,_=parent.native_step(state,edge,coin)
                nxt[(after,score)]+=weight/F(2*len(family.seams))
        law=nxt
    score=F();mass=F()
    for (state,past_score),weight in law.items():
        for a,b in product((0,1),repeat=2):
            mid,_,d1=parent.native_step(state,future[0],a)
            _,_,d2=parent.native_step(mid,future[1],b)
            score+=weight*past_score*d1*d2/4;mass+=weight/4
    if mass!=1:raise ArithmeticError('complete cross-probe control lost probability')
    return score


def snapshot_program(level,occupied,attempts):
    terminal=tree.program(level=level,total=occupied);r.s.integer(attempts)
    if attempts>source.MAX_ATTEMPTS-len(terminal['commands']):
        raise ValueError('body plus terminal tree exceeds native attempt admission')
    body=[{'op':'repair','selector':{'kind':'uniform'},'law':'E','repeat':attempts}] if attempts else []
    return source.program(federation={'level':level,'total':occupied},
                          biased=False,commands=body+terminal['commands'],capacity=1,overflow='export_oldest')


def snapshot_specification(value):
    value=source.admit(value)
    if value['carrier']['kind']!='federation' or not value['commands']:
        raise ValueError('the actual fair federation and terminal detector required')
    level=value['carrier']['level'];occupied=value['carrier']['total']
    body=value['commands'][0]['repeat'] if value['commands'][0]['selector']['kind']=='uniform' else 0
    if source.canonical(value)!=source.canonical(snapshot_program(level,occupied,body)):
        raise ValueError('the exact pure-E body and terminal native tree required')
    return value,cone.source_family(value),occupied,body


def decode_snapshots(value,records):
    value,family,occupied,body=snapshot_specification(value)
    if type(records) is not tuple or len(records)!=body+len(tree.tree(family)):
        raise ValueError('the complete body and final detector tape required')
    for mark in records[:body]:
        if type(mark) is not source.Record or mark.selector!='uniform' or mark.law!='E':
            raise ValueError('live original uniform E records required')
    decoded=tree.decode(family,records[body:],occupied)
    if not decoded:return ()
    initial_terminal,post_detector=decoded[0]
    path=tree.reverse_currents(family,initial_terminal,records[:body])
    if any(x.bit_count()!=occupied for x in path):raise ArithmeticError('reconstructed path changed total occupation')
    return {'body_masks':path,'post_detector_mask':post_detector}


def trace_snapshots(value,initial,edges,coins):
    """One original complete kernel word, not an expanded-history execution."""
    value,family,occupied,body=snapshot_specification(value);ports=len(family.weights)
    if (body>4 or type(initial) is not int or initial<0 or initial.bit_length()>ports
            or initial.bit_count()!=occupied or type(edges) is not tuple or len(edges)!=body
            or any(type(e) is not tuple or len(e)!=2 or any(type(x) is not int for x in e)
                   or e not in family.seams for e in edges)
            or type(coins) is not tuple or len(coins)!=body+len(tree.tree(family))
            or any(type(x) is not int or x not in (0,1) for x in coins)):
        raise ValueError('bounded exact body, native initial mask, seams and complete coin word required')
    command={'op':'repair','selector':{'kind':'uniform'},'law':'E','repeat':1}
    commands=[command]*body+(value['commands'][1:] if body else value['commands'])
    current=initial;density={(initial,0,initial,0):F(1)};memory=archive=();marks=[];slots=[];path=[initial]
    for index,(cmd,coin) in enumerate(zip(commands,coins,strict=True)):
        edge=edges[index] if index<body else tuple(cmd['selector']['seam'])
        following,j,d=tree.native.image(current,edge,coin,ports)
        candidates=[(slot,out) for slot,out in enumerate(source._outcomes(value,cmd,density,memory,archive))
                    if out.record.kind=='repair' and out.record.seam==edge and out.record.law=='E'
                    and out.record.coin==coin and out.record.label==d and out.record.current==j]
        if len(candidates)!=1:raise ArithmeticError('complete native word lost its original address')
        slot,out=candidates[0];slots.append(slot);marks.append(out.record)
        density,memory,archive=out.density,out.memory,out.archive;current=following
        if index<body:path.append(current)
    decoded=decode_snapshots(value,tuple(marks))
    if decoded!={'body_masks':tuple(path),'post_detector_mask':current}:
        raise ArithmeticError('final-tree record reconstruction disagrees with original kernels')
    weight=F(1,2*len(family.seams))**body*F(1,2)**len(tree.tree(family))
    if source.trace(density)!=weight:raise ArithmeticError('complete word probability changed')
    if not tree._same_records(archive+memory,tuple(marks)):
        raise ArithmeticError('actual exported archive and FIFO lost the record tape')
    return {'program_digest':source.digest(value),'records':tuple(marks),'slots':tuple(slots),
            'decoded':decoded,'density':density,'archive':archive,'memory':memory,'weight':weight,
            'status':'COMPLETE_RAW_WORD_QUERY'}
