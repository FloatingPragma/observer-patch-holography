"""Endpoint readback from the complete body tape and the native E tree."""
from fractions import Fraction as F
from dataclasses import fields
from model.native_sector_readout import readout as detector
from model.native_local_response import cone
from . import geometry

source=detector.source


def same(left,right):
    """Type-preserving equality for exact source objects, including dictionary keys."""
    if type(left) is not type(right):return False
    if type(left) is source.Record:
        return all(same(getattr(left,f.name),getattr(right,f.name)) for f in fields(left))
    if type(left) is dict:
        return len(left)==len(right) and all(any(same(k,rk) and same(v,rv) for rk,rv in right.items())
                                               for k,v in left.items())
    if type(left) in (tuple,list):
        return len(left)==len(right) and all(same(a,b) for a,b in zip(left,right,strict=True))
    if type(left) in (int,str,bool,F,type(None)):return left==right
    return False


def program(level,attempts):
    if type(attempts) is not int or attempts<0:raise ValueError('a finite body count required')
    terminal=detector.program(level=level,total=1)
    if attempts+len(terminal['commands'])>min(source.MAX_NODES,source.MAX_ATTEMPTS):
        raise ValueError('body and terminal tree exceed the existing native admission limits')
    body=[{'op':'repair','selector':{'kind':'uniform'},'law':'E','repeat':1} for _ in range(attempts)]
    value=source.program(federation={'level':level,'total':1},biased=False,
                         commands=body+terminal['commands'],capacity=1,overflow='export_oldest')
    return value


def specification(value):
    value=source.admit(value)
    c=value['carrier']
    if c['kind']!='federation' or c['total']!=1 or c['level']>2:
        raise ValueError('an evaluated native one-unit federation required')
    family=cone.source_family(value);length=len(detector.tree(family));body=len(value['commands'])-length
    if body<0 or source.canonical(value)!=source.canonical(program(c['level'],body)):
        raise ValueError('the fixed uniform body and complete terminal E tree must bind exactly')
    return value,family,body


def decode(value,records):
    value,family,body=specification(value)
    if type(records) is not tuple or len(records)!=len(value['commands']):
        raise ValueError('every original body and detector record required')
    for record in records[:body]:
        if (type(record) is not source.Record or record.kind!='repair' or record.selector!='uniform'
                or record.law!='E' or record.seam not in family.seams):
            raise ValueError('a live original uniformly selected E body record required')
    decoded=detector.decode(family,records[body:],1)
    if not decoded:return ()
    before,after=decoded[0]
    path=detector.reverse_currents(family,before,records[:body])
    if any(mask.bit_count()!=1 for mask in path):raise ValueError('decoded path left the one-unit sector')
    return {'body_masks':path,'body_terminal_mask':before,'post_detector_mask':after}


def query(value,slots,density):
    """One unnormalized word of the complete original kernel, with FIFO export."""
    value,_,_=specification(value)
    if type(slots) is not tuple or len(slots)!=len(value['commands']):
        raise ValueError('a full native raw slot word required')
    state=source._validate_density(density,value['carrier']);memory=archive=();records=[]
    for command,slot in zip(value['commands'],slots,strict=True):
        outcomes=source._outcomes(value,command,state,memory,archive)
        if type(slot) is not int or not 0<=slot<len(outcomes):raise ValueError('an original raw slot required')
        selected=outcomes[slot];state,memory,archive=selected.density,selected.memory,selected.archive
        records.append(selected.record)
    return {'program_digest':source.digest(value),'slots':slots,'records':tuple(records),'density':state,
            'archive':archive,'memory':memory,'status':'COMPLETE_RAW_WORD_QUERY'}


def trace(value,initial,edges,coins):
    value,family,body=specification(value);ports=len(family.weights)
    if (type(initial) is not int or not 0<=initial<ports or type(edges) is not tuple
            or len(edges)!=body or any(type(edge) is not tuple or len(edge)!=2
                or any(type(v) is not int for v in edge) or edge not in family.seams for edge in edges)
            or type(coins) is not tuple or len(coins)!=len(value['commands'])
            or any(type(c) is not int or c not in (0,1) for c in coins)):
        raise ValueError('an actual initial port, native body seams and complete coin word required')
    current=1<<initial;state={(current,0,current,0):F(1)};memory=archive=();slots=[];records=[];path=[current]
    for t,(command,coin) in enumerate(zip(value['commands'],coins,strict=True)):
        edge=edges[t] if t<body else tuple(command['selector']['seam'])
        target,j,d=detector.native.image(current,edge,coin,ports)
        outcomes=source._outcomes(value,command,state,memory,archive)
        matches=[(slot,out) for slot,out in enumerate(outcomes) if out.record.kind=='repair'
                 and out.record.seam==edge and out.record.law=='E' and out.record.coin==coin
                 and out.record.label==d and out.record.current==j]
        if len(matches)!=1:raise ArithmeticError('literal body/detector branch is not unique')
        slot,out=matches[0];slots.append(slot);records.append(out.record)
        state,memory,archive=out.density,out.memory,out.archive;current=target
        if t<body:path.append(current)
    decoded=decode(value,tuple(records))
    if decoded['body_masks']!=tuple(path) or decoded['post_detector_mask']!=current:
        raise ArithmeticError('record-only path differs from original endpoint execution')
    expected=F(1,2*len(family.seams))**body*F(1,2)**len(detector.tree(family))
    if source.trace(state)!=expected:raise ArithmeticError('native complete word weight changed')
    return {'program_digest':source.digest(value),'slots':tuple(slots),'records':tuple(records),'density':state,
            'archive':archive,'memory':memory,'decoded':decoded,'weight':expected}


def authenticate(value,result):
    if (type(result) is not dict or set(result)!={'program_digest','slots','records','density','archive','memory','decoded','weight'}
            or result.get('program_digest')!=source.digest(source.admit(value))):
        raise ValueError('trace must bind the exact admitted program')
    if type(result['weight']) is not F:raise ValueError('an exact unnormalized word probability required')
    recovered=decode(value,result['records'])
    if not recovered:raise ValueError('a live decoded trace required for authentication')
    initial=recovered['body_masks'][0]
    actual=query(value,result['slots'],{(initial,0,initial,0):F(1)})
    for key in ('records','density','archive','memory'):
        if not same(actual[key],result[key]):raise ValueError('literal source query differs from supplied trace')
    if not same(recovered,result['decoded']) or source.trace(actual['density'])!=result['weight']:
        raise ValueError('decoded path or complete probability changed')
    return True


def shortest_path(graph,sender,receiver):
    if any(type(v) is not int or not 0<=v<graph.ports for v in (sender,receiver)):
        raise ValueError('actual source endpoints required')
    parents={sender:None};queue=[sender]
    for v in queue:
        if v==receiver:break
        for w in graph.neighbours[v]:
            if w not in parents:parents[w]=v;queue.append(w)
    if receiver not in parents:raise ArithmeticError('native source disconnected')
    path=[receiver]
    while path[-1]!=sender:path.append(parents[path[-1]])
    return tuple(reversed(path))


def arrival_witness(level,receiver=0):
    graph=geometry.graph(level);distance={receiver:0};queue=[receiver]
    for v in queue:
        for w in graph.neighbours[v]:
            if w not in distance:distance[w]=distance[v]+1;queue.append(w)
    sender=max(distance,key=lambda v:(distance[v],v));path=shortest_path(graph,sender,receiver)
    edges=tuple(tuple(sorted((a,b))) for a,b in zip(path,path[1:]));coins=tuple(int(b==min(a,b)) for a,b in zip(path,path[1:]))
    mask=1<<sender
    for index,(edge,coin) in enumerate(zip(edges,coins,strict=True)):
        if (receiver in edge)!=(index==len(edges)-1):raise ArithmeticError('an early receiver tick occurred')
        mask,_,_=detector.native.image(mask,edge,coin,graph.ports)
        if mask!=1<<path[index+1]:raise ArithmeticError('native placement failed the path witness')
    return {'level':level,'sender':sender,'receiver':receiver,'graph_distance':len(edges),
            'edges':edges,'coins':coins,'receiver_ticks':1,
            'specified_complete_body_word_probability':F(1,2*len(graph.seams))**len(edges),
            'expanded_execution_claimed':False}


def certificate():
    value=program(0,2);family=cone.source_family(value)
    path=(2,10,0);edges=tuple(tuple(sorted(e)) for e in zip(path,path[1:]));coins=(0,1)+(0,)*len(detector.tree(family))
    result=trace(value,path[0],edges,coins);authenticate(value,result)
    stopped=source.execute(value,{(1<<path[0],0,1<<path[0],0):F(1)})
    if stopped.status!='RESOURCE_STOP':raise ArithmeticError('expanded native evaluator boundary changed')
    aggregate={}
    for out in stopped.branches:
        for (i,_,j,_),weight in out.density.items():
            if i!=j:raise ArithmeticError('basis preparation created an undeclared coherence')
            aggregate[i]=aggregate.get(i,F())+weight
    graph=geometry.graph(0);M=len(graph.seams)
    expected={1<<p:F(1,2*M) for p in graph.neighbours[path[0]]}
    expected[1<<path[0]]=1-F(len(graph.neighbours[path[0]]),2*M)
    if aggregate!=expected:raise ArithmeticError('the original uniform native prefix differs from I-L/(2M)')
    crossing=graph.inter[0];cross_value=program(0,1)
    cross=trace(cross_value,crossing[0],(crossing,),(0,)+(0,)*len(detector.tree(family)))
    authenticate(cross_value,cross)
    x,denominator=geometry.coordinates(0)
    norm=geometry.normalized_square(tuple(a-b for a,b in zip(x[crossing[0]],x[crossing[1]],strict=True)),denominator)
    if geometry.sign(norm)<=0:raise ArithmeticError('the cross-cell readout has zero comparison displacement')
    arrivals=[]
    for ell in (0,1,2,3):
        witness=arrival_witness(ell)
        if witness['graph_distance']<=2**ell:raise ArithmeticError('finite receiver-distance control collapsed')
        arrivals.append({'level':ell,'positive_complete_first_tick_word_checked':True,
                         'distance_exceeds':2**ell,'receiver_ticks':1,
                         'body_word_probability':'(2M)^(-D)',
                         'expanded_execution_claimed':False})
    return {'actual_body_plus_detector':{'program_digest':source.digest(value),'body_attempts':2,
              'detector_attempts':len(detector.tree(family)),'decoded_body_ports':[m.bit_length()-1 for m in result['decoded']['body_masks']],
              'word_probability':str(result['weight']),'retained_records':len(result['records']),
              'exported_records':len(result['archive']),'FIFO_slots':len(result['memory']),
              'all_body_records_used_for_readback':True,'original_kernel_query_authenticated':True,
              'actual_expanded_prefix_matches_one_unit_transition':True,
              'expanded_status':stopped.status,'expanded_completed_attempts':stopped.attempted_steps},
            'actual_crossing_body_plus_detector':{'program_digest':source.digest(cross_value),
                'squared_comparison_displacement':geometry.display(norm),'weight':str(cross['weight']),
                'body_path_decoded_from_complete_records':True,'original_kernel_query_authenticated':True},
            'first_receiver_tick_witnesses':arrivals,
            'initial_occupation_disclosure_used_by_decoder':False,
            'local_receiver_records_alone_reconstruct_whole_body':False,
            'complete_continued_quantum_state_replaced_by_decoded_classical_path':False}
