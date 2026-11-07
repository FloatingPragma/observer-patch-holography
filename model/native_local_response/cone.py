"""Exact backward-cone reduction with complete local raw records.

Remote raw values are traced by an explicit observation map. They remain
in the source archive; a controller reading them is outside this reduction.
"""
from collections import defaultdict
from fractions import Fraction as F
from model.native_program import interpreter as source
from model.weighted_source_family import family as native


def subset(vertices,ports):
    if (type(vertices) is not tuple or vertices!=tuple(sorted(set(vertices)))
            or any(type(v) is not int or not 0<=v<ports for v in vertices)):
        raise ValueError('a canonical source vertex subset required')
    return vertices


def source_family(value):
    """Read an admitted binary source; no weighted federation is introduced."""
    value=source.admit(value);carrier=value['carrier']
    if carrier['kind']=='local' and carrier['shell']=='H66':
        parameters=source._parameters(value['parameters'],carrier)[0]
        return native.local_family(parameters)
    if carrier['kind']=='federation' and carrier['level']<=source.MAX_EVALUATED_LEVEL:
        return native.federation_family(carrier['level'],(1,)*12)
    raise ValueError('an actual binary H66 or evaluated fair federation source required')


def audit_program(value):
    value=source.admit(value);family=source_family(value)
    attempts=sum(command['repeat'] for command in value['commands'])
    if value['memory']['overflow']!='export_oldest' and value['memory']['capacity']<attempts:
        raise ValueError('intrinsic FIFO termination cannot be erased by a spatial reduction')
    for command in value['commands']:
        if command['op']=='repair' and command['selector']['kind'] not in ('fixed','uniform'):
            raise ValueError('record-adaptive selector requires its actual controller dependency cone')
    return {'attempts':attempts,'vertices':len(family.weights),'seams':len(family.seams),
            'selector_reads_FIFO_values':False,'raw_archive_preserved_by_source':True,
            'local_observation_traces_unselected_raw_values':True,
            'host_expansion_completion_asserted':False}


def backward(schedule,final,observed=()):
    """A conservative exact-support cone, not a minimality claim.

    Observed command endpoints are roots at their event times. An unobserved
    event is summed over its complete raw alphabet before it can be removed.
    """
    if type(schedule) is not tuple or any(type(e) is not tuple or len(e)!=2 or
            any(type(v) is not int or v<0 for v in e) or e[0]>=e[1] for e in schedule):
        raise ValueError('a finite ordered endpoint schedule required')
    ports=1+max((*final,*(v for e in schedule for v in e)),default=0)
    subset(final,ports)
    if type(observed) is not tuple or observed!=tuple(sorted(set(observed))) or any(
            type(t) is not int or not 0<=t<len(schedule) for t in observed):
        raise ValueError('canonical observed command indices required')
    support=set(final);kept=[];layers=[tuple(sorted(support))]
    for t in reversed(range(len(schedule))):
        edge=schedule[t]
        if t in observed or support.intersection(edge):
            kept.append(t);support.update(edge)
        layers.append(tuple(sorted(support)))
    return {'initial_support':tuple(sorted(support)),
            'observation_roots':tuple(sorted(set(final)|{v for t in observed for v in schedule[t]})),
            'retained_commands':tuple(reversed(kept)),
            'reverse_supports':tuple(layers),
            'observed_commands':observed}


def partial_trace(density,window,ports):
    subset(window,ports);mask=sum(1<<i for i in window);result=defaultdict(F)
    for (i,a,j,b),value in density.items():
        if (i&~mask)==(j&~mask):
            result[i&mask,a,j&mask,b]+=value
    return {key:value for key,value in result.items() if value}


def fixed_instrument(family,density,edge,*,law='mixture',c_weight=F(1,2),repair=F(1)):
    """All25 literal fixed-seam slots on global binary masks."""
    if type(family) is not native.Family or type(density) is not dict:
        raise ValueError('a native graph adapter and sparse source/reference operator required')
    if edge not in family.seams or law not in ('C','E','mixture'):
        raise ValueError('an existing fixed native seam and law required')
    if any(type(x) is not F or not 0<=x<=1 for x in (c_weight,repair)):
        raise ValueError('exact source probabilities required')
    for key,value in density.items():
        if (type(key) is not tuple or len(key)!=4 or any(type(i) is not int or i<0 for i in key)
                or max(key[0],key[2]).bit_length()>len(family.weights)
                or type(value) not in (int,F)):
            raise ValueError('bounded binary masks and exact rational matrix entries required')
    mix=c_weight if law=='mixture' else F(law=='C')
    rows=[(source.Record('hold'),{k:(1-repair)*v for k,v in density.items() if (1-repair)*v})]
    for kind,weight in (('C',mix),('E',1-mix)):
        for coin in (0,1):
            decoding=source.circuit.record_values(kind,coin)
            for raw in range(1<<source.circuit.RECORD_BITS[kind]):
                label=decoding.get(raw)
                current=None if label is None else label if kind=='C' else (label-(2*coin-1)*(abs(label)%2))//2
                record=source.Record('repair','fixed',edge,kind,coin,raw,label,current)
                coefficient=repair*weight*family.coin(edge,coin)
                rows.append((record,native.apply_branch(family,density,edge,kind,coin,label,coefficient)))
    return tuple(rows)


def instrument(family,density,schedule,final,observed=(),*,reduce=False,law='mixture'):
    """Complete observed CQ channel, summing only explicitly unread outcomes.

    Keys are original (command index, native raw slot, full native Record)
    tuples, retaining impossible observed records as zero blocks.
    """
    if type(reduce) is not bool:
        raise ValueError('explicit reduction flag required')
    info=backward(schedule,final,observed)
    if any(e not in family.seams for e in schedule):
        raise ValueError('the schedule contains a nonnative seam')
    ports=len(family.weights);subset(final,ports)
    initial=partial_trace(density,info['initial_support'],ports) if reduce else density
    frames={():initial}
    for t,edge in enumerate(schedule):
        if reduce and t not in info['retained_commands']:
            continue
        following={}
        for history,rho in frames.items():
            for slot,(record,block) in enumerate(fixed_instrument(family,rho,edge,law=law)):
                key=history+((t,slot,record),) if t in observed else history
                output=following.setdefault(key,defaultdict(F))
                for index,coefficient in block.items():output[index]+=coefficient
        frames={key:{i:v for i,v in block.items() if v} for key,block in following.items()}
    return {key:partial_trace(rho,final,ports) for key,rho in frames.items()}


def source_masks(value,density):
    c=source.admit(value)['carrier']
    if c['kind']=='federation':return dict(density)
    if c!={'kind':'local','shell':'H66'}:
        raise ValueError('a binary native carrier required')
    masks=tuple(sum(n<<k for k,n in enumerate(x)) for x in source.local.basis())
    return {(masks[i],a,masks[j],b):v for (i,a,j,b),v in density.items()}


def observe_history(value,history,final,observed):
    """Actual archive/FIFO extraction; resource stops cannot become success."""
    audit_program(value)
    if type(history) is not source.History or history.program_digest!=source.digest(source.admit(value)):
        raise ValueError('the returned history must bind the exact admitted source program')
    if history.status!='COMPLETE':
        raise ValueError('a returned host stop or intrinsic timeout must be retained separately')
    if history.attempted_steps!=sum(command['repeat'] for command in value['commands']):
        raise ValueError('every admitted command must have completed before observation')
    if any(command['op']!='repair' or command['repeat']!=1 or
           command['selector']['kind']!='fixed' for command in value['commands']):
        raise ValueError('this extraction adapter requires literal unit fixed-seam commands')
    family=source_family(value);frames={}
    subset(final,len(family.weights))
    subset(observed,len(value['commands']))
    addresses={t:{outcome.record:slot for slot,outcome in enumerate(
        source._outcomes(value,value['commands'][t],{}))} for t in observed}
    for branch in history.branches:
        records=branch.archive+branch.memory
        if len(records)!=history.attempted_steps or branch.halted:
            raise ValueError('complete actual record histories required')
        cut=max(0,len(records)-value['memory']['capacity'])
        if (branch.archive,branch.memory)!=(records[:cut],records[cut:]):
            raise ValueError('the actual archive/FIFO partition cannot be relabeled')
        key=[]
        for t in observed:
            record=records[t]
            if record not in addresses[t]:raise ValueError('an observed record has no original address')
            key.append((t,addresses[t][record],record))
        output=frames.setdefault(tuple(key),defaultdict(F))
        reduced=partial_trace(source_masks(value,branch.density),final,len(family.weights))
        for index,coefficient in reduced.items():output[index]+=coefficient
    return {key:{i:v for i,v in block.items() if v} for key,block in frames.items()}
