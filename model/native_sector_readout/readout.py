"""Complete tree records on a conserved-number native binary source.

The total is public source specification, not a hidden state query. Sparse
raw-word queries retain all records without enumerating the global carrier.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import prod
from dataclasses import asdict

from model.native_program import interpreter as source
from model.native_local_response import cone
from model.weighted_source_family import family as native


@lru_cache(maxsize=16)
def tree(family):
    if type(family) is not native.Family:
        raise ValueError('an authenticated finite native graph adapter required')
    adjacency=[[] for _ in family.weights]
    for a,b in family.seams:
        adjacency[a].append(b);adjacency[b].append(a)
    seen={0};queue=[0];edges=[]
    for a in queue:
        for b in sorted(adjacency[a]):
            if b not in seen:
                seen.add(b);queue.append(b);edges.append(tuple(sorted((a,b))))
    if len(seen)!=len(family.weights):raise ValueError('connected binary graph required')
    return tuple(edges)


def program(*,level=None,total=None):
    if level is None:
        if total is not None and not (type(total) is int and total==2):
            raise ValueError('the actual local member has total two')
        options={}
    else:
        if type(level) is not int or not 0<=level<=source.MAX_EVALUATED_LEVEL:
            raise ValueError('only the evaluated native federation geometry is attached')
        if type(total) is not int or not 0<=total<=240*4**level:
            raise ValueError('read a valid conserved total from the source specification')
        options={'federation':{'level':level,'total':total},'biased':False}
    seed=source.program(**options)
    edges=tree(cone.source_family(seed))
    commands=[{'op':'repair','selector':{'kind':'fixed','seam':list(edge)},
               'law':'E','repeat':1} for edge in edges]
    return source.program(commands=commands,capacity=1,overflow='export_oldest',**options)


def specification(value):
    value=source.admit(value);family=cone.source_family(value)
    total=2 if value['carrier']['kind']=='local' else value['carrier']['total']
    expected=[{'op':'repair','selector':{'kind':'fixed','seam':list(edge)},
               'law':'E','repeat':1} for edge in tree(family)]
    if source.canonical(value['commands'])!=source.canonical(expected) or value['parameters']['repair_probability']!='1':
        raise ValueError('the complete fixed E tree with repair probability one is required')
    if value['memory']['overflow']!='export_oldest':
        raise ValueError('complete archive export is required; intrinsic halts are distinct')
    return value,family,total


@lru_cache(maxsize=16384)
def alphabet(family,edge):
    rows=cone.fixed_instrument(family,{},edge,law='E')
    if len(rows)!=25:raise ArithmeticError('the literal fixed raw alphabet changed')
    return tuple((record,family.coin(edge,record.coin) if record.law=='E' else F())
                 for record,_ in rows)


def _slots(slots,size):
    if type(slots) is not tuple or len(slots)!=size or any(
            type(k) is not int or not 0<=k<25 for k in slots):
        raise ValueError('one original raw slot for every tree command required')


def _same_record(left,right):
    return (type(left) is source.Record and type(right) is source.Record
            and all(type(getattr(left,k)) is type(getattr(right,k)) for k in vars(right))
            and (left.seam is None or all(type(x) is int for x in left.seam))
            and source.canonical(asdict(left))==source.canonical(asdict(right)))


def _same_records(left,right):
    return (type(left) is tuple and type(right) is tuple and len(left)==len(right)
            and all(_same_record(a,b) for a,b in zip(left,right,strict=True)))


def decode(family,records,total=None):
    """Return all supported (initial,terminal) masks, without enumeration.

    A known total gives zero or one pair. Without it the only two-column
    ambiguity is the all-empty/all-full subspace, whose coherence is retained.
    """
    edges=tree(family);ports=len(family.weights)
    if type(records) is not tuple or len(records)!=len(edges):
        raise ValueError('the complete original tree record word required')
    if total is not None and (type(total) is not int or not 0<=total<=ports):
        raise ValueError('a declared integer total or explicit unknown total required')
    shift=[0]*ports;offset={0:0}
    for edge,record in zip(edges,records,strict=True):
        if type(record) is not source.Record or not any(_same_record(record,r) for r,_ in alphabet(family,edge)):
            raise ValueError('the record is not an original address of this command')
        if record.law!='E' or record.label not in (-1,0,1):return ()
        a,b=edge;delta=record.label-shift[a]+shift[b]
        if a in offset:offset[b]=offset[a]-delta
        elif b in offset:offset[a]=offset[b]+delta
        else:raise ArithmeticError('tree traversal lost its known endpoint')
        shift[a]-=record.current;shift[b]+=record.current
    baselines=(F(total-sum(offset.values()),ports),) if total is not None else (F(),F(1))
    pairs=[]
    for baseline in baselines:
        initial=tuple(baseline+offset[i] for i in range(ports))
        terminal=tuple(initial[i]+shift[i] for i in range(ports))
        if any(x not in (0,1) for x in initial+terminal):continue
        # Check the entire native trajectory, not only binary endpoints.
        mask=sum(int(x)<<i for i,x in enumerate(initial));current=mask
        for edge,record in zip(edges,records,strict=True):
            after,j,d=native.image(current,edge,record.coin,ports)
            if (d,j)!=(record.label,record.current):break
            current=after
        else:
            target=sum(int(x)<<i for i,x in enumerate(terminal))
            if current!=target:raise ArithmeticError('record displacement lost an endpoint')
            pairs.append((mask,target))
    return tuple(pairs)


def branch(family,slots,total=None):
    edges=tree(family);_slots(slots,len(edges))
    if total is not None and (type(total) is not int or not 0<=total<=len(family.weights)):
        raise ValueError('a declared integer total or explicit unknown total required')
    rows=tuple(alphabet(family,e)[slot] for e,slot in zip(edges,slots,strict=True))
    records=tuple(row[0] for row in rows);weight=prod((row[1] for row in rows),start=F(1))
    mapping=decode(family,records,total) if weight else ()
    return weight,mapping,records


def apply_branch(family,slots,density,total=None):
    """One exact rational real component, without normalization.

    Arbitrary complex operators use separate real and imaginary calls and
    their complex-linear sum, as in the owning native interpreter.
    """
    ports=len(family.weights)
    if type(density) is not dict:raise ValueError('a sparse exact operator required')
    for key,value in density.items():
        if (type(key) is not tuple or len(key)!=4 or any(type(x) is not int or x<0 for x in key)
                or max(key[0],key[2]).bit_length()>ports or type(value) not in (int,F)):
            raise ValueError('binary masks, reference indices and exact rational entries required')
        if total is not None and (key[0].bit_count()!=total or key[2].bit_count()!=total):
            raise ValueError('input does not belong to the declared number sector')
    weight,pairs,_=branch(family,slots,total);mapping=dict(pairs)
    return {(mapping[i],a,mapping[j],b):weight*v for (i,a,j,b),v in density.items()
            if i in mapping and j in mapping and v and weight}


def trace_word(value,mask,coins):
    """Produce a concrete endpoint trace; the independent query authenticates it."""
    value,family,total=specification(value);edges=tree(family);ports=len(family.weights)
    if type(mask) is not int or mask<0 or mask.bit_length()>ports or mask.bit_count()!=total:
        raise ValueError('a source mask in the actual conserved sector required')
    if type(coins) is not tuple or len(coins)!=len(edges) or any(type(c) is not int or c not in (0,1) for c in coins):
        raise ValueError('one literal native placement bit per command required')
    current=mask;slots=[];records=[];weight=F(1)
    for edge,coin in zip(edges,coins,strict=True):
        following,j,d=native.image(current,edge,coin,ports)
        choices=[(i,r,w) for i,(r,w) in enumerate(alphabet(family,edge))
                 if w and r.coin==coin and r.label==d and r.current==j]
        if len(choices)!=1:raise ArithmeticError('endpoint execution lost its original E record')
        slot,record,w=choices[0];slots.append(slot);records.append(record);weight*=w;current=following
    if decode(family,tuple(records),total)!=((mask,current),):
        raise ArithmeticError('record-only reconstruction differs from endpoint execution')
    cut=max(0,len(records)-value['memory']['capacity'])
    return {'program_digest':source.digest(value),'initial_mask':mask,'terminal_mask':current,
            'slots':tuple(slots),'records':tuple(records),'weight':weight,
            'archive':tuple(records[:cut]),'memory':tuple(records[cut:])}


def source_density(value,density):
    """Convert only the original H66 mask chart; federation masks stay literal."""
    if value['carrier']['kind']=='federation':return density
    states=tuple(sum(n<<j for j,n in enumerate(row)) for row in source.local.basis())
    inverse={mask:i for i,mask in enumerate(states)}
    try:return {(inverse[i],a,inverse[j],b):v for (i,a,j,b),v in density.items()}
    except KeyError as e:raise ValueError('mask is outside the actual H66 chart') from e


def continue_raw(value,slots,density):
    """One unnormalized complete original-kernel query, retaining zero words."""
    value,family,total=specification(value);_slots(slots,len(tree(family)))
    expected=apply_branch(family,slots,density,total)
    actual=source._validate_density(source_density(value,density),value['carrier'])
    memory=archive=()
    for command,slot in zip(value['commands'],slots,strict=True):
        out=source._outcomes(value,command,actual,memory,archive)[slot]
        actual,memory,archive=out.density,out.memory,out.archive
    actual=cone.source_masks(value,actual)
    if actual!=expected:raise ArithmeticError('literal continuation and decoded full-sector branch differ')
    return {'program_digest':source.digest(value),'slots':slots,'density':actual,
            'memory':memory,'archive':archive,'attempts':len(slots),'status':'COMPLETE_RAW_WORD_QUERY'}


def authenticate_trace(value,trace):
    """Bind every field to the admitted program and independent original kernel."""
    if type(trace) is not dict or set(trace)!={'program_digest','initial_mask','terminal_mask','slots','records','weight','archive','memory'}:
        raise ValueError('the exact complete trace fields required')
    if trace['program_digest']!=source.digest(source.admit(value)):
        raise ValueError('the detector trace names another source program')
    value,family,total=specification(value)
    weight,pairs,records=branch(family,trace['slots'],total)
    if (type(trace['weight']) is not F or trace['weight']!=weight or not weight
            or type(trace['initial_mask']) is not int or type(trace['terminal_mask']) is not int
            or pairs!=((trace['initial_mask'],trace['terminal_mask']),)
            or not _same_records(trace['records'],records)):
        raise ValueError('trace state, weight or records differ from the complete word')
    result=continue_raw(value,trace['slots'],{(trace['initial_mask'],0,trace['initial_mask'],0):F(1)})
    if not _same_records(trace['archive'],result['archive']) or not _same_records(trace['memory'],result['memory']):
        raise ValueError('the actual FIFO/archive partition differs')
    return True


def reverse_currents(family,terminal,records):
    """Recover chronological masks from terminal readback and live C/E records.

    This checks native endpoint semantics. The caller must additionally bind
    the schedule, selector probabilities and FIFO to its actual source program.
    No initial configuration or auxiliary position is an input.
    """
    if type(family) is not native.Family:raise ValueError('an exact native graph adapter required')
    ports=len(family.weights)
    if type(terminal) is not int or terminal<0 or terminal.bit_length()>ports or type(records) is not tuple:
        raise ValueError('a binary terminal mask and literal complete record tuple required')
    current=terminal;reverse=[current]
    for record in reversed(records):
        if _same_record(record,source.Record('hold')):
            reverse.append(current);continue
        if (type(record) is not source.Record or record.kind!='repair'
                or record.selector not in ('fixed','uniform','feedback') or record.seam not in family.seams
                or record.law not in ('C','E') or type(record.coin) is not int or record.coin not in (0,1)
                or type(record.raw) is not int or type(record.label) is not int or type(record.current) is not int
                or any(type(x) is not int for x in record.seam)):
            raise ValueError('a live literal C/E record or hold required')
        if source.circuit.record_values(record.law,record.coin).get(record.raw)!=record.label:
            raise ValueError('the recorded raw code does not decode to its label')
        a,b=record.seam;j=record.current
        previous_a=((current>>a)&1)+j;previous_b=((current>>b)&1)-j
        if previous_a not in (0,1) or previous_b not in (0,1):
            raise ValueError('reverse current is outside the binary source')
        previous=(current&~((1<<a)|(1<<b)))|(previous_a<<a)|(previous_b<<b)
        target,actual_j,d=native.image(previous,record.seam,record.coin,ports)
        if target!=current or actual_j!=j or record.label!=(j if record.law=='C' else d):
            raise ValueError('reverse current fails the literal native endpoint update')
        current=previous;reverse.append(current)
    return tuple(reversed(reverse))


def no_reset_bound(value):
    """An unchanged-input branch obstructs exact finite repair-only reset."""
    value=source.admit(value);family=cone.source_family(value)
    attempts=sum(c['repeat'] for c in value['commands'])
    repairs=sum(c['repeat'] for c in value['commands'] if c['op']=='repair')
    minimum=min(family.coin(edge,c) for edge in family.seams for c in (0,1))
    repair=F(value['parameters']['repair_probability'])
    factor=1-repair+repair*minimum
    return {'attempts':attempts,'repair_attempts':repairs,'minimum_coin':str(minimum),
            'per_repair_no_motion_lower':str(factor),
            'no_motion_probability_lower_formula':f'({factor})^{repairs}',
            'no_motion_probability_lower':str(factor**repairs) if repairs<=512 else None,
            'exact_reset_on_nontrivial_sector_possible':False,
            'record_adaptive_native_selectors_covered':True,
            'coherent_services_or_postselection_covered':False}
