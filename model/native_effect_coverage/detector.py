"""A fixed E spanning-tree program reads the initial H66 configuration.

The decoder uses only retained records. It corrects each incoming contrast
for earlier recorded displacement; native measuring back-action is retained.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from itertools import product

from model.native_program import interpreter
from model.local_repair_family import local_bias as local, modular_clock as clock
from model.observer import native

DIMENSION=66
PORTS=12


@lru_cache(maxsize=1)
def tree():
    seen={0};queue=[0];edges=[]
    for a in queue:
        for edge in local.SEAMS:
            if a not in edge:
                continue
            b=edge[0] if edge[1]==a else edge[1]
            if b not in seen:
                seen.add(b);queue.append(b);edges.append(edge)
    if seen!=set(range(PORTS)) or len(edges)!=PORTS-1:
        raise ArithmeticError('native seams do not supply the declared connected tree')
    return tuple(edges)


def program(edges=None):
    edges=tree() if edges is None else tuple(edges)
    value=interpreter.program(commands=[
        {'op':'repair','selector':{'kind':'fixed','seam':list(edge)},'law':'E','repeat':1}
        for edge in edges])
    expected={'port_weights':[25,25]+[8]*10,'seam_weights':[1]*30,
              'c_weight':'1/2','repair_probability':'1','exploration_probability':'1/2'}
    if value['parameters']!=expected or value['memory']!={'capacity':1,'overflow':'export_oldest'}:
        raise ArithmeticError('the authoritative weighted native member changed')
    return value


@dataclass(frozen=True)
class Branch:
    record:interpreter.Record
    weight:F
    pairs:tuple[tuple[int,int],...]


@lru_cache(maxsize=30)
def kernel(edge):
    value=program((edge,));marker={(i,i,i,i):F(1) for i in range(DIMENSION)}
    rows=[]
    for outcome in interpreter._outcomes(value,value['commands'][0],marker):
        record=outcome.record
        probability=local.coin_probability(clock.parameters(),edge,record.coin) if record.law=='E' else F()
        pairs=[]
        for (j,i,k,l),coefficient in outcome.density.items():
            if j!=k or i!=l or coefficient!=probability:
                raise ArithmeticError('the literal native map lost its weighted partial injection')
            pairs.append((i,j))
        rows.append(Branch(record,probability,tuple(sorted(pairs))))
    if len(rows)!=25 or sum(bool(row.pairs) for row in rows)!=6:
        raise ArithmeticError('the fixed E command lost a complete native address')
    return tuple(rows)


def decode(records,edges=None):
    """Reconstruct initial/final occupations using no hidden source state."""
    edges=tree() if edges is None else tuple(edges)
    if type(records) is not tuple or len(records)!=len(edges) or len(edges)!=11:
        raise ValueError('the complete eleven-record tree experiment is required')
    shift=[0]*PORTS;difference=[]
    for edge,record in zip(edges,records,strict=True):
        if (type(record) is not interpreter.Record or record.kind!='repair'
                or record.selector!='fixed' or record.seam!=edge or record.law!='E'
                or type(record.coin) is not int or record.coin not in (0,1)
                or record not in tuple(row.record for row in kernel(edge) if row.pairs)):
            raise ValueError('a live literal fixed-E record is required')
        a,b=edge
        difference.append(record.label-shift[a]+shift[b])
        shift[a]-=record.current;shift[b]+=record.current
    offsets={0:0}
    pending=list(zip(edges,difference,strict=True))
    while pending:
        next_pending=[]
        for (a,b),delta in pending:
            if a in offsets:
                wanted=offsets[a]-delta
                if b in offsets and offsets[b]!=wanted:
                    raise ValueError('inconsistent recorded tree contrasts')
                offsets[b]=wanted
            elif b in offsets:
                offsets[a]=offsets[b]+delta
            else:
                next_pending.append(((a,b),delta))
        if len(next_pending)==len(pending):
            raise ValueError('a connected source tree is required')
        pending=next_pending
    baseline=F(2-sum(offsets.values()),PORTS)
    initial=tuple(baseline+offsets[i] for i in range(PORTS))
    final=tuple(x+j for x,j in zip(initial,shift,strict=True))
    if (any(v not in (0,1) for v in initial+final) or sum(initial)!=2 or sum(final)!=2):
        raise ValueError('the retained word is not supported on the binary total-two source')
    return tuple(map(int,initial)),tuple(map(int,final))


def path_mapping(slots,edges=None):
    """Exact complete raw-word map; zero words have an empty mapping.

    This defines all 25**n outcomes lazily, including impossible prefixes.
    It is the partial-injection composition of the interpreter kernel.
    """
    edges=tree() if edges is None else tuple(edges)
    if (type(slots) is not tuple or len(slots)!=len(edges)
            or any(type(k) is not int or not 0<=k<25 for k in slots)):
        raise ValueError('one original raw address per literal command is required')
    mapping={i:i for i in range(DIMENSION)};weight=F(1);records=[]
    for edge,slot in zip(edges,slots,strict=True):
        row=kernel(edge)[slot];step=dict(row.pairs);weight*=row.weight
        mapping={i:step[j] for i,j in mapping.items() if j in step}
        records.append(row.record)
    if not weight:
        mapping={}
    return weight,mapping,tuple(records)


def branch(slots):
    weight,mapping,records=path_mapping(slots)
    if len(mapping)>1:
        raise ArithmeticError('the complete tree word did not isolate its initial source column')
    if mapping:
        initial,final=decode(records)
        states=local.basis()
        if mapping!={states.index(initial):states.index(final)}:
            raise ArithmeticError('record-only decoding disagrees with actual branch composition')
    return weight,mapping,records


def apply_path(slots,density,edges=None):
    """Complex-linear full source/reference branch action, without normalization."""
    weight,mapping,_=path_mapping(slots,edges)
    return {(mapping[i],a,mapping[j],b):weight*coefficient
            for (i,a,j,b),coefficient in density.items() if i in mapping and j in mapping and coefficient}


def run(initial,coins):
    """Execute a concrete deterministic native run; only decode sees records."""
    if type(initial) is not int or not 0<=initial<DIMENSION:
        raise ValueError('an original H66 basis index is required')
    if type(coins) is not tuple or len(coins)!=11 or any(type(c) is not int or c not in (0,1) for c in coins):
        raise ValueError('eleven original placement bits are required')
    states=local.basis();state=states[initial];records=[];slots=[];weight=F(1)
    for edge,coin in zip(tree(),coins,strict=True):
        before=state
        state,(actual_coin,current)=native.native(state,edge,coin)
        contrast=before[edge[0]]-before[edge[1]]
        candidates=[(k,row) for k,row in enumerate(kernel(edge))
                    if row.pairs and row.record.coin==actual_coin and row.record.label==contrast
                    and row.record.current==current]
        if len(candidates)!=1:
            raise ArithmeticError('native endpoint record has no unique retained E address')
        slot,row=candidates[0];records.append(row.record);slots.append(slot);weight*=row.weight
    if decode(tuple(records))!=(states[initial],state):
        raise ArithmeticError('native execution failed record-only reconstruction')
    return {'initial':initial,'final':states.index(state),'slots':tuple(slots),
            'records':tuple(records),'weight':weight,
            'archive':tuple(records[:-1]),'FIFO':tuple(records[-1:])}


def certificate():
    value=program();states=local.basis();checked=0
    for edge in tree():
        command=program((edge,))
        literal=interpreter._outcomes(command,command['commands'][0],{})
        if tuple(row.record for row in kernel(edge))!=tuple(row.record for row in literal):
            raise ArithmeticError('the source certificate lost an original raw address')
        for i,state in enumerate(states):
            for coin in (0,1):
                output,(_,current)=native.native(state,edge,coin)
                contrast=state[edge[0]]-state[edge[1]]
                candidates=[row for row in kernel(edge) if (i,states.index(output)) in row.pairs
                            and row.record.coin==coin and row.record.label==contrast
                            and row.record.current==current]
                if len(candidates)!=1 or candidates[0].weight!=local.coin_probability(clock.parameters(),edge,coin):
                    raise ArithmeticError('source endpoint and interpreter branch disagree')
                checked+=1
    signatures={tuple(x[a]-x[b] for a,b in tree()) for x in states}
    if len(signatures)!=DIMENSION:
        raise ArithmeticError('tree contrasts do not separate all source configurations')
    # A column has exactly one outcome for each independent placement bit.
    for edge in tree():
        for i in range(DIMENSION):
            weight=sum(row.weight for row in kernel(edge) if any(x==i for x,_ in row.pairs))
            if weight!=1:
                raise ArithmeticError('a native source column is incomplete')
    return {'source_program':value,'program_digest':interpreter.digest(value),
            'tree_edges':[list(e) for e in tree()], 'native_attempts':11,
            'raw_slots_per_attempt':25,'complete_raw_words':25**11,
            'nonzero_rank_one_words':DIMENSION*2**11,
            'complete_basis_projectors':DIMENSION,
            'source_kernel_columns_checked':checked,
            'initial_contrast_formula':'delta_initial=delta_record-shift[a]+shift[b]',
            'source_shift_update':'shift[a]-=current; shift[b]+=current',
            'initial_offset_formula':'n_i=(2-sum(offsets))/12+offset_i',
            'fine_branch':'sqrt(product coin probabilities)*|terminal><initial|',
            'coarse_effect':'|initial><initial|',
            'raw_update_is_Lueders':False,'arbitrary_untouched_reference':True,
            'native_timeout_probability':'0','law':'actual weighted E',
            'FIFO_slots':1,'exported_native_records':10,'retained_native_records':1,
            'native_record_storage_bits_with_terminal_status':value['declared_resources']['record_storage_bits_upper'],
            'decoder_additional_classical_working_bits_upper':1024,
            'decoder_working_storage_description':'64 signed16-bit cells, excluding immutable source chart and retained record archive',
            'generic_expanded_interpreter_finishes_eleven_steps':False,
            'owning_evaluator':'lazy exact partial-injection factorization of all complete raw words'}


def continue_raw(slots,density):
    """Stream one complete raw word through the existing continuation kernel.

    This is an unnormalized instrument query, not a random-outcome sampler or
    the expanded execute() API. It retains zero words and the actual FIFO.
    """
    value=program();path_mapping(slots)  # Validate the full address domain.
    current=interpreter._validate_density(density,value['carrier'])
    memory=archive=()
    for command,slot in zip(value['commands'],slots,strict=True):
        choices=interpreter._outcomes(value,command,current,memory,archive)
        outcome=choices[slot]
        current,memory,archive=outcome.density,outcome.memory,outcome.archive
    return {'density':current,'memory':memory,'archive':archive,
            'attempts':11,'program_digest':interpreter.digest(value),
            'status':'COMPLETE_RAW_WORD_QUERY'}


@lru_cache(maxsize=1)
def population_transition():
    """Coarse initial/terminal joint channel; fine tapes remain in branch()."""
    import sympy as s
    result=s.eye(DIMENSION)
    for edge in tree():
        step=s.zeros(DIMENSION)
        for row in kernel(edge):
            for i,j in row.pairs:
                step[j,i]+=s.Rational(row.weight.numerator,row.weight.denominator)
        result=step*result
    if s.ones(1,DIMENSION)*result!=s.ones(1,DIMENSION) or any(x<0 for x in result):
        raise ArithmeticError('the complete endpoint marginal lost probability')
    return s.ImmutableMatrix(result)


def aggregated_output(density,mode='raw',target=None):
    """Exact coarse instrument on arbitrary source/reference matrix units.

    Keys retain inferred initial label and actual terminal label. Restored
    and prepared outputs are ideal targets for the separately compiled words.
    """
    if mode not in ('raw','restore','prepare'):
        raise ValueError('a declared detector output mode is required')
    if mode=='prepare' and (type(target) is not int or not 0<=target<DIMENSION):
        raise ValueError('a declared source basis preparation target is required')
    rho=interpreter._validate_density(density,program()['carrier'])
    transition=population_transition();output={}
    for i in range(DIMENSION):
        diagonal={(a,b):v for (x,a,y,b),v in rho.items() if x==y==i}
        for j in range(DIMENSION):
            probability=F(transition[j,i])
            destination=j if mode=='raw' else i if mode=='restore' else target
            output[i,j]={(destination,a,destination,b):probability*v
                         for (a,b),v in diagonal.items() if probability*v}
    return output
