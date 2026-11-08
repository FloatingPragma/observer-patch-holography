"""A complete local E word and an exact finite mixing certificate."""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from math import comb
import numpy as np

from model.native_sector_readout import readout
from model.native_program import interpreter as source

FAMILY = readout.native.local_family(source.local.Parameters(port_weights=(1,)*12))
SEAMS = FAMILY.seams
TREE = readout.tree(FAMILY)
PORTS = 12
FULL = (1 << PORTS)-1
DENOMINATOR = 1 << 60
MAX_BLOCKS = (source.MAX_NODES-len(TREE))//60


def integer(value, low=0, high=None):
    if type(value) is not int or value<low or (high is not None and value>high):
        raise ValueError('a literal integer in the stated domain required')
    return value


def edges(blocks):
    integer(blocks,0,MAX_BLOCKS)
    return TREE + SEAMS*(2*blocks)


def commands(blocks,cell=0):
    integer(cell)
    return [{'op':'repair','selector':{'kind':'fixed','seam':[12*cell+a,12*cell+b]},
             'law':'E','repeat':1} for a,b in edges(blocks)]


def program(*,level=0,total=2,cell=0,blocks=32):
    integer(level);integer(total);integer(cell)
    return source.program(biased=False,federation={'level':level,'total':total},
                          commands=commands(blocks,cell),capacity=1,overflow='export_oldest')


def budget(blocks):
    schedule=edges(blocks)
    return {'blocks':blocks,'tree_attempts':len(TREE),'mixing_attempts':60*blocks,
            'attempts':len(schedule),'AST_nodes':len(schedule),
            'half_diamond_error_after_archive_marginal':str(F(1,2**blocks)),
            'maximum_admitted_blocks':MAX_BLOCKS,'all_local_numbers_retained':True}


@lru_cache(maxsize=13)
def masks(total):
    integer(total,0,PORTS)
    return readout.native.sector_masks(PORTS,total)


@lru_cache(maxsize=13)
def _block_bytes(total):
    """Columns are distributions; all arithmetic is exact signed int64.

    At step j every entry and column sum is at most 2**j. The final exponent
    is 60, below the signed63-bit bound. No floating operation is performed.
    """
    states=masks(total);index={mask:i for i,mask in enumerate(states)}
    value=np.eye(len(states),dtype=np.int64)
    for a,b in SEAMS*2:
        swap=[index[mask^((1<<a)|(1<<b))] if ((mask>>a)^(mask>>b))&1 else i
              for i,mask in enumerate(states)]
        value=value+value[swap,:]
    if (np.any(value<0) or not np.all(value.sum(axis=0)==DENOMINATOR)
            or not np.all(value.sum(axis=1)==DENOMINATOR)):
        raise ArithmeticError('the exact native doubly stochastic block changed')
    return len(value),value.astype('<i8').tobytes(order='C')


def block_numerators(total):
    size,payload=_block_bytes(total)
    return np.frombuffer(payload,dtype='<i8').reshape((size,size))


def minorization():
    rows=[]
    for total in range(13):
        value=block_numerators(total)
        mass=sum(map(int,value.min(axis=1)))
        if 2*mass<DENOMINATOR:
            raise ArithmeticError('the two-sweep native minorization is below one half')
        rows.append({'total':total,'dimension':len(value),'denominator':DENOMINATOR,
                     'minorization_numerator':mass,'minorization_mass':str(F(mass,DENOMINATOR)),
                     'matrix_sha256':sha256(value.astype('<i8').tobytes(order='C')).hexdigest()})
    return rows


def validate_density(density):
    if type(density) is not dict:raise ValueError('an exact sparse source/reference operator required')
    for key,value in density.items():
        if (type(key) is not tuple or len(key)!=4 or any(type(i) is not int or i<0 for i in key)
                or max(key[0],key[2])>FULL or type(value) not in (int,F)):
            raise ValueError('binary local masks, reference indices and rational entries required')
    return density


def reference_channel(density):
    """All sectors, with the complete two-dimensional constant corner intact."""
    validate_density(density);result=defaultdict(F)
    for (x,a,y,b),value in density.items():
        if x in (0,FULL) and y in (0,FULL):result[x,a,y,b]+=value
        elif x==y:
            total=x.bit_count();coefficient=F(value,comb(PORTS,total))
            for z in masks(total):result[z,a,z,b]+=coefficient
    return {key:value for key,value in result.items() if value}


def unread(density,schedule):
    """Exact marginal of all 25 raw addresses per supplied local E command."""
    validate_density(density)
    if type(schedule) is not tuple or any(e not in SEAMS for e in schedule):
        raise ValueError('literal existing internal seams required')
    for edge in schedule:
        following=defaultdict(F)
        for _,block in readout.cone.fixed_instrument(FAMILY,density,edge,law='E'):
            for key,value in block.items():following[key]+=value
        density={key:value for key,value in following.items() if value}
    return density


def local_trace(mask,blocks,coins):
    """A recorded full word without reading a source-state attribute at runtime."""
    integer(mask,0,FULL);schedule=edges(blocks)
    if type(coins) is not tuple or len(coins)!=len(schedule) or any(type(c) is not int or c not in (0,1) for c in coins):
        raise ValueError('one literal native coin per command required')
    start=mask;slots=[];records=[]
    for edge,coin in zip(schedule,coins,strict=True):
        following,current,contrast=readout.native.image(mask,edge,coin,PORTS)
        matches=[(i,r) for i,(r,w) in enumerate(readout.alphabet(FAMILY,edge))
                 if w and (r.coin,r.label,r.current)==(coin,contrast,current)]
        if len(matches)!=1:raise ArithmeticError('literal raw E address is not unique')
        slot,record=matches[0];slots.append(slot);records.append(record);mask=following
    return {'initial_mask':start,'terminal_mask':mask,'slots':tuple(slots),
            'records':tuple(records),'weight':F(1,2**len(schedule))}


def source_query(*,level=0,total=2,cell=0,blocks=0,mask=3,outside=0,coins=None):
    """Authenticate a lazy unnormalized raw word against the original kernel."""
    if level>source.MAX_EVALUATED_LEVEL:
        raise ValueError('literal kernel controls use the existing evaluated levels')
    value=program(level=level,total=total,cell=cell,blocks=blocks)
    if coins is None:coins=tuple(i%2 for i in range(len(edges(blocks))))
    trace=local_trace(mask,blocks,coins);shift=12*cell
    integer(outside)
    if outside & (FULL<<shift):raise ValueError('outside mask overlaps the local cell')
    start=(mask<<shift)|outside;finish=(trace['terminal_mask']<<shift)|outside
    density=source._validate_density({(start,0,start,0):F(1)},value['carrier'])
    memory=archive=()
    for command,slot in zip(value['commands'],trace['slots'],strict=True):
        outputs=source._outcomes(value,command,density,memory,archive)
        if len(outputs)!=25:raise ArithmeticError('the complete raw alphabet changed')
        chosen=outputs[slot];density,memory,archive=chosen.density,chosen.memory,chosen.archive
    if density!={(finish,0,finish,0):trace['weight']}:
        raise ArithmeticError('local word disagrees with the admitted federation kernel')
    if len(archive)!=len(edges(blocks))-1 or len(memory)!=1:
        raise ArithmeticError('complete FIFO/archive accounting changed')
    return {'program_digest':source.digest(value),'local_initial_mask':mask,
            'local_terminal_mask':trace['terminal_mask'],'attempts':len(edges(blocks)),
            'raw_weight':str(trace['weight']),'exported_records':len(archive),
            'retained_records':len(memory),'all25_addresses_per_attempt':True}


def archive_separation(total):
    integer(total,1,11)
    return F(comb(12,total)-1,comb(12,total))


def phase_blindness_control():
    """Actual local raw records do not access occupation-basis phases."""
    plus={(1,0,1,0):F(1,2),(2,0,2,0):F(1,2),(1,0,2,0):F(1,2),(2,0,1,0):F(1,2)}
    minus={key:(-value if key[0]!=key[2] else value) for key,value in plus.items()}
    checked=0
    for edge in SEAMS:
        for law in ('C','E'):
            left=readout.cone.fixed_instrument(FAMILY,plus,edge,law=law)
            right=readout.cone.fixed_instrument(FAMILY,minus,edge,law=law)
            for (record,a),(other,b) in zip(left,right,strict=True):
                if record!=other:raise ArithmeticError('raw labels differ')
                tr=lambda x:sum(value for (i,r,j,s),value in x.items() if i==j and r==s)
                if tr(a)!=tr(b):raise ArithmeticError('a native raw effect reads a diagonal phase')
                checked+=1
    return {'actual_raw_effect_controls':checked,'orthogonal_phase_inputs':True,
            'record_only_independent_service_response_identical':True,
            'proper_code_growing_source_load_derived':False}


def retained_tree_reference(slots,density):
    """One retained tree-record block of the limiting preparation instrument."""
    return reference_channel(readout.apply_branch(FAMILY,slots,validate_density(density)))
