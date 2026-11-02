"""Literal native templates, including the full-domain signed field repair."""
from fractions import Fraction as F
from functools import lru_cache
import math
from ..observer.source_scalar import compiler
from ..certified_reference import certificate as cert
from ..common_source import native_routing
from ..defect_field_bridge import full_domain


def emit(template):
    for kind,targets,angle in template:
        value=float(angle.value)*math.pi if type(angle) is cert.PiMultiple else float(angle)
        yield compiler.pulse(kind,value,list(targets))


def remap(template,wires):return [(k,tuple(wires[j] for j in t),a) for k,t,a in template]
def xgate(bit):return [('I',(bit,),cert._pi(-1,2)),('X',(bit,),cert._pi(1,2))]


def encode(sites):
    if type(sites) is not int or sites<1:raise ValueError('positive number of spin sites required')
    return [event for s in range(sites) for event in xgate(2*s+1)+cert.cnot(2*s,2*s+1)]


def spin_rotation(axis,angle,site=0,sites=1):
    if axis not in 'XYZ' or type(site) is not int or not 0<=site<sites:raise ValueError('valid local spin axis/site required')
    label=['I']*(2*sites);label[2*site]=axis
    if axis!='Z':label[2*site+1]='X'
    return cert.pauli(''.join(label),angle)


def exchange(time,coupling=F(1),sites=2,seam=(0,1)):
    if type(time) is not F or type(coupling) is not F or coupling<=0:raise ValueError('exact time and positive exact coupling required')
    if type(sites) is not int or sites<2 or type(seam) is not tuple or len(seam)!=2 or not 0<=seam[0]<seam[1]<sites:
        raise ValueError('two valid increasing spin sites required')
    word=[];a=time*coupling/4
    for axis in 'XYZ':
        label=['I']*(2*sites)
        for s in seam:
            label[2*s]=axis
            if axis!='Z':label[2*s+1]='X'
        word+=cert.pauli(''.join(label),-a)
    word += [('I',(2*seam[0],),a)]
    return word


def probe(axis='Z',sites=2,site=0,pointer=None):
    if pointer is None:pointer=2*sites
    if axis not in 'XYZ' or not 0<=site<sites or pointer<2*sites:raise ValueError('spin axis and separate pointer required')
    basis=[]
    if axis=='X':basis=spin_rotation('Y',cert._pi(-1,4),site,sites)
    if axis=='Y':basis=spin_rotation('X',cert._pi(1,4),site,sites)
    return cert.inverse(basis)+cert.cnot(2*site,pointer)+basis


def fswap(left):
    # SWAP * CZ, including the scalar needed because native J(-pi)=-CZ.
    return cert.cnot(left,left+1)+cert.cnot(left+1,left)+cert.cnot(left,left+1)+[
        ('J',(left,left+1),cert._pi(-1)),('I',(left,),cert._pi(1))]


def route_pair(sites,seam):
    a,b=seam;labels=list(range(2*sites));word=[]
    for offset in range(2):
        old=labels.index(2*b+offset);target=2*a+2+offset
        while old>target:
            word+=fswap(old-1);labels[old-1],labels[old]=labels[old],labels[old-1];old-=1
    return word


@lru_cache(maxsize=2)
def repair_edges(law):
    if law not in ('C','E'):raise ValueError('C or E repair required')
    pairs=[]
    for coin in range(2):
        for mask in range(16):
            raw,out,sign,_=full_domain.branch_action(mask,law,coin,sites=2,flavors=2,seam=(0,1))
            if sign!=1:raise ArithmeticError('adjacent two-site local hopping sign changed')
            pairs.append((mask+(coin<<4),out+(coin<<4)+(raw<<5)))
    permutation=native_routing.complete_permutation(tuple(pairs))
    return native_routing.permutation_edges(permutation)


def controlled_flip(target,pattern):
    # Complete nine-bit reflection, with shared Gray-code parity phases.
    wires=[j for j in range(9) if j!=target]+[target]
    before=[('X',(j,),cert._pi(1,2)) for j in range(9) if pattern>>j&1]
    core=[('Y',(8,),cert._pi(1,4)),('I',(0,),cert._pi(-1,512))]
    for pivot in range(9):
        count=1<<pivot
        for j in range(count):
            core.append(('Z',(pivot,),cert._pi(1 if pivot==8 else -1,512)))
            x=(j^(j>>1))^(((j+1)%count)^(((j+1)%count)>>1))
            if x:core+=cert.cnot(x.bit_length()-1,pivot)
    core.append(('Y',(8,),cert._pi(-1,4)))
    return before+remap(core,wires)+cert.inverse(before)


def repair(law,sites=2,seam=(0,1)):
    if type(sites) is not int or sites<2 or type(seam) is not tuple or len(seam)!=2 or not 0<=seam[0]<seam[1]<sites:
        raise ValueError('valid ordered field seam required')
    a,b=seam;route=route_pair(sites,seam)
    wires=list(range(2*a,2*a+4))+list(range(2*sites,2*sites+5))
    core=[('Y',(4,),cert._pi(-1,4))]
    for target,pattern in repair_edges(law):core+=controlled_flip(target,pattern)
    return route+remap(core,wires)+cert.inverse(route)


def pulse_counts():
    return {'exchange_bond':len(exchange(F(1))),
            'C_repair_adjacent_sites':len(repair('C')),
            'E_repair_adjacent_sites':len(repair('E')),
            'fermionic_swap':len(fswap(0)),
            'extra_pulses_per_intervening_site':88,
            'single_Z_probe':len(probe('Z')),'single_X_probe':len(probe('X')),
            'single_Y_probe':len(probe('Y')),'encode_per_site':len(encode(1))}
