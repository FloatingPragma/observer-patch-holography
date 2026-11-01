"""Native charge instrument and binary clock for the original full cube H."""
from fractions import Fraction as F
from functools import lru_cache
import math

from model.observer.source_scalar import compiler
from model.certified_reference import certificate as cert, continuation as tmpl
from model.protected_particle_readout import execute
from model.common_source.chiral_candidate import gaussian_source as gaussian
from model.cube_eta_protection.cube import GEOMETRY


FIELD_BITS=16
CLOCK_BITS=2
M=1<<CLOCK_BITS
TAU=F(1)
STEPS=(16,32)


def number(angle):
    return float(angle.value)*math.pi if type(angle) is cert.PiMultiple else float(angle)


def emit(template):
    return tuple(compiler.pulse(kind,number(angle),list(targets)) for kind,targets,angle in template)


def block(template,bits,targets):
    word=emit(template)
    return execute.compile_block(execute.key(word),bits),targets


def hadamard_template(bit=0):
    return [('H',(bit,),cert._pi(1,2)),('I',(bit,),cert._pi(-1,2))]


def charge_template():
    # Read the input occupation before the odd unitary. The two output
    # columns are exactly c_1^dagger and c_1, including Jordan-Wigner Z_0.
    return cert.cnot(1,2)+cert.pauli('ZXI',cert._pi(1,2))+[('I',(0,),cert._pi(-1,2))]


def charge_blocks():return (block(charge_template(),3,(0,1,16)),)


def argument(z):
    real,imag=F(z.real),F(z.imag)
    if real and imag:raise ValueError('the source links lie on coordinate axes')
    if real:return cert._pi(0 if real>0 else 1),abs(real)
    if imag:return cert._pi(1 if imag>0 else -1,2),abs(imag)
    return cert._pi(0),F()


def controlled_hopping_template(z,duration):
    angle,magnitude=argument(z)
    result=tmpl.phase(0,-angle)
    for label in ('XX','YY'):
        result+=cert.pauli(label+'I',-duration*magnitude/4)
        result+=cert.pauli(label+'Z',duration*magnitude/4)
    return result+tmpl.phase(0,angle)


def controlled_contact_template(duration):
    result=[]
    for mask in range(8):
        label=''.join('Z' if mask>>j&1 else 'I' for j in range(3))
        result+=cert.pauli(label,-duration*F(GEOMETRY.interaction)*(-1)**mask.bit_count()/8)
    return result


@lru_cache(maxsize=32)
def controlled_field_blocks(duration,steps,clock_bit):
    if type(duration) is not F or duration<=0 or type(steps) is not int or steps<1:
        raise ValueError('exact positive duration and positive steps required')
    if type(clock_bit) is not int or not 16<=clock_bit<24:raise ValueError('external clock bit required')
    half=duration/(2*steps)
    swap_word=gaussian.local_gate('fswap').word
    swap=execute.compile_block(execute.key(swap_word),2)
    terms=tuple(term for group in GEOMETRY.groups for term in group)
    hops={z:execute.compile_block(execute.key(emit(controlled_hopping_template(z,half))),3)
          for _,_,z in terms}
    contact=execute.compile_block(execute.key(emit(controlled_contact_template(half))),3)
    result=[]
    def hop(p,q,z):
        route=tuple(range(q-1,p,-1))
        result.extend((swap,(j,j+1)) for j in route)
        result.append((hops[z],(p,p+1,clock_bit)))
        result.extend((swap,(j,j+1)) for j in reversed(route))
    for _ in range(steps):
        for p,q,z in terms:hop(p,q,z)
        for _ in range(2):
            result.extend((contact,(2*x,2*x+1,clock_bit)) for x in range(GEOMETRY.volume))
        for p,q,z in reversed(terms):hop(p,q,z)
    return tuple(result)


def ground_phase_word(duration):
    # E_G=-sqrt(3)/4 is an independently proved source value, not a fitted
    # calibration. exp(+it E_G |1><1|) fixes the requested excitation origin.
    angle=-float(duration)*math.sqrt(3)/8
    return (compiler.pulse('I',angle,[0]),compiler.pulse('Z',-angle,[0]))


def controlled_phase_template(phase):
    return sum((cert.pauli(label,cert.PiMultiple(phase*sign/4))
                for label,sign in (('II',1),('ZI',-1),('IZ',-1),('ZZ',1))),[])


def fourier_blocks(bits=CLOCK_BITS):
    if type(bits) is not int or not 1<=bits<=8:raise ValueError('one through eight clock bits required')
    result=[]
    # Positive-sign Fourier transform, so e^{-ij*tau*E} peaks at +E.
    for high in reversed(range(bits)):
        result.append(block(hadamard_template(),1,(16+high,)))
        for low in reversed(range(high)):
            result.append(block(controlled_phase_template(F(1,1<<(high-low))),2,(16+low,16+high)))
    for low in range(bits//2):
        high=bits-1-low
        result.append(block(cert.cnot(0,1)+cert.cnot(1,0)+cert.cnot(0,1),2,(16+low,16+high)))
    return tuple(result)


def clock_blocks():
    result=[block(hadamard_template(),1,(16+j,)) for j in range(CLOCK_BITS)]
    for j,steps in enumerate(STEPS):
        duration=TAU*(1<<j)
        result.extend(controlled_field_blocks(duration,steps,16+j))
        word=ground_phase_word(duration)
        result.append((execute.compile_block(execute.key(word),1),(16+j,)))
    result.extend(fourier_blocks())
    return tuple(result)
