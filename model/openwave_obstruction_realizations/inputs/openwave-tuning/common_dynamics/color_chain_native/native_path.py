"""Finite lazy M6 programs for the colour-chain observation experiment.

The frozen CAR bridge supplies the exact rational M6 compiler and pointer
services only. The SU(3) Hamiltonian and Gauss embedding are this amendment.
"""
from fractions import Fraction as F
from pathlib import Path
from functools import lru_cache
import importlib.util,sys
import source

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
PROVIDER=ROOT/'openwave-tuning/next/CHAIN5/species_native_bridge'
spec=importlib.util.spec_from_file_location('color_frozen_car_bridge',PROVIDER/'bridge.py')
car=importlib.util.module_from_spec(spec);spec.loader.exec_module(car)
native=car.native;fw=car.fw

def remap(word,wires):
    for k,t,a in word:yield k,tuple(wires[j] for j in t),a

def pair_word(left,right,angle,wires,controls=()):
    """exp(i angle (|left><right|+h.c.)); angle is exact (pi,alpha,rational)."""
    b=len(wires);controls=tuple(controls);diff=left^right
    assert 0<=left<right<1<<b
    for mask in range(1<<b):
        ny=(mask&diff).bit_count()
        if ny%2:continue
        sign=(-1)**(ny//2+(left&mask).bit_count())
        lab=''.join(('Y' if mask>>j&1 else 'X') if diff>>j&1 else ('Z' if mask>>j&1 else 'I') for j in range(b))
        for cm in range(1<<len(controls)):
            cs=(-1)**sum(v for j,(q,v) in enumerate(controls) if cm>>j&1)
            label=lab+''.join('Z' if cm>>j&1 else 'I' for j in range(len(controls)))
            coef=F(2*sign*cs,1<<(b+len(controls)))
            yield from remap(native.pauli(label,tuple(coef*x for x in angle)),tuple(wires)+tuple(q for q,v in controls))

def diagonal_word(values,time,wires,controls=()):
    """exp(-i time diag(values)); all padding columns are specified."""
    b=len(wires);assert len(values)==1<<b
    for mask in range(1<<b):
        c=-sum((v*(-1)**((j&mask).bit_count()) for j,v in enumerate(values)),F())/(1<<b)
        lab=''.join('Z' if mask>>j&1 else 'I' for j in range(b))
        for cm in range(1<<len(controls)):
            cs=(-1)**sum(v for j,(q,v) in enumerate(controls) if cm>>j&1)
            label=lab+''.join('Z' if cm>>j&1 else 'I' for j in range(len(controls)))
            yield from remap(native.pauli(label,tuple(c*cs*x/(1<<len(controls)) for x in time)),tuple(wires)+tuple(q for q,v in controls))

def evolution_word(d,time,subdivisions,precision,wires=None,control=None,field=F()):
    assert isinstance(subdivisions,int) and subdivisions>=1
    wires=tuple(range(d['bits'])) if wires is None else tuple(wires)
    controls=() if control is None else ((control,1),)
    rounded=source.rounding(d,precision,field);dt=tuple(x/subdivisions for x in time)
    for step in range(subdivisions):
        yield from diagonal_word(rounded['diagonal'],dt,wires,controls)
        for left,right,c in rounded['pairs']:
            yield from pair_word(left,right,tuple(-c*x for x in dt),wires,controls)

def had(bit):return remap(native.primitive_templates()['hadamard'],(bit,))

def swap(a,b):
    yield from native.cnot(a,b);yield from native.cnot(b,a);yield from native.cnot(a,b)

def qft_word(wires,inverse=False):
    word=[]
    for j in reversed(range(len(wires))):
        word.extend(had(wires[j]))
        for k in reversed(range(j)):
            word.extend(remap(native.diagonal([F(),F(),F(),F(1,1<<(j-k))]),(wires[k],wires[j])))
    for j in range(len(wires)//2):word.extend(swap(wires[j],wires[-1-j]))
    return native.negate_template(word) if inverse else word

def mixing_word(data,mixing):
    """Bell purification of I_D/D. Mixing registers remain as traced spectators."""
    for a,b in zip(data,mixing):
        yield from had(a);yield from native.cnot(a,b)

def qpe_unitary(d,phase,subdivisions,precision):
    for bit in phase:yield from had(bit)
    for j,bit in enumerate(phase):
        # Positive phase U=exp(+2 pi i H/32); evolution_word uses exp(-itH).
        yield from evolution_word(d,(F(-(1<<j),16),F(),F()),subdivisions[j],precision,control=bit)
    yield from qft_word(phase,True)

def layout(d,m=12):
    b=d['bits'];data=tuple(range(b));mix=tuple(range(b,2*b));phase=tuple(range(2*b,2*b+m))
    return dict(data=data,mixing=mix,phase=phase,clock=2*b+m,work=(2*b+m+1,2*b+m+2),memory_bits=2*b+m+3)

def pointer_word(bit,work):return fw.pointer_word(bit,*work)

def pulse_events(word,bits):
    for k,t,a in word:yield native.certified_pulse(k,t,a,bits)

def program(config,law='C',clock_axis='X',field=F(),time=F(),separated_label=None):
    """Typed finite controller, including all complements and first failures.

    Runtime executes literal reset/read instruments at the control events.
    It retains all phase words and stops only on an explicitly recorded branch.
    No large accepted program is claimed executed by providing this iterator.
    """
    if law not in ['C','E'] or clock_axis not in ['X','Y']:raise ValueError('fixed apparatus required')
    if not 0<=F(time)<=F(config['time_max_exact']) or abs(F(field))>F(config['field_max_exact']):
        raise ValueError('probe lies outside the certified observation window')
    import budget
    if config!=budget.configuration(config['L'],F(config['kappa']),config.get('weakup',1),F(config['time_max_exact']),F(config['field_max_exact'])):
        raise ValueError('the complete exact budget/configuration must be regenerated without edits')
    d=source.data(config['L'],F(config['kappa']),config.get('weakup',1));info=layout(d,config['phase_bits'])
    bits=config['precision_bits'];reset=config['reset_blocks_per_register']
    def resets(wires):
        for wire in wires:yield dict(instrument='literal_native_reset',wire=wire,blocks=reset,retain='all first-success fine tapes and timeout; stop only on recorded timeout')
    def read(bit,tag):
        yield from resets(info['work'])
        yield from pulse_events(pointer_word(bit,info['work']),bits)
        yield dict(instrument='literal_native_C_or_E',law=law,wire=bit,work=info['work'],tag=tag,
                   decode_function='decode_native_bit(law, coin, label)',
                   retain='all coin/current or coin/energy fine records and first export; exact decoded bit retained')
    yield from resets((*info['data'],*info['mixing'],info['clock']))
    if separated_label is None:
        yield from pulse_events(mixing_word(info['data'],info['mixing']),bits)
        for repetition in range(config['filter_repetitions']):
            yield from resets(info['phase'])
            yield from pulse_events(qpe_unitary(d,info['phase'],config['qpe_subdivisions'],bits),bits)
            for j,bit in enumerate(info['phase']):yield from read(bit,('phase',repetition,j))
        yield dict(instrument='exact_classical_filter',phase_bits=config['phase_bits'],cutoff_numerator=59,cutoff_denominator=512,
                   minimum_low_votes=(config['filter_repetitions']+1)//2,
                   repetitions=config['filter_repetitions'],retain='every full phase tape; rejected filter is a terminal preparation-failure branch, never renormalized')
    else:
        if not 0<=separated_label<d['dimension'] or d['span'][separated_label]!=d['L']-1:raise ValueError('actual maximal-span source codeword required')
        for j in range(d['bits']):
            if separated_label>>j&1:yield from pulse_events(remap(native.primitive_templates()['flip'],(j,)),bits)
    yield from pulse_events(had(info['clock']),bits)
    # External controlled preparation queries do not increment the subsequent
    # free-evolution clock. The declared source experiment starts here at tau=0.
    n=max(1,(abs(F(time))*F(config['subdivisions_per_unit_time'])).__ceil__())
    yield from pulse_events(evolution_word(d,(F(),F(),F(time)),n,bits,field=F(field)),bits)
    yield from pulse_events([('Z',(info['clock'],),native.angle(0,0,-F(time)/16))],bits)
    if clock_axis=='Y':yield from pulse_events([('Z',(info['clock'],),native.angle(F(1,4)))],bits)
    yield from pulse_events(had(info['clock']),bits)
    for bit in info['data']:yield from read(bit,('terminal_code',bit))
    yield from read(info['clock'],'clock')
    yield dict(instrument='exact_classical_decode',rule='valid code -> full orbital label and cell-occupation pattern; every padded label -> invalid-code outcome',
               spectators='mixing registers traced in the observed law; external quantum references retained',conditioning=False)

def filter_acceptance(phase_words,m=12,r=31):
    if len(phase_words)!=r or any(type(y)is not int or not 0<=y<1<<m for y in phase_words):
        raise ValueError('complete finite phase record required')
    return sum(512*y<59*(1<<m) for y in phase_words)>=(r+1)//2

def decode(d,label):
    if type(label)is not int or not 0<=label<d['code_dimension']:raise ValueError('complete code label required')
    return {'outcome':'occupation','orbitals':d['labels'][label],'pattern':d['patterns'][label]} if label<d['dimension'] else {'outcome':'invalid_code','label':label}

def decode_native_bit(law,coin,label):
    if coin not in (0,1) or law not in ('C','E') or label not in ((-1,0,1) if law=='C' else (-2,-1,0,1)):
        raise ValueError('complete proper-code native pointer record required')
    if law=='C':return int(2*label+2*coin-1<0)
    # The non-ready E=0 outcomes have zero ideal probability but are retained
    # under native errors; their decoder is defined without consulting J.
    return int(label<0) if label else 1-coin
