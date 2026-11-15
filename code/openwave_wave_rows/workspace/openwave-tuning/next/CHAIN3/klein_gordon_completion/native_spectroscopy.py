"""Native controlled-walk spectroscopy and complete finite resource bounds.

The controlled processor uses the same rational M6 controls, proper-code
services and native pointer as the frozen source. No controlled-unitary
service or terminal projector is supplied as an extra primitive.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import importlib.util
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.dont_write_bytecode = True
sys.path[:0] = [str(ROOT/'openwave-tuning/next/CHAIN3/larger_walk_connection'),
               str(ROOT/'openwave-tuning/next/CHAIN2/refinement_apparatus')]
import larger_walk as walk
import native_word
import apparatus as b
a = walk.adapter


def scale(angle, coefficient):
    return tuple(F(coefficient)*x for x in angle)


def mapped_pauli(label, angle, wires):
    return [(kind, tuple(wires[t] for t in targets), v)
            for kind, targets, v in a.pauli(label, angle)]


def conditional_pauli(label, angle, wires, control):
    """exp(i theta (1-Zc) P/2); every scalar factor is retained."""
    assert control not in wires
    return (mapped_pauli(label, scale(angle, F(1,2)), wires)
            + mapped_pauli(label+'Z', scale(angle, F(-1,2)), (*wires, control)))


def controlled_event(kind, targets, angle, control):
    if control in targets:
        raise ValueError('fresh control outside the source word required')
    if kind in 'XYZI':
        return conditional_pauli(kind, angle, targets, control)
    if kind == 'H':
        # H = exp(-i pi Y/8) Z exp(i pi Y/8). The two outer
        # factors cancel exactly when the control is zero.
        return ([('Y', targets, a.angle(F(1,8)))]
                + conditional_pauli('Z', angle, targets, control)
                + [('Y', targets, a.angle(F(-1,8)))])
    if kind == 'J':
        # The native proper-code J generator is diag(1,1,1,0).
        result = []
        for label, coefficient in (('II',F(3,4)),('ZI',F(1,4)),
                                   ('IZ',F(1,4)),('ZZ',F(-1,4))):
            result += conditional_pauli(label, scale(angle,coefficient), targets, control)
        return result
    raise ValueError('unknown native source primitive')


def controlled_word(word, control):
    for kind, targets, angle in word:
        yield from controlled_event(kind, tuple(targets), angle, control)


def preparation_template(size, mode):
    """Literal QFT preparation of one integer momentum and uniform port/event."""
    q, d = b.member(size)
    if len(mode) != 3 or any(type(k) is not int or not 0 <= k < size for k in mode):
        raise ValueError('three Fourier mode labels in the finite register required')
    word = b.remap(a.uniform_template(), [0,1,2,3])
    word += b.remap(a.uniform_template(), [4,5,6,7])
    for axis, k in enumerate(mode):
        wires = list(range(8+(2-axis)*q, 8+(3-axis)*q))
        for j, wire in enumerate(wires):
            if k >> j & 1:
                word += b.remap(a.primitive_templates()['flip'], [wire])
        word += b.fourier(wires)
    return word


def reset_plan(size, tolerance):
    memory = a.layout(size)['memory_bits']
    blank, ready = memory+2, 2
    base, zero, orient = b.reset_base()
    repetitions = 1
    while (blank+ready)*base**repetitions > F(tolerance):
        repetitions += 1
    cap = max(zero['native_attempts_max'], orient['native_attempts_max'])
    return {'blank_carriers':blank, 'ready_carriers':ready,
            'blocks_per_register':repetitions, 'base_failure_exact':str(base),
            'timeout_bound_exact':str((blank+ready)*base**repetitions),
            'native_attempts_max':(blank+ready)*repetitions*cap,
            'fine_reset_records_retained':True, 'postselection':False}


@lru_cache(None)
def controlled_primitive_templates():
    return {name: list(controlled_word(word, 3))
            for name, word in a.primitive_templates().items()}


def certificate(size, steps, mass=F(1), mode=(1,0,0), precision=None, bandwidth=F(7)):
    """Uniform family bound without constructing the large route graph.

    The ideal controlled walk includes the inactive reference branch. The
    same full-state continuum bound applies on its active branch.
    """
    q,d = b.member(size)
    h = F(1,size)
    if type(steps) is not int or steps < 0:
        raise ValueError('nonnegative integer number of completed macros required')
    walk.mass_source.parameters(h,mass)
    resource = a.resource_upper_bound(size)
    memory = a.layout(size)['memory_bits']
    control, clock = memory, memory+1
    prep = preparation_template(size,mode)+b.H(control)+b.H(clock)
    coin = list(controlled_word(native_word.coin_template(h,mass), control))
    primitives = controlled_primitive_templates()
    # A harmless common address renaming fixes bounds independent of size.
    gate_count = resource['reversible_gates_upper']
    bits = precision or max(64,(resource['native_pulses_upper']*size*size).bit_length()+48)
    while True:
        dp = a.template_bound(prep,bits)
        dw = a.template_bound(coin,bits)+gate_count*max(a.template_bound(v,bits) for v in primitives.values())
        analysis = (b.H(control)+b.H(clock)
                    + [('Z',(control,),a.angle(F(1,4))),('Z',(clock,),a.angle(F(1,4)))])
        tick = [('Z',(clock,),a.angle(0,0,-h/16))]
        dr = a.template_bound(analysis,bits)+2*a.template_bound(b.detector_template(),bits)+steps*a.template_bound(tick,bits)
        if precision or max(dp,dw,dr) <= h*h:
            break
        bits *= 2
    reset = reset_plan(size,h*h)
    native = dp+steps*dw+dr+F(reset['timeout_bound_exact'])
    continuum = (walk.mass_source.continuum_bound(h,mass,bandwidth,steps)
                 if h*bandwidth <= F(1,2) else None)
    limit = F(continuum['uncapped_vector_bound_exact']) if continuum else None
    source_pulses = len(coin)+gate_count*max(map(len,primitives.values()))
    all_pulses = len(prep)+steps*(source_pulses+1)+len(analysis)+2*len(b.detector_template())
    return {'size':size, 'steps':steps, 'mass_exact':str(mass), 'mode':list(mode),
            'h_exact':str(h), 'source_time_in_tau0_exact':str(steps*h),
            'precision_bits':bits, 'source_memory_bits':memory,
            'control_wire':control, 'clock_wire':clock, 'total_memory_bits':memory+6,
            'delta_preparation_exact':str(dp), 'delta_controlled_step_exact':str(dw),
            'delta_clock_and_readout_exact':str(dr), 'reset':reset,
            'native_vs_finite_complete_TV_upper_exact':str(min(F(1),native)),
            'native_vs_finite_uncapped_exact':str(native),
            'finite_vs_continuum_vector_upper_exact':str(limit) if continuum else None,
            'native_vs_continuum_complete_TV_upper_exact':str(min(F(1),native+limit)) if continuum else None,
            'native_vs_continuum_uncapped_exact':str(native+limit) if continuum else None,
            'squared_frequency_error_upper_exact':continuum['squared_frequency_error_upper_exact'] if continuum else None,
            'controlled_source_pulses_upper':source_pulses,
            'whole_experiment_pulses_upper':all_pulses, 'load_export_services_upper':2*all_pulses,
            'large_word_executed':False,
            'complete_failures':'All first M6 export complements, reset attempts/timeouts, native C/E labels and retained source poststates.',
            'services':'Declared coherent proper-code load/export; finite native reset and pointer. No ideal controlled-unitary or terminal projector service.'}


def native_program(size, steps, mass, mode=(1,0,0), quadrature='X', clock_quadrature='X', precision=96):
    """Finite lazy program; large cases are not represented as executed."""
    if quadrature not in ('X','Y') or clock_quadrature not in ('X','Y'):
        raise ValueError('X/Y readout settings required')
    memory = a.layout(size)['memory_bits']
    control,clock = memory,memory+1
    h = F(1,size)
    yield {'operation':'native_configuration_reset_batch','plan':reset_plan(size,h*h),
           'blank_wires':list(range(memory+2)), 'ready_pairs':[[memory+2,memory+3],[memory+4,memory+5]]}
    for kind,targets,angle in preparation_template(size,mode)+b.H(control)+b.H(clock):
        yield {'operation':'M6',**a.certified_pulse(kind,targets,angle,precision)}
    for _ in range(steps):
        for kind,targets,angle in controlled_word(native_word.template_stream(size,h,mass),control):
            yield {'operation':'M6',**a.certified_pulse(kind,targets,angle,precision)}
        yield {'operation':'M6',**a.clock_tick_pulse(size,h,precision,clock)}
    for wire,setting,aux in ((control,quadrature,[memory+2,memory+3]),(clock,clock_quadrature,[memory+4,memory+5])):
        analysis=([('Z',(wire,),a.angle(F(1,4)))] if setting=='Y' else [])+b.H(wire)
        for kind,targets,angle in analysis+b.remap(b.detector_template(),[*aux,wire]):
            yield {'operation':'M6',**a.certified_pulse(kind,targets,angle,precision)}
        yield {'operation':'native_C_or_E','seam':[0,1],'auxiliary_wires':aux,
               'retained':['coin','current','label','output_configuration','source_poststate'],
               'read_bit':wire}
