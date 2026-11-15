"""One complete massive field/position/clock experiment on frozen sources.

Source and apparatus words retain their native export/reset/commit records.
The clock occupies a distinct wire from every coarse-position flag.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
for directory in ('particles', 'vector_field', 'clock'):
    sys.path.insert(0, str(HERE.parent / directory))
import massive_coupled as massive
import particle_records as particle
import native_field
import phase_clock
from model.observer.source_scalar import compiler
from model.native_interface import pointer

CANDIDATE = 'native-local12-mass-fermion-plaquette-v3-proposed'
SOURCE_WIRES = 3203
FLAG_WIRES = tuple(range(3203, 3210))
CLOCK_WIRE = 3210
AUX_WIRES = (3211, 3212)


def admit(*, candidate=CANDIDATE, modulus=16, charge=1, steps=3,
          h=F(1, 8), mass=F(1), clock_wire=CLOCK_WIRE,
          flag_base=SOURCE_WIRES, auxiliaries=AUX_WIRES):
    if candidate != CANDIDATE or type(modulus) is not int or modulus != 16:
        raise ValueError('the frozen massive Z16 source is required')
    if type(charge) is not int or charge != 1:
        raise ValueError('this instrument is a single-particle massive-source join')
    if type(steps) is not int or not 0 <= steps <= 3:
        raise ValueError('zero through three complete macros required')
    if type(h) is not F or h != F(1, 8) or type(mass) is not F or mass not in (F(0), F(1)):
        raise ValueError('the frozen h=1/8, m=0 or1 development members are required')
    if type(flag_base) is not int or flag_base < SOURCE_WIRES:
        raise ValueError('flags must follow the complete source')
    flags = tuple(range(flag_base, flag_base + 7))
    if type(clock_wire) is not int or clock_wire < SOURCE_WIRES or clock_wire in flags:
        raise ValueError('clock aliases the source or position flags')
    if (type(auxiliaries) is not tuple or len(auxiliaries) != 2
            or any(type(w) is not int for w in auxiliaries)
            or len(set(auxiliaries)) != 2
            or any(w < SOURCE_WIRES or w in flags or w == clock_wire for w in auxiliaries)):
        raise ValueError('two fresh pointer-code wires required')
    return {'candidate': candidate, 'modulus': modulus, 'charge': charge,
            'steps': steps, 'h': h, 'mass': mass, 'clock_wire': clock_wire,
            'flags': flags, 'auxiliaries': auxiliaries}


def reset_operation(wire, target, phase):
    return {'operation': 'native_configuration_reset', 'phase': phase,
            'wire_or_auxiliary_code': wire, 'target_code': target,
            'sector': 'full78 total-two native configurations',
            'controller': 'configuration.policy(4,32,full=True,mode=reset)',
            'timeout': 'stop, retain fine tape, final carrier state and last completed macro'}


def pointer_operations(target, auxiliaries, law):
    if law not in ('C', 'E'):
        raise ValueError('native C or E law required')
    yield reset_operation(list(auxiliaries), pointer.READY_CODE, 'pointer_ready')
    for pulse in compiler.remap(pointer._local_word(), [*auxiliaries, target]):
        yield {'operation': 'pulse', 'phase': 'pointer', 'pulse': pulse}
    yield {'operation': 'native_commit', 'law': law, 'seam': list(pointer.SEAM),
           'auxiliary_code_wires': list(auxiliaries), 'work_target': target,
           'retained': 'coin, original C/E label, current, final native carrier',
           'decode': 'bit=int(2*current+2*coin-1<0)'}


def macro_word(h=F(1,8), mass=F(1), clock_wire=CLOCK_WIRE):
    """Chronological native mass kick, complete matter/field V, then tick."""
    admit(h=h, mass=mass)
    if type(clock_wire) is not int or clock_wire < SOURCE_WIRES:
        raise ValueError('fresh clock wire required')
    kick = massive.fock_kick_word(h, mass)
    for address in range(128):
        yield from compiler.remap(kick, [16*address+j for j in range(16)])
    yield from native_field.complete_word(.17, [.41] + [0.] * 207)
    yield from phase_clock.clock_tick_word(clock_wire, h, SOURCE_WIRES)


def complete_program(*, setting='X', law='C', **parameters):
    """A lazy full word/instrument plan; no terminal success normalization."""
    p = admit(**parameters)
    if setting not in ('X', 'Y') or law not in ('C', 'E'):
        raise ValueError('configured X/Y quadrature and native C/E law required')
    for wire in (*range(SOURCE_WIRES), *p['flags'], p['clock_wire']):
        yield reset_operation(wire, 0, 'blank')
    for pulse in massive.composition.preparation_word(1):
        yield {'operation': 'pulse', 'phase': 'source_preparation', 'pulse': pulse}
    local = phase_clock.local_words(0, setting, p['h'])
    for pulse in compiler.remap(local[:2], [p['clock_wire']]):
        yield {'operation': 'pulse', 'phase': 'clock_preparation', 'pulse': pulse}
    for step in range(p['steps']):
        for pulse in macro_word(p['h'], p['mass'], p['clock_wire']):
            yield {'operation': 'pulse', 'phase': 'mass_field_clock_macro',
                   'macro': step+1, 'pulse': pulse}
        yield {'operation': 'macro_completed', 'count': step+1,
               'declared_time_over_tau0': str((step+1)*p['h'])}
    for control, target in particle.coarse_position_calls(p['flags'][0]):
        for pulse in compiler.cnot(control, target):
            yield {'operation': 'pulse', 'phase': 'position_compute', 'pulse': pulse}
    for target in p['flags']:
        yield from pointer_operations(target, p['auxiliaries'], law)
    for control, target in particle.coarse_position_calls(p['flags'][0], undo=True):
        for pulse in compiler.cnot(control, target):
            yield {'operation': 'pulse', 'phase': 'position_uncompute', 'pulse': pulse}
    for target in range(2048, 3200):
        yield from pointer_operations(target, p['auxiliaries'], law)
    for pulse in compiler.remap(local[2:], [p['clock_wire']]):
        yield {'operation': 'pulse', 'phase': 'clock_analysis', 'pulse': pulse}
    yield from pointer_operations(p['clock_wire'], p['auxiliaries'], law)


def error_contract(*, setting='X', **parameters):
    p = admit(**parameters)
    base = massive.error_contract(p['h'], p['mass'], p['steps'])
    reset = particle.reset_bound(1, 1)
    clock = phase_clock.exact_word_bound(p['steps'], setting, p['h'])
    total = F(base['complete_source_preparation_evolution_and_coarse_field_readout_TV_exact']) + F(reset['timeout_upper_exact']) + clock
    return {'complete_TV_upper_exact': str(total),
            'complete_TV_upper': math.nextafter(float(total), math.inf),
            'source_and_position_field_bound_exact': base['complete_source_preparation_evolution_and_coarse_field_readout_TV_exact'],
            'additional_clock_native_word_and_pointer_error_exact': str(clock),
            'additional_clock_reset_bound_exact': reset['timeout_upper_exact'],
            'blank_native_resets': SOURCE_WIRES + 8,
            'ready_native_resets': 7 + 1152 + 1,
            'native_terminal_attempts': 7 + 1152 + 1,
            'source_clock_tensor_rule': 'V_m tensor exp(-i h Z/16), one tick after the complete mass/field macro',
            'reference_domain': 'arbitrary initialized source references; arbitrary coherent/correlated raw reset inputs through rank-one fine reset Kraus and conditional union bounds',
            'clock_preparation': 'declared separate plus state; source-correlated input is not a clock calibration',
            'comparison_boundary': 'private reset labels/tapes retained in the full instrument and marginalized only in comparison to the ideal source/clock law',
            'not_included': ['spectral controlled-V probe', 'physical seconds calibration',
                             'finite-shot phase confidence', 'two-particle position interpretation']}


def program_summary(**parameters):
    counts = {}; phases = {}; digest = hashlib.sha256(); max_wire = -1
    for operation in complete_program(**parameters):
        kind = operation['operation']; counts[kind] = counts.get(kind, 0) + 1
        if kind == 'pulse':
            phase = operation['phase']; phases[phase] = phases.get(phase, 0)+1
            max_wire = max(max_wire, *operation['pulse']['targets'])
        digest.update(json.dumps(operation, sort_keys=True, separators=(',', ':')).encode())
        digest.update(b'\n')
    return {'counts': counts, 'pulse_counts_by_phase': phases,
            'stream_sha256': digest.hexdigest(), 'max_wire': max_wire}


def block_summary(*, setting='X', law='C', **parameters):
    """Exact ordered repeated blocks; existing source words stay factored."""
    p=admit(**parameters)
    if setting not in ('X','Y') or law not in ('C','E'):
        raise ValueError('configured quadrature and native law required')
    mass_word=massive.fock_kick_word(p['h'],p['mass'])
    preparation=massive.composition.preparation_word(1)
    field_receipt=json.loads((HERE.parent/'vector_field/NATIVE_FIELD.json').read_text())
    composition=json.loads((HERE.parent/'vector_field/COMPOSITION.json').read_text())
    field_pulses=composition['source_pulses_per_matter_current_layer']+field_receipt['oneface_field_step_native_pulses']
    assert field_pulses==1209426
    local=phase_clock.local_words(0,setting,p['h'])
    calls=particle.coarse_position_calls(p['flags'][0])
    assert len(calls)==7168 and all(0<=a<2048 and b in p['flags'] for a,b in calls)
    mappings=[[16*address+j for j in range(16)] for address in range(128)]
    assert sorted(wire for group in mappings for wire in group)==list(range(2048))
    field_files={name:hashlib.sha256((HERE.parent/'vector_field'/name).read_bytes()).hexdigest()
                 for name in ('MANIFEST.json','FIELD_PROGRAM.json','NATIVE_FIELD.json','COMPOSITION.json','native_field.py')}
    program={
        'form':'ordered repeated blocks; complete_program exposes the full primitive expansion',
        'candidate':CANDIDATE,'h':str(p['h']),'mass':str(p['mass']),
        'operations':[
            {'op':'native_blank_resets','wires':[*range(SOURCE_WIRES),*p['flags'],p['clock_wire']],
             'target_code':0,'controller':'configuration policy(4,32,full78); retain every first failure'},
            {'op':'source_preparation','native_word':preparation},
            {'op':'clock_preparation','native_word':compiler.remap(local[:2],[p['clock_wire']])},
            {'op':'repeat_complete_macro','times':p['steps'],'ordered_body':[
                {'op':'remap_onsite_mass_word','native_word':mass_word,'ordered_wire_maps':mappings},
                {'op':'frozen_full_matter_field_word','source_pins':field_files,
                 'generator':'native_field.complete_word(.17,[.41]+[0.]*207,clock_spacing=None)',
                 'primitive_pulses':field_pulses},
                {'op':'native_clock_tick','native_word':phase_clock.clock_tick_word(p['clock_wire'],p['h'],SOURCE_WIRES)}],
             'on_completion':'record complete macro count and time=count*tau0*h'},
            {'op':'position_compute','ordered_CNOT_calls':calls,'local_native_CNOT':compiler.cnot(0,1)},
            {'op':'native_pointer_reads','targets':list(p['flags']),'auxiliaries':list(p['auxiliaries']),
             'local_native_word':pointer._local_word(),'native_law':law,'ready_reset_before_each':True},
            {'op':'position_uncompute','reverse_previous_CNOT_calls':True},
            {'op':'native_pointer_reads','targets':list(range(2048,3200)),
             'auxiliaries':list(p['auxiliaries']),'local_native_word':pointer._local_word(),
             'native_law':law,'ready_reset_before_each':True},
            {'op':'clock_analysis','native_word':compiler.remap(local[2:],[p['clock_wire']])},
            {'op':'native_pointer_reads','targets':[p['clock_wire']],
             'auxiliaries':list(p['auxiliaries']),'local_native_word':pointer._local_word(),
             'native_law':law,'ready_reset_before_each':True}],
        'source_failure_policy':'every pulse is its complete native loading/response/export instrument; retain complements and stop with last completed macro',
        'decoder':'all128 position-parity outcomes, all288 Z16 fields and fine native clock records; retain padding/failure labels; no success normalization'}
    phase_counts={'source_preparation':len(preparation),'clock_preparation':2,
                  'mass_field_clock_macro':p['steps']*(128*len(mass_word)+field_pulses+1),
                  'position_compute':3*len(calls),'position_uncompute':3*len(calls),
                  'pointer':6*1160,'clock_analysis':len(local[2:])}
    encoded=json.dumps(program,sort_keys=True,separators=(',',':')).encode()
    (HERE/'BLOCK_PROGRAM.json').write_bytes(encoded+b'\n')
    return {'program_form':program['form'],'counts':{'pulse':sum(phase_counts.values()),
            'native_configuration_reset':3211+1160,'native_commit':1160,'macro_completed':p['steps']},
            'pulse_counts_by_phase':phase_counts,'ordered_block_program_sha256':hashlib.sha256(encoded).hexdigest(),
            'max_wire':max(*p['auxiliaries'],p['clock_wire'],*p['flags'],SOURCE_WIRES-1),
            'repetition_proof':'Exact ordered concatenation; the same frozen mass/V/tick body repeats without an intervening hidden state reset. The block program preserves every scalar phase and current word.'}
