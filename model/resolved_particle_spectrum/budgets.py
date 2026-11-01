"""Exact full-history error and compact native-word resource certificates."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache

from ..protected_particle_readout import circuits as c, bounds as previous
from ..certified_reference import certificate as cert, continuation as templates
from ..common_source.chiral_candidate import gaussian_source as gaussian
from ..interacting_field import source
from .factorization import STEP


STEPS_PER_INTERVAL = 2**19
TIMES = tuple(39*j for j in range(8))


@lru_cache(maxsize=1)
def certificate():
    blocks = c.field_binding(STEP, 1)
    counts = Counter(instrument.tape for instrument, _ in blocks)
    mapping = {c.execute.key(gaussian.local_gate('fswap').word): templates.swap()}
    for group in c.GEOMETRY.groups:
        for _, _, z in group:
            mapping[c.execute.key(gaussian.local_gate('hopping', z, float(STEP/2)).word)] = templates.hopping(z, STEP/2)
    _, contact = source.contact(float(STEP/2), c.GEOMETRY.interaction)
    mapping[c.execute.key(contact['word'])] = templates.contact(STEP/2)
    error = F(); rows = []
    for tape, count in counts.items():
        word = [dict(kind=k, angle=a, targets=list(targets)) for k, a, targets in tape]
        bound = previous.local(word, mapping[tape]); error += count*bound
        rows.append({'word_sha256': cert.word_digest(word), 'repetitions_per_step': count,
                     'complete_local_source_error': str(bound)})
    tape = c.flatten(blocks)
    if len(tape) != 10560: raise ArithmeticError('original source step pulse count changed')
    slots = sum(6-(1 << len(event['targets'])) for event in tape)
    if slots != 38976: raise ArithmeticError('original source step abort count changed')
    old = previous.source_certificate()
    static = [F(old['preparation_complete_state_error'])+F(old['pointer_preparation_error'])+F(x)
              for x in old['analyzer_source_errors']]
    schedule = []
    for j, time in enumerate(TIMES):
        steps = j*STEPS_PER_INTERVAL
        splitting = F(629, 256)*F(time**3, steps**2) if steps else F()
        bounds = [steps*error+splitting+x for x in static]
        if max(bounds) >= F(1, 65536): raise ArithmeticError('native long-time systematic budget exhausted')
        schedule.append({'time': time, 'steps_per_half': steps//2, 'total_steps': steps,
                         'field_native_pulses': 10560*steps,
                         'all_native_pulses_per_setting': 10560*steps+5443,
                         'field_abort_slots': slots*steps,
                         'all_source_abort_slots_per_setting': slots*steps+18664,
                         'complete_repair_records': 130, 'final_pointer_outcomes': 6,
                         'full_Fock_Strang_error': str(splitting),
                         'complete_repeated_field_source_error': str(steps*error),
                         'native_systematic_probability_bounds': [str(x) for x in bounds]})
    return {'step_duration': str(STEP), 'source_step_sha256': cert.word_digest(tape),
            'source_step_native_pulses': len(tape), 'source_step_abort_slots': slots,
            'complete_source_step_error': str(error), 'local_source_templates': rows,
            'schedule': schedule,
            'native_history_encoding': 'initial source; repeat(step,steps_per_half); five complete C seams; repeat(step,steps_per_half); five complete E seams; analyzer; six pointer exports',
            'native_abort_record_address': '(half,step_index,pulse_index,complement_slot); first failure terminates, all zero slots retained',
            'all_long_native_pulses_individually_executed': False,
            'all_long_native_failure_weights_individually_evaluated': False,
            'full_history_norm_bounds_include_all_native_exports': True,
            'outside_repair_continuations_are_not_postselected_or_renormalized': True}
