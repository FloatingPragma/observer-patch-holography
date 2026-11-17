"""Exact complete-source error and positive-ground spectral probability bounds."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from pathlib import Path

from model.certified_reference import certificate as cert, continuation as templates, interval as iv
from model.common_source.chiral_candidate import gaussian_source as gaussian
from model.interacting_field import source
from model.full_fock_field.local_simulation import split_bound
from model.cube_eta_protection.cube import GEOMETRY
from . import preparation as p, circuits as c


TIME = F(1)
STEPS = 32
SEGMENTS = 2


def local(word, template):
    offset, angle, variables = cert.consume(word, 0, template)
    if offset != len(word) or variables: raise ArithmeticError('unaccounted native control')
    return angle+cert.controls_bound(word, {})


@lru_cache(maxsize=1)
def certificate():
    preparation = F(p.certificate()['complete_preparation_state_error'])
    pointer = local(c.pointer_word(), c.pointer_template())
    analyzer = [local(c.pointer_word(phase), c.pointer_template(phase)) for phase in c.PHASES]
    duration = TIME/SEGMENTS; half = duration/(2*STEPS)
    plans = {c.execute.key(gaussian.local_gate('fswap').word): templates.swap()}
    for group in GEOMETRY.groups:
        for _, _, z in group:
            plans[c.execute.key(gaussian.local_gate('hopping', z, float(half)).word)] = templates.hopping(z, half)
    _, contact = source.contact(float(half), GEOMETRY.interaction)
    plans[c.execute.key(contact['word'])] = templates.contact(half)
    blocks = c.field_blocks(duration, STEPS)
    counts = Counter(instrument.tape for instrument, _ in blocks)
    field_error, rows = F(), []
    for tape, count in counts.items():
        word = [dict(kind=k, angle=a, targets=list(t)) for k, a, t in tape]
        error = local(word, plans[tape]); field_error += count*error
        rows.append({'native_word_sha256':cert.word_digest(word), 'pulses':len(word),
                     'blocks_per_segment':count, 'complete_operator_error':str(error)})
    splitting = SEGMENTS*split_bound(GEOMETRY.sides, duration, STEPS, refined=True)
    errors = [preparation+pointer+SEGMENTS*field_error+x+splitting for x in analyzer]
    if max(errors) >= F(3, 5000):
        raise ArithmeticError('declared complete Ramsey accuracy failed')
    return {'preparation_complete_state_error':str(preparation),
            'pointer_preparation_operator_error':str(pointer),
            'analyzer_operator_errors':[str(x) for x in analyzer],
            'field_segment_operator_error':str(field_error),
            'field_word_certificates':rows, 'full_Fock_Strang_error':str(splitting),
            'complete_native_state_and_probability_errors':[str(x) for x in errors],
            'field_duration':str(TIME), 'field_segments':SEGMENTS, 'steps_per_segment':STEPS,
            'source_pulse_counts':{'preparation':p.certificate()['tree_pulses']+p.certificate()['frame_pulses'],
                                 'pointer':len(c.pointer_word()),
                                 'field_per_segment':sum(len(i.tape) for i, _ in blocks),
                                 'analyzers':[len(c.pointer_word(phase)) for phase in c.PHASES]},
            'arbitrary_reference_field_operator_error':True,
            'initialized_blank_preparation_domain':True,
            'all_native_export_and_repair_outcomes_retained':True,
            'hardware_or_host_floating_accuracy_proved':False}


def spectral_rows():
    import model
    root = Path(model.__file__).resolve().parent.parent
    path = root/'model/retarded_matter_response/results.json'
    value = json.loads(path.read_text())
    if value.get('all_pass') is not True: raise ArithmeticError('unproved charged spectral source')
    for name, digest in value['source_bindings'].items():
        if hashlib.sha256((root/name).read_bytes()).hexdigest() != digest:
            raise ArithmeticError('charged spectral source changed: '+name)
    rows = value['resolved_local_spectrum']['rows']
    if len(rows) != 78 or value['resolved_local_spectrum']['unresolved_local_probe_weight'] != '0':
        raise ArithmeticError('complete local spectral measure required')
    return rows


def characteristic_interval(duration):
    if type(duration) is not F or not 0 <= duration <= 2:
        raise ValueError('exact comparison time in[0,2] required')
    real, imag = iv.Interval.point(0), iv.Interval.point(0)
    for row in spectral_rows():
        lo, hi = map(F, row['gap_above_ground_interval'])
        low, high = map(F, row['normalized_weight_interval'])
        midpoint, radius = duration*(lo+hi)/2, duration*(hi-lo)/2
        sine, cosine = iv.sine_cosine(midpoint)
        extra = iv.Interval.point(radius).hi
        sine = iv.Interval(sine.lo-extra, sine.hi+extra)
        cosine = iv.Interval(cosine.lo-extra, cosine.hi+extra)
        weight = iv.Interval(iv.Interval.point(low).lo, iv.Interval.point(high).hi)
        real += weight*cosine; imag -= weight*sine
    return {'real':real.json(), 'imag':imag.json()}


def probability_intervals():
    response = characteristic_interval(TIME)
    errors = list(map(F, certificate()['complete_native_state_and_probability_errors']))
    result = []
    for index, axis in enumerate(('real', 'imag')):
        lo, hi = map(F, response[axis])
        ideal = ((1+lo)/2, (1+hi)/2) if index == 0 else ((1-hi)/2, (1-lo)/2)
        actual = max(F(), ideal[0]-errors[index]), min(F(1), ideal[1]+errors[index])
        result.append({'phase_over_pi':str(c.PHASES[index]),
                       'ideal_pointer_zero_interval':[str(x) for x in ideal],
                       'literal_complete_pointer_zero_interval':[str(x) for x in actual]})
    if F(result[1]['literal_complete_pointer_zero_interval'][0]) <= F(1, 2):
        raise ArithmeticError('interacting ground quadrature not separated from free control')
    return result
