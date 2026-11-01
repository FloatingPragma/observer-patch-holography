"""Independent literal source and all-Fock field error certificates."""
from collections import Counter
from fractions import Fraction as F

from ..certified_reference import certificate as cert, continuation as templates, interval as iv
from ..cube_eta_protection import source as preparation
from ..full_fock_field.local_simulation import split_bound
from ..common_source.chiral_candidate import gaussian_source as gaussian
from ..interacting_field import source
from . import circuits as c


TIME = F(1)
SEGMENTS = 2
STEPS = 32


def local(word, template):
    offset, angles, variables = cert.consume(word, 0, template)
    if offset != len(word) or variables: raise ArithmeticError('unconsumed literal source request')
    return angles+cert.controls_bound(word, {})


def source_certificate():
    tree, copy = preparation.tree_word(), preparation.copy_word()
    amplitude = cert.amplitude_certificate(tree, preparation.reference_values(), bits=8)
    if amplitude['offset'] != len(tree): raise ArithmeticError('incomplete source preparation tree')
    copy_template = [event for j in range(8) for event in cert.cnot(2*j, 2*j+1)]
    prep = amplitude['tree_error']+amplitude['fixed_angle_error']+cert.controls_bound(tree, {})+local(copy, copy_template)
    pointer = local(c.pointer_prepare(), c.pointer_template())
    analyzer = [local(c.pointer_analyzer(phase), c.pointer_template(phase)) for phase in (F(), F(1, 2))]
    duration = TIME/SEGMENTS; half = duration/(2*STEPS)
    template_table = {c.execute.key(gaussian.local_gate('fswap').word): templates.swap()}
    for group in c.GEOMETRY.groups:
        for _, _, z in group:
            template_table[c.execute.key(gaussian.local_gate('hopping', z, float(half)).word)] = templates.hopping(z, half)
    _, contact = source.contact(float(half), c.GEOMETRY.interaction)
    template_table[c.execute.key(contact['word'])] = templates.contact(half)
    blocks = c.field_binding(duration, STEPS)
    counts = Counter(instrument.tape for instrument, _ in blocks)
    field_error = F()
    local_rows = []
    for tape, count in counts.items():
        word = [dict(kind=k, angle=a, targets=list(t)) for k, a, t in tape]
        error = local(word, template_table[tape])
        field_error += count*error
        local_rows.append({'word_sha256': cert.word_digest(word), 'pulses': len(word),
                           'blocks_per_field_segment': count, 'complete_local_source_error': str(error)})
    splitting = SEGMENTS*split_bound(c.GEOMETRY.sides, duration, STEPS, refined=True)
    source_only = [prep+pointer+SEGMENTS*field_error+error for error in analyzer]
    return {'preparation_complete_state_error': str(prep), 'pointer_preparation_error': str(pointer),
            'local_field_word_certificates': local_rows,
            'one_field_segment_source_error': str(field_error),
            'analyzer_source_errors': [str(x) for x in analyzer],
            'source_only_complete_state_errors': [str(x) for x in source_only],
            'full_Fock_Strang_error': str(splitting),
            'complete_native_state_errors': [str(splitting+x) for x in source_only],
            'probability_error_bounds': [str(splitting+x) for x in source_only],
            'field_duration': str(TIME), 'field_segments': SEGMENTS, 'steps_per_segment': STEPS,
            'field_error_bound_is_all_input_and_arbitrary_reference': True,
            'state_preparation_bound_requires_initialized_blank': True,
            'floating_block_residuals_used_as_theorem_premises': False,
            'physical_hardware_accuracy_established': False}


def characteristic_interval(duration):
    """Exact trig/residue enclosure of the known source spectral response."""
    from ..particle_identification import vacuum
    if type(duration) is not F or not 0 <= duration <= 2: raise ValueError('exact local comparison time in[0,2] required')
    real, imag = iv.Interval.point(F(1, 16)), iv.Interval.point(0)
    for row in vacuum.certificate()['spectral_rows']:
        lo, hi = map(F, row['omega_interval'])
        low, high = map(F, row['positive_normalized_probe_weight_interval'])
        # Interval angle variation is controlled by the global unit Lipschitz
        # constant of sine/cosine, independently of midpoint evaluations.
        midpoint, radius = duration*(lo+hi)/2, duration*(hi-lo)/2
        sine, cosine = iv.sine_cosine(midpoint)
        delta = iv.Interval.point(radius).hi
        sine = iv.Interval(sine.lo-delta, sine.hi+delta)
        cosine = iv.Interval(cosine.lo-delta, cosine.hi+delta)
        weight = iv.Interval(iv.Interval.point(low).lo, iv.Interval.point(high).hi)*F(15, 16)
        real += weight*cosine
        imag -= weight*sine
    return {'real': real.json(), 'imag': imag.json()}
