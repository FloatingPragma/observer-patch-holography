"""Record-only admission facade for the byte-frozen row1 integrity producer.

The frozen AU5 schema requires final_loads. Reconstruct it separately for
every hypothesized initial column and use only replay-compatible synthetic
states, then union the compatible columns. No supplied final_loads is read.
No prior over columns is assumed. This does not authenticate the E channel.
"""
from copy import deepcopy
from fractions import Fraction


def _finish(claim, G, region, retries, finals, *, producer):
    if not isinstance(claim, dict) or not isinstance(claim.get('meter'), dict):
        return {'admitted': False, 'reasons': ['I0: claim and pointer envelope required'],
                'compatible_pointer_initials': []}
    cleaned = deepcopy(claim)
    # Replace, do not validate or use, the untrusted diagnostic state.
    cleaned['meter'].pop('final_loads', None)
    compatible = set()
    observable = None
    rejections = set()
    for final_loads in finals:
        candidate = deepcopy(cleaned)
        candidate['meter']['final_loads'] = list(final_loads)
        verdict = producer.admit(candidate, G, region, retries)
        if not verdict['admitted']:
            rejections.update(verdict['reasons'])
            continue
        compatible.update(verdict['compatible_pointer_initials'])
        projection = {k: v for k, v in verdict.items() if k != 'compatible_pointer_initials'}
        if observable is not None and projection != observable:
            raise AssertionError('Compatible columns disagree on observable result')
        observable = projection
    if observable is None:
        return {'admitted': False, 'reasons': sorted(rejections),
                'compatible_pointer_initials': []}
    conditional_weight = Fraction(1, 2**observable['record_attempts'])
    return {**observable, 'compatible_pointer_initials': sorted(compatible),
            'compatible_initial_indicator': [int(x in compatible) for x in range(12)],
            'fine_record_weight_conditional_on_compatible_initial': str(conditional_weight),
            'pointer_initial_classical_effect': [str(conditional_weight if x in compatible else 0)
                                                 for x in range(12)],
            'probability_interpretation': 'Effect entries are complete emitted fine-record likelihoods conditional on each initial column; an explicit prior gives their weighted sum. No initial prior is assumed.'}


def admit_record_only_naive(claim, G, region, retries, *, producer):
    """Slow existential reference, retained for representative cross-checks."""
    return _finish(claim, G, region, retries,
                   [tuple(int(v == x) for v in range(12)) for x in range(12)], producer=producer)


def admit_record_only(claim, G, region, retries, *, producer):
    """Same existential union, with impossible final placeholders pruned."""
    try:
        if not isinstance(claim, dict) or not isinstance(claim.get('meter'), dict):
            raise ValueError('I0: claim and pointer envelope required')
        prepared = producer.preparation(G)
        ports = producer.region_set(region)
        latches = producer.observe_source(prepared, claim, claim.get('declared_total'))
        command = producer.lift(sum(latches['L1']['counters'][v] for v in ports) % 6, sum(prepared))
        meter = deepcopy(claim['meter'])
        meter.pop('final_loads', None)
        records = meter.get('records')
        if not isinstance(records, list):
            raise ValueError('I0: pointer tape required')
        for record in records:
            if (not isinstance(record, dict) or not isinstance(record.get('seam'), (list, tuple))
                or tuple(record['seam']) not in producer.EDGES
                or type(record.get('coin')) is not int or record['coin'] not in (0, 1)):
                raise ValueError('I0: malformed pointer event')
        finals = set()
        reasons = set()
        for initial in range(12):
            loads = tuple(int(v == initial) for v in range(12))
            for record in records:
                loads, _ = producer.strict.native(loads, tuple(record['seam']), record['coin'])
            candidate = deepcopy(meter)
            candidate['final_loads'] = list(loads)
            errors, _ = producer.strict.validate_pointer_tape(candidate, initial, command, retries)
            if errors:
                reasons.update(errors)
            else:
                finals.add(loads)
        if not finals:
            return {'admitted': False, 'reasons': ['I4: no record-only initial pointer column validates']+sorted(reasons),
                    'compatible_pointer_initials': []}
        return _finish(claim, G, region, retries, sorted(finals), producer=producer)
    except (ValueError, TypeError, KeyError, IndexError, OverflowError) as error:
        return {'admitted': False, 'reasons': [str(error)], 'compatible_pointer_initials': []}
