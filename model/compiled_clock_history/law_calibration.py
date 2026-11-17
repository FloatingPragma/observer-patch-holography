"""Propagate complete repair calibration into the compiled clock experiment.

The nine probabilities identify a common-current instrument on the 78-state
carrier. Restriction to its invariant 66-state binary sector and common
classical record coarse-graining are contractions. A stopped adaptive
history has an exact expected exposure to this instrument; its error need
not be multiplied by the worst-branch call count.

This is a conditional source prediction. Counting synthetic executions is
not physical calibration, and raw binary-array maps need not preserve trace.
"""
from fractions import Fraction as F
import hashlib
from pathlib import Path

from ..observer import quantum_mixing
from ..repair_identification import calibration, core
from ..record_laws import family
from ..quantum_observer_correspondence import construction as ideal
from . import history, receipt, source
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins


ROOT = Path(__file__).resolve().parents[2]


def nominal_table():
    """Read the executed C instrument without constructing unrelated controls."""
    return core.validate_table(family.extract_table(family.INSTRUMENTS['C']()))


def exposure(cap=2):
    """Expected number of native repairs, with every timeout retained.

    A capped incident wait has success q=1-(1-p)^h and mean q/p.
    After a first success there is one feedback repair and another wait.
    Thus mu=(1+q)q/p+q, independent of the state and admissible record law.
    This parameter counts native attempts, not source pulses or seconds.
    """
    history._arguments(0, cap, False)
    seams = quantum_mixing.seams()
    p = F(sum(bool(set(edge).intersection(source.UNION)) for edge in seams), len(seams))
    if p != ideal.P or p != F(2, 5):
        raise ArithmeticError('the compiled observer clock domain changed')
    q = 1-(1-p)**cap
    return {'incident_probability': p, 'first_success_probability': q,
            'complete_probability': q*q, 'timeout_probability': 1-q*q,
            'scheduled_attempts_mean': (1+q)*q/p,
            'feedback_attempts_mean': q, 'native_attempts_mean': (1+q)*q/p+q,
            'native_attempts_maximum': 2*cap+1}


def complete_bound(per_attempt_error, *, cap=2, gamma=1, reference=False):
    """Reference-uniform bound to any calibrated ideal instrument history.

    Telescope the chronological complete CP instruments. Before each call,
    its live branch has a state-independent scalar trace weight determined
    by the seam scheduler. Sum the per-call error times these weights.
    Identical CPTP intervening controls and record coarse-graining contract
    the trace norm. This includes the conditional feedback repair.

    The raw comparison is nominal compiled C versus the calibrated ideal
    history, by a triangle through ideal C. It does not execute the raw
    compiler with a different native table. Its raw
    binary maps are CP but are not asserted TP; their discrepancy cannot be
    capped at one, and an event weight can differ by twice the half norm.
    """
    history._arguments(gamma, cap, reference)
    error = calibration._probability(per_attempt_error)
    budget = exposure(cap)
    law = min(F(1), budget['native_attempts_mean']*error)
    # Emitted floating arrays can differ across supported platforms. Recheck
    # their exact residual, then retain the shared certified envelope rather
    # than that platform's residual as a supposedly portable exact constant.
    actual_compiled = receipt.bound(cap=cap, gamma=gamma, reference=reference)
    compiled = receipt.ERROR
    if not 0 <= actual_compiled <= compiled:
        raise ArithmeticError('emitted arrays exceed the retained compilation envelope')
    return {**budget, 'per_attempt_half_diamond_upper': error,
            'ideal_complete_half_trace_upper': law,
            'worst_branch_ideal_upper': min(F(1), (2*cap+1)*error),
            'raw_compilation_half_trace_upper': compiled,
            'raw_to_calibrated_ideal_half_trace_upper': law+compiled,
            'raw_event_weight_error_upper': law+2*compiled}


def from_tubes(centers, radii, *, cap=2, gamma=1, reference=False):
    """Bound every admissible table in fixed tubes relative to nominal C.

    The same homogeneous seam law is used at every relabelled seam; probing
    one seam does not establish that homogeneity. The nominal C instrument
    must itself satisfy the tubes. This supplies
    one feasible witness; it neither exhausts the feasible family nor
    chooses a law by projecting noisy estimates onto the CP cone.
    """
    history._arguments(gamma, cap, reference)
    nominal = nominal_table()
    tube = calibration.analyze_tubes(centers, radii, {'C': nominal}, 2*cap+1)
    if tube['compatible_locked_candidates'] != ['C']:
        raise ValueError('the nominal C instrument is outside the calibration tubes')
    probabilities = core.predictions(nominal)
    difference = F(0)
    for value, bounds in zip(probabilities, tube['probability_intervals'], strict=True):
        if value.b or value.c or value.d:
            raise ArithmeticError('the nominal C quadratures ceased to be rational')
        lower, upper = map(F, bounds)
        difference = max(difference, value.a-lower, upper-value.a)
    result = complete_bound(core.coordinate_error_bound(difference),
                            cap=cap, gamma=gamma, reference=reference)
    return {'calibration': tube, 'nominal_probability_difference_upper': difference,
            'prediction_budget': result, 'physical_law_selected': False,
            'scope': 'any CP table in these nine tubes, common current records, '
                     'same homogeneous relabelled seam law, identical ideal controls, '
                     'all timeouts, arbitrary untouched reference'}


def from_counts(evidence, design, service_errors, *, cap=2, gamma=1, reference=False):
    """Use the strict complete-attempt interface and retain its sampling caveats."""
    history._arguments(gamma, cap, reference)
    checked = calibration.analyze_counts(evidence, design, {'C': nominal_table()},
                                         2*cap+1, service_errors)
    result = from_tubes(tuple(map(F, checked['centers'])), tuple(map(F, checked['radii'])),
                        cap=cap, gamma=gamma, reference=reference)
    result['count_evidence'] = checked['count_evidence']
    return result


def event_interval(nominal_probability, ideal_complete_error):
    """Probability range for any fixed event of the complete normalized output."""
    p = calibration._probability(nominal_probability)
    error = calibration._probability(ideal_complete_error)
    return max(F(0), p-error), min(F(1), p+error)


def conditional_interval(nominal_joint, nominal_selected, ideal_complete_error):
    """A conditional event needs a positive lower bound on selected mass.

    No success-only renormalization of calibration or evidence is performed.
    The numerator and denominator can each move by the complete-law error;
    the resulting interval is conservative and retains the rare-event cost.
    """
    joint, selected, error = map(calibration._probability,
                                (nominal_joint, nominal_selected, ideal_complete_error))
    if joint > selected or selected <= error:
        raise ValueError('conditioning needs joint <= selected and selected > complete error')
    lower, upper = event_interval(joint, error)
    return lower/min(F(1), selected+error), min(F(1), upper/(selected-error))


def _json(value):
    if type(value) is F:
        return str(value)
    if type(value) is dict:
        return {key: _json(item) for key, item in value.items()}
    if type(value) in (tuple, list):
        return [_json(item) for item in value]
    return value


def run():
    centers = tuple(value.a for value in core.predictions(nominal_table()))
    radii = (F(1, 10000),)*9
    rows = []
    for gamma in (0, 1):
        calibrated = from_tubes(centers, radii, gamma=gamma)
        outputs = ideal.execute(gamma=gamma, cap=2)
        if sum(map(ideal.trace, outputs.values()), F(0)) != 1:
            raise ArithmeticError('the nominal complete ideal law lost probability')
        success = sum((ideal.trace(state) for (_, status, _), state in outputs.items()
                       if status == 'complete'), F(0))
        joint = sum((ideal.trace(state) for (_, status, bit), state in outputs.items()
                     if status == 'complete' and bit == 1), F(0))
        if success != exposure()['complete_probability']:
            raise ArithmeticError('the independent source execution changed its stopping law')
        error = calibrated['prediction_budget']['ideal_complete_half_trace_upper']
        rows.append({'gamma': gamma, 'forecast_status': 'conditional_source_forecast',
                     'event': 'both incident waits complete and terminal detector bit is one',
                     'nominal_joint_probability': joint, 'nominal_complete_probability': success,
                     'joint_probability_interval': event_interval(joint, error),
                     'success_conditioned_interval': conditional_interval(joint, success, error),
                     **calibrated})
    paths = set(core.SOURCE_PATHS) | set(ideal.SOURCE_PATHS) | set(source.bindings()) | {
        'model/predictions/__init__.py', 'model/record_laws/__init__.py',
        'model/repair_identification/__init__.py',
        'model/repair_identification/calibration.py',
        'model/compiled_clock_history/law_calibration.py',
        'model/compiled_clock_history/history.py', 'model/compiled_clock_history/receipt.py'}
    contrast = rows[1]['nominal_joint_probability']-rows[0]['nominal_joint_probability']
    uncertainty = sum((row['prediction_budget']['ideal_complete_half_trace_upper']
                       for row in rows), F(0))
    return _json({'schema': 'oph.compiled_clock_history.law_calibration.v1',
                  'source_sha256': bound_pins({path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
                                    for path in sorted(paths)}),
                  'executions': rows, 'all_pass': True,
                  'intervention_contrast': {
                      'event': 'complete and detector one, gamma one minus gamma zero',
                      'nominal_probability_difference': contrast,
                      'calibrated_interval': (max(F(-1), contrast-uncertainty),
                                              min(F(1), contrast+uncertainty)),
                      'same_law_and_calibration_for_both_settings': True},
                  'calibration_provenance': 'declared probability tubes; no measured counts',
                  'physical_dictionary_validated': False,
                  'unconditional_continuum_or_lab_clock_claim': False})


if __name__ == '__main__':
    import argparse
    from tools import strict_json
    from ..observer.scalar_response.experiment import canonical, compare_receipt
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare_receipt(result, strict_json.loads(args.check.read_text(encoding='utf-8')))
        print('Complete calibrated repair-to-clock observation bounds reproduced.')
    else:
        print(canonical(result), end='')
