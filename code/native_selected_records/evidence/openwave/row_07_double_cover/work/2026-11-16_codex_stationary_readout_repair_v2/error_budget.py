"""Explicit finite service-error and sign-assay bounds; no fitted thresholds."""
import math
from fractions import Fraction as F

SERVICE_CAP=2**22
ACCEPTANCE_FLOOR=F(3,100)


def budget(k,turns,calibrations=2):
    # A service is one specified native/inverse access, coherent memory gate,
    # or complete record/controller operation, not an unspecified hardware gate.
    q=256*2**k;frames=96*turns;filters=192*turns
    terms={'calibration_services_max':9*37568*calibrations,
        'preparation_services_max':6000,
        'rotation_services_max':2*filters*(9*21+8),
        'marker_services_max':9*frames*16,
        'two_counters_and_one_largest_ramsey_probe':10*(2*(q-1)+128*2**k),
        'counter_choice_qft_analyzer_services_max':2048}
    total=sum(terms.values())
    if total>SERVICE_CAP:raise RuntimeError('Declared service envelope exceeded')
    return {'terms':terms,'total_upper':total,'uniform_upper':SERVICE_CAP,
        'scope':'One complete trial and one chosen Ramsey axis/time, conditional on admitted entry; all failed prefixes included in worst-case work.'}


def envelope(bits=40,copies=1):
    delta=F(1,2**bits)
    epsilon=(2*copies*SERVICE_CAP+1)*delta
    p=ACCEPTANCE_FLOOR**copies
    if epsilon>=p:raise RuntimeError('Error envelope cannot preserve acceptance')
    return {'per_service_operator_error_upper':str(delta),'initial_trace_error_upper':str(delta),
        'copies':copies,'unconditional_trace_error_upper':str(epsilon),
        'ideal_acceptance_lower':str(p),'implemented_acceptance_lower':str(p-epsilon),
        'conditional_binary_probability_error_upper':str(epsilon/p),
        'conditional_visibility_error_upper':str(2*epsilon/p)}


def sign_assay():
    # Both-port interference conditions on the pre-analyzer herald, whose
    # probability exceeds the selected-clock lower bound. No sign-selective
    # conditioning is used to estimate visibility.
    error=float(F(envelope()['conditional_visibility_error_upper']))
    contrast=.95-error
    accepted=64;attempts=8192;conditions=224
    sign_failure=conditions*math.exp(-accepted*contrast**2/2)
    p=float(F(envelope()['implemented_acceptance_lower']))
    mean=attempts*p
    shortage=conditions*math.exp(-(mean-accepted)**2/(2*mean))
    if sign_failure+shortage>=1e-6:raise RuntimeError('Finite sign assay bound failed')
    return {'accepted_records_per_condition':accepted,'raw_attempt_cap_per_condition':attempts,
        'conditions_in_union':conditions,'wrong_sign_probability_upper':sign_failure,
        'maximum_complete_batch_service_calls':conditions*attempts*SERVICE_CAP,
        'insufficient_heralds_probability_upper':shortage,'union_upper':sign_failure+shortage,
        'meaning':'Prospective finite sign assay of the model, not sampled laboratory data or a high-precision energy measurement.'}
