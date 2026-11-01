"""Source-executed faithful-reference clock with exact inverse enclosures.

The whole field Hamiltonian and this clock share lambda_j=1+j/4. The
measured time is accumulated action time, not source runtime or pulse count.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
import math

import numpy as np

from ..common_source.mass_response import history as instruments, model
from ..observer.source_scalar import compiler

REFERENCE_ENERGY = F(1, 4)  # GeV, one explicit calibration for every field
PHASE_LIMIT = F(3, 2)       # below pi/2; no phase alias chosen from a mass
SERVICE_ERROR = F(1, 10**8)


def _positive(value, name, zero=False):
    if isinstance(value, (bool, str)) or not isinstance(value, (int, float, F)):
        raise ValueError(f'{name} must be a finite real calibration')
    try:
        finite = math.isfinite(float(value))
    except (ValueError, OverflowError):
        finite = False
    if not finite:
        raise ValueError(f'{name} must be finite')
    result = F(value)
    if result < 0 or (not zero and result == 0):
        raise ValueError(f'{name} is outside its positive domain')
    return result


def drive(current):
    if type(current) is not int or current not in (-1, 0, 1):
        raise ValueError('one actual native current -1, 0 or 1 is required')
    return 1+F(current, 4)


def calibration(energy=REFERENCE_ENERGY):
    energy = _positive(energy, 'clock energy')
    reference = model.ideal_clock_calibration()
    populations = tuple(F(value) for value in reference['populations'])
    if populations != (F(1, 6), F(5, 6)):
        raise ArithmeticError('native faithful clock reference changed')
    return {'reference_energy_GeV': str(energy), 'populations': list(map(str, populations)),
            'contrast': str(populations[0]-populations[1]),
            'reference_modular_gap': 'log(5)',
            'inverse_temperature_GeV_inverse': 'log(5)/Omega_ref',
            'reference_identity': '-log(rho_ref)=constant I+(log(5)/Omega_ref) H_clock',
            'clock_Hamiltonian': 'H_clock=Omega_ref Z/2',
            'phase_window': ['0', str(PHASE_LIMIT)],
            'time_convention': 'tau_effective=sum_j lambda_j delta_tau_j; the same lambda_j multiplies the entire field Hamiltonian and both clock hands',
            'host_runtime_or_pulse_count_is_elapsed_time': False,
            'energy_calibration_is_species_dependent': False}


@lru_cache(maxsize=1)
def preparation():
    """Prepare two faithful clock hands once, retaining both purifiers."""
    p = tuple(F(value) for value in calibration()['populations'])
    amplitudes = np.sqrt(np.array(list(map(float, p))))
    word = []
    for clock, purifier in ((0, 2), (1, 3)):
        word += compiler.remap(instruments.preparation_word(amplitudes), [clock])
        word += compiler.cnot(clock, purifier)
        word += instruments._hadamard(clock)
    blank = np.zeros((16, 1), complex); blank[0, 0] = 1
    output, cost = instruments.execute_word(word, 4, blank)
    state = output[:, 0].reshape(4, 4).T.copy()  # clock work x purifier columns
    H = np.array([[1, 1], [1, -1]], complex)/math.sqrt(2)
    one = H@np.diag(amplitudes)
    ideal = np.kron(one, one)
    error = float(np.linalg.norm(state-ideal))
    if error > float(SERVICE_ERROR)/8:
        raise ArithmeticError('clock source preparation exceeds its numerical admission')
    return state, ideal, tuple(word), cost, error


def advance(state, current, duration, energy=REFERENCE_ENERGY):
    """One common-action interval; no reset and no change of reference."""
    duration, energy = _positive(duration, 'interval', zero=True), _positive(energy, 'clock energy')
    phase = energy*duration*drive(current)
    word = tuple(compiler.pulse('Z', -float(phase)/2, [bit]) for bit in (0, 1))
    output, cost = instruments.execute_word(word, 2, state)
    return output, word, cost


def read(state):
    """Two commuting X/Y hand reads, with all four output states retained."""
    first = instruments.pauli_read(state, 0, 'X')
    outputs, costs, words = {}, [], []
    for x, branch in enumerate(first.outputs):
        second = instruments.pauli_read(branch, 1, 'Y')
        costs.append(second.resources); words.append(second.word)
        for y, output in enumerate(second.outputs):
            outputs[x, y] = output
    if costs[0] != costs[1] or words[0] != words[1]:
        raise ArithmeticError('clock read cost unexpectedly depends on the first outcome')
    return outputs, (first.resources, costs[0]), (first.word, words[0])


def _resource_sum(rows):
    keys = set.union(*(set(row) for row in rows))
    result = {key: (max if key.startswith('maximum_') else sum)(row.get(key, 0) for row in rows) for key in sorted(keys)
            if all(type(row.get(key, 0)) in (int, float) for row in rows)
            and key not in ('work_qubits', 'reference_columns')}
    result['work_preparations'] = 1
    result['total_clock_qubits_including_purifiers_and_fresh_pointers'] = 6
    return result


def history(currents, duration, energy=REFERENCE_ENERGY):
    currents = tuple(currents)
    duration, energy = _positive(duration, 'interval', zero=True), _positive(energy, 'clock energy')
    effective = duration*sum((drive(j) for j in currents), F())
    if energy*effective > PHASE_LIMIT:
        raise ValueError('the declared clock phase window is exceeded; no alias is selected')
    initial, ideal_initial, prep_word, prep_cost, preparation_error = preparation()
    state, words, costs = initial.copy(), [], [prep_cost]
    for current in currents:
        state, word, cost = advance(state, current, duration, energy)
        words.extend(word); costs.append(cost)
    outputs, read_costs, read_words = read(state)
    costs.extend(read_costs)
    probabilities = {key: float(np.vdot(value, value).real) for key, value in outputs.items()}
    phase = float(energy*effective)
    c = float(F(calibration(energy)['contrast']))
    expected = {(x, y): (1+(-1)**x*c*math.cos(phase))*(1+(-1)**y*c*math.sin(phase))/4
                for x in (0, 1) for y in (0, 1)}
    if max(abs(probabilities[key]-expected[key]) for key in probabilities) > float(SERVICE_ERROR)/4:
        raise ArithmeticError('complete source clock disagrees with its calibrated circle')
    # The same source word acts on every work-input/reference column.
    U = compiler.execute(words, 2)
    complete, _, _ = read(U)
    phi = np.diag([np.exp(-.5j*phase), np.exp(.5j*phase)])
    ideal_U = np.kron(phi, phi)
    X = np.array([[0,1],[1,0]], complex); Y = np.array([[0,-1j],[1j,0]], complex)
    error = 0.
    for (x, y), K in complete.items():
        ideal_K = np.kron((np.eye(2)+(-1)**y*Y)/2,
                          (np.eye(2)+(-1)**x*X)/2)@ideal_U
        error += float(np.linalg.norm(K-ideal_K, 'fro')**2)
    error = math.sqrt(error)
    if error+preparation_error > float(SERVICE_ERROR):
        raise ArithmeticError('full clock instrument exceeds its source error admission')
    return {'outputs': outputs, 'probabilities': probabilities, 'calibration': calibration(energy),
            'work_kraus': complete, 'prepared_clock_and_purifiers': initial.copy(),
            'source_cost': _resource_sum(costs), 'effective_time_for_comparison_only': str(effective),
            'phase_for_comparison_only': str(energy*effective),
            'source_full_input_columns': 4, 'retained_purifier_columns': 4,
            'complete_instrument_frobenius_error': error, 'preparation_error': preparation_error,
            'numerical_service_error_bound': str(SERVICE_ERROR),
            'fresh_clock_hands': 2, 'fresh_purifiers': 2, 'new_clock_pointer_bits': 2,
            'preparations': 1, 'clock_resets': 0}


def _sin_bounds(x, terms=24):
    total = F()
    for k in range(terms):
        total += (-1)**k*x**(2*k+1)/math.factorial(2*k+1)
    next_term = (-1)**terms*x**(2*terms+1)/math.factorial(2*terms+1)
    return min(total, total+next_term), max(total, total+next_term)


def _cos_bounds(x, terms=24):
    total = F()
    for k in range(terms):
        total += (-1)**k*x**(2*k)/math.factorial(2*k)
    next_term = (-1)**terms*x**(2*terms)/math.factorial(2*terms)
    return min(total, total+next_term), max(total, total+next_term)


def _asin_bracket(value, steps=60):
    if value <= 0:
        return F(), F()
    lo, hi = F(), PHASE_LIMIT
    for _ in range(steps):
        middle = (lo+hi)/2
        left, right = _sin_bounds(middle)
        if right < value:
            lo = middle
        elif left > value:
            hi = middle
        else:
            break
    return lo, hi


def _acos_bracket(value, steps=60):
    if value >= 1:
        return F(), F()
    lo, hi = F(), PHASE_LIMIT
    for _ in range(steps):
        middle = (lo+hi)/2
        left, right = _cos_bounds(middle)
        if left > value:
            lo = middle
        elif right < value:
            hi = middle
        else:
            break
    return lo, hi


def infer_time(probabilities, *, energy=REFERENCE_ENERGY, probability_error=SERVICE_ERROR):
    """Record-only enclosure: receives no duration, history or matter energies.

    The input is the four retained clock outcome probabilities conditional
    on a known native record. probability_error bounds either binary clock
    marginal, including any explicitly propagated conditioning error.
    """
    energy = _positive(energy, 'clock energy')
    error = _positive(probability_error, 'probability error', zero=True)
    if set(probabilities) != {(0,0),(0,1),(1,0),(1,1)}:
        raise ValueError('all four clock outcomes are required')
    p = {key: _positive(value, 'outcome probability', zero=True) for key, value in probabilities.items()}
    if any(value > 1+error for value in p.values()) or abs(sum(p.values())-1) > 4*error:
        raise ValueError('clock probabilities do not have the declared mass')
    px, py = p[0,0]+p[0,1], p[0,0]+p[1,0]
    contrast = F(calibration(energy)['contrast'])
    sine = sorted(((2*(py-error)-1)/contrast, (2*(py+error)-1)/contrast))
    cosine = sorted(((2*(px-error)-1)/contrast, (2*(px+error)-1)/contrast))
    cap_sine = _sin_bounds(PHASE_LIMIT)
    if sine[1] < 0 or sine[0] > cap_sine[1]:
        raise ValueError('clock sine record lies outside the declared phase window')
    lo = _asin_bracket(max(F(), sine[0]))[0]
    hi = PHASE_LIMIT if sine[1] >= cap_sine[0] else _asin_bracket(sine[1])[1]
    c_lo = _cos_bounds(hi)[0]; c_hi = _cos_bounds(lo)[1]
    if cosine[1] < c_lo or cosine[0] > c_hi:
        raise ValueError('the X/Y clock records are incompatible with one calibrated hand phase')
    if lo > hi:
        raise ValueError('empty clock phase interval')
    return {'phase_interval': [str(lo), str(hi)],
            'effective_time_GeV_inverse_interval': [str(lo/energy), str(hi/energy)],
            'inverse_uses_matter_energies_or_program_duration': False,
            'inverse_method': 'exact rational alternating Taylor enclosures and 60 bisections on [0,3/2]',
            'probability_error_bound': str(error)}


def infer_energy_gap(ramsey_probability, clock_probabilities, *,
                     energy=REFERENCE_ENERGY, field_error=SERVICE_ERROR,
                     clock_error=SERVICE_ERROR):
    """One scale for every species: Delta=measured_field_phase/measured_time.

    The explicitly chosen window is0<=Delta*tau<=3/2. Neither the supplied
    field energies nor its programmed duration are arguments to this inverse.
    """
    p = _positive(ramsey_probability, 'Ramsey probability', zero=True)
    error = _positive(field_error, 'field probability error', zero=True)
    time = infer_time(clock_probabilities, energy=energy, probability_error=clock_error)
    tlo, thi = map(F, time['effective_time_GeV_inverse_interval'])
    if tlo <= 0:
        raise ValueError('the measured clock interval cannot resolve a positive duration')
    coslo, coshi = 2*(p-error)-1, 2*(p+error)-1
    caplo, caphi = _cos_bounds(PHASE_LIMIT)
    if coshi < caplo or coslo > 1:
        raise ValueError('field Ramsey record is outside the declared phase window')
    lo = _acos_bracket(min(F(1), coshi))[0]
    hi = PHASE_LIMIT if coslo <= caphi else _acos_bracket(coslo)[1]
    return {'energy_gap_GeV_interval': [str(lo/thi), str(hi/tlo)],
            'field_phase_interval': [str(lo), str(hi)],
            'clock': time, 'species_specific_scale_used': False,
            'source_energies_or_program_duration_used_by_inverse': False,
            'phase_window': '0 <= Delta*tau_effective <= 3/2',
            'absolute_common_energy_offset_inferred': False}


def sampling_envelope(shots, cell_error=F(1, 1000), divisions=256):
    """Simultaneous finite-shot bound for four clock cells and one Ramsey event.

    For independent repetitions of the same conditional native branch,
    Hoeffding and a union bound give 2*events*exp(-2*shots*cell_error**2).
    exp(-x) <= (1+x/divisions)**(-divisions) makes the reported upper
    bound rational and outward. No probability estimate is called one shot.
    """
    events = 5
    if any(type(x) is not int or x <= 0 for x in (shots, divisions)):
        raise ValueError('positive integer shot, event and subdivision counts required')
    error = _positive(cell_error, 'frequency error')
    if error >= 1:
        raise ValueError('frequency error must be less than one')
    failure = min(F(1), 2*events/(1+2*shots*error**2/divisions)**divisions)
    return {'independent_repetitions_per_conditional_branch': shots,
            'individual_frequency_error': str(error),
            'clock_binary_marginal_error': str(2*error),
            'four_cell_clock_total_variation_error': str(2*error),
            'field_Ramsey_binary_error': str(error),
            'simultaneous_failure_probability_upper': str(failure),
            'assumptions': 'independent repetitions of the same prepared conditional native branch; add these sampling errors to the separate source-service errors',
            'observed_shots_in_this_receipt': 0}


def si_calibration():
    """Exact SI definitions plus a rational Machin enclosure of pi."""
    def atan_bounds(x, terms=24):
        total = sum(((-1)**k*x**(2*k+1)/(2*k+1) for k in range(terms)), F())
        extra = (-1)**terms*x**(2*terms+1)/(2*terms+1)
        return min(total,total+extra), max(total,total+extra)
    a, b = atan_bounds(F(1,5)), atan_bounds(F(1,239))
    pi_lo, pi_hi = 16*a[0]-4*b[1], 16*a[1]-4*b[0]
    h, charge = F('6.62607015e-34'), F('1.602176634e-19')
    factor = h/(2*10**9*charge)
    return {'source': 'https://www.bipm.org/en/measurement-units/si-defining-constants',
            'Planck_h_J_s_exact': str(h), 'elementary_charge_C_exact': str(charge),
            'hbar_GeV_s_interval': [str(factor/pi_hi), str(factor/pi_lo)],
            'formula': 't_seconds=hbar_GeV_s*tau_effective_GeV_inverse',
            'pi_enclosure': 'Machin pi=16 atan(1/5)-4 atan(1/239), exact alternating rational bounds',
            'role': 'unit conversion and specified reference-energy calibration; not a fit of the model to a physical laboratory clock'}
