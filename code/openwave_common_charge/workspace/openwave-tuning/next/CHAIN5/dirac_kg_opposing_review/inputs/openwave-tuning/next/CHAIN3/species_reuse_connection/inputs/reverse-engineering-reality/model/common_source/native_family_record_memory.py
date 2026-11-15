"""Complete C/E-record family memory with an explicitly supplied feedback service.

The native events are C=(coin0,raw1) and E=(coin0,raw4 or raw5).
All128 source labels and all12 old-memory inputs remain in the instrument.
A C event rotates the port label; a joint event then applies E's rotation.
Other records, and proposed core destinations, reset to port2.

The classical memory law and its positive projector average are finite results.
The additional pulse exp(-ih X tensor P_new) and fixed Pauli detectors are
declared logical services, not literal compiled native pulse histories. Their
coefficients are dimensionless. No fermion mass, generation identification,
physical clock calibration or coherent U(1) phase follows from this module.
"""
from __future__ import annotations

from functools import lru_cache
import math

import numpy as np
import sympy as s

from . import native_gauge_object as objects
from . import native_repair_mass as feedback
from . import native_repair_mass_readout as fixed_readout

CORE = (0, 1)
RESET = 2
READOUT_STEPS = (1 / 8, 1 / 4)
BURNIN = 8
NUMERICAL_CONTRAST_ALLOWANCE = 1e-12
ADMISSION = 2e-10
RAW_LABELS = tuple((cc, cr, ec, er) for cc in (0, 1) for cr in range(4)
                   for ec in (0, 1) for er in range(8))


def _number(value):
    try:
        value = s.Rational(str(value)) if isinstance(value, (float, np.floating)) else s.sympify(value)
    except (s.SympifyError, TypeError, ValueError, OverflowError) as exc:
        raise ValueError('a finite real event factor with no free parameters is required') from exc
    if (not isinstance(value, s.Expr) or value.free_symbols
            or value.is_number is not True or value.is_real is not True
            or value.is_finite is not True):
        raise ValueError('a finite real event factor with no free parameters is required')
    return value


def _probability(value, *, strict_upper=False):
    value = _number(value)
    # Exact probabilities must not cross a domain boundary through float
    # rounding or underflow. Unknown orderings do not certify the domain.
    complement = 1 - value
    upper_admitted = complement.is_positive if strict_upper else complement.is_nonnegative
    if value.is_nonnegative is not True or upper_admitted is not True:
        raise ValueError('event probability outside its declared domain')
    return value


@lru_cache(maxsize=1)
def geometry():
    data = objects.geometry()
    rotations = []
    for a, b in ((2, 3), (2, 10)):
        candidates = [g for g in data['rotations'] if objects.permutation_order(g) == 2
                      and g[a] == b and g[b] == a]
        if len(candidates) != 1:
            raise ArithmeticError('native seam half-turn is not unique')
        rotations.append(candidates[0])
    return data, tuple(rotations)


def _destination(old, c_event, e_event, reset, gC, gE):
    new = gC[old] if c_event else reset
    if c_event and e_event:
        new = gE[new]
    return reset if new in CORE else new


def _validated_records(records):
    records = tuple(records)
    if len(records) != 128 or tuple(tuple(row[:4]) for row in records) != RAW_LABELS:
        raise ValueError('all128 ordered raw C/E labels, including zero labels, are required')
    result = []
    for row in records:
        if len(row) != 6:
            raise ValueError('complete native record tuples are required')
        cc, cr, ec, er, advance, weight = row
        if type(advance) is not bool or advance != ((cc, cr) == (0, 1)):
            raise ValueError('C advance flag differs from the actual raw event')
        weight = _probability(weight)
        result.append((cc, cr, ec, er, advance, weight))
    if s.simplify(sum(row[-1] for row in result) - 1) != 0:
        raise ValueError('complete raw weights must sum to one without conditioning')
    return tuple(result)


def transition_from_records(records, *, reset=RESET):
    """Construct the full12-input transition from an explicitly supplied record law.

    This validates the native label/event interface. It does not certify that
    replacement weights came from a physical instrument. architecture() binds
    the weights to the actual unchanged native C/E source.
    """
    if type(reset) is not int or reset not in set(range(12)) - set(CORE):
        raise ValueError('a non-core reset port is required')
    records = _validated_records(records)
    _, (gC, gE) = geometry()
    T = s.zeros(12)
    labels = []
    for cc, cr, ec, er, advance, weight in records:
        for old in range(12):
            new = _destination(old, advance, ec == 0 and er in (4, 5), reset, gC, gE)
            T[new, old] += weight
            labels.append((cc, cr, ec, er, old, new))
    return s.ImmutableMatrix(T), tuple(labels)


def _projectors(conjugate):
    data, _ = geometry()
    B = data['frame']
    if conjugate:
        B = B.subs(s.sqrt(5), -s.sqrt(5))
    return tuple(s.ImmutableMatrix(objects._clean(B[i, :].T * B[i, :] /
                 (B[i, :] * B[i, :].T)[0])) for i in range(12))


def _memory(T, p, q, reset, conjugate):
    if reset == 2:
        # Solve first in a=p(1-q),b=pq. The support is2,3,9,10, but no
        # transient or zero-probability input is removed from T/instrument.
        a, b = p * (1 - q), p * q
        heavy = 1 / (1 + a + a * b / (1 - a))
        stationary = s.zeros(12, 1)
        stationary[2], stationary[3] = heavy, a * heavy
        stationary[10] = a * b * heavy / (1 - a * a)
        stationary[9] = a * stationary[10]
    else:
        A = s.Matrix(s.eye(12) - T)
        A[11, :] = s.ones(1, 12)
        rhs = s.zeros(12, 1)
        rhs[11] = 1
        stationary = A.inv() * rhs
    # Rational native factors simplify exactly; preserve symbolic interaction
    # factors without a global trigonometric simplification/approximation.
    if p.is_Rational and q.is_Rational:
        stationary = stationary.applyfunc(s.cancel)
    projectors = _projectors(conjugate)
    Y = sum((stationary[k] * projectors[k] for k in range(12)), s.zeros(3))
    if p.is_Rational and q.is_Rational:
        Y = objects._clean(Y)
    _, (gC, gE) = geometry()
    return {'T': s.ImmutableMatrix(T), 'stationary': s.ImmutableMatrix(stationary),
            'projectors': projectors, 'Y': s.ImmutableMatrix(Y), 'p': p, 'q': q,
            'gC': gC, 'gE': gE, 'reset': reset, 'conjugate': conjugate,
            'event_factor_domain': '0<=p<1,0<=q<=1; unique stationary reset chain',
            'feedback_service_declared': True,
            'physical_fermion_mass_map_derived': False}


def architecture_from_factors(p, q, reset=RESET, conjugate=False):
    """Reuse the memory with explicitly supplied independent event factors.

    The caller owns their interaction-instrument ancestry and independence
    proof. Exact SymPy expressions are retained. Float inputs are represented
    by their decimal rational values, not promoted to exact source results.
    """
    p, q = _probability(p, strict_upper=True), _probability(q)
    if type(reset) is not int or reset not in set(range(12)) - set(CORE):
        raise ValueError('a non-core reset port is required')
    if type(conjugate) is not bool:
        raise ValueError('the conjugate-triplet convention must be explicit')
    _, (gC, gE) = geometry()
    T = s.zeros(12)
    for old in range(12):
        for c, e, weight in ((False, False, 1-p), (True, False, p*(1-q)), (True, True, p*q)):
            T[_destination(old, c, e, reset, gC, gE), old] += weight
    result = _memory(T, p, q, reset, conjugate)
    result['event_factor_origin'] = 'explicit caller-supplied independent factors; not derived by this constructor'
    return result


@lru_cache(maxsize=1)
def architecture():
    records = _validated_records(feedback.native_records())
    p = sum(row[-1] for row in records if row[4])
    q = sum(row[-1] for row in records if row[2] == 0 and row[3] in (4, 5))
    joint = sum(row[-1] for row in records if row[4] and row[2] == 0 and row[3] in (4, 5))
    if (p, q, joint) != (s.Rational(1, 12), s.Rational(1, 6), s.Rational(1, 72)):
        raise ArithmeticError('actual native C/E event weights changed')
    T, labels = transition_from_records(records)
    result = _memory(T, p, q, RESET, False)
    if T * result['stationary'] != result['stationary'] or sum(result['stationary']) != 1:
        raise ArithmeticError('exact stationary record memory failed')
    result.update(records=records, labels=labels, coincidence=joint,
                  event_factor_origin='all128 unchanged native C/E raw-record weights')
    return result


def kernel_diagnostics(data=None):
    """Uncalibrated dimensionless coefficient diagnostics; no empirical packet."""
    data = architecture() if data is None else data
    Y = data['Y']
    eigenvalues = np.linalg.eigvalsh(np.asarray(Y.evalf(30), float))
    return {'trace_exact': str(s.simplify(s.trace(Y))),
            'second_elementary_symmetric_exact': str(s.factor((s.trace(Y)**2-s.trace(Y*Y))/2)),
            'determinant_exact': str(s.factor(Y.det())),
            'characteristic_polynomial_exact': str(Y.charpoly(s.Symbol('x')).as_expr()),
            'eigenvalues_diagnostic': eigenvalues.tolist(),
            'ratios_diagnostic': (eigenvalues/eigenvalues[0]).tolist() if eigenvalues[0] > 0 else None,
            'eigenvalues_are_dimensionless_coupling_coefficients': True,
            'physical_pole_masses_inferred': False}


def feedback_unitaries(step, data=None):
    """Declared logical projector pulses on Nambu2 tensor family3."""
    if type(step) not in (int, float) or not math.isfinite(step):
        raise ValueError('finite dimensionless coherent step required')
    data = architecture() if data is None else data
    X = np.array([[0, 1], [1, 0]], complex)
    return np.asarray([np.eye(6)+(math.cos(step)-1)*np.kron(np.eye(2), np.asarray(P, complex))
                       -1j*math.sin(step)*np.kron(X, np.asarray(P, complex)) for P in data['projectors']])


def raw_instrument_factors(step):
    """All1536 factors K=sqrt(w)|new><old| tensor U_new, including zeros.

    The small6-by6 factor and explicit old/new labels specify each full72-by72
    Kraus operator without allocating a mostly zero1536-by72-by72 array.
    Old memory is an explicit discarded record, so no memory coherences are
    silently retained by this classical-quantum instrument.
    """
    data = architecture()
    U = feedback_unitaries(step, data)
    factors = np.asarray([math.sqrt(float(data['records'][i//12][-1])) * U[row[5]]
                          for i, row in enumerate(data['labels'])])
    return factors, np.asarray(data['labels'], np.int64)


def joint_channel(step, data=None):
    """Full432-dimensional CQ superoperator, not independent memory samples."""
    data = architecture() if data is None else data
    T = np.asarray(data['T'].evalf(30), float)
    U = feedback_unitaries(step, data)
    return np.block([[T[new, old]*np.kron(U[new], U[new].conj())
                      for old in range(12)] for new in range(12)])


def averaging_bound(steps, time=1, *, initial_memory_stationary=False, data=None):
    """Reference-safe correlated-memory bound in dimensionless model time.

    T=(1-p)R+pS with S stochastic even with forbidden-core resets, hence
    ||T eta||1<=p||eta||1 on zero-sum operator-valued memory blocks. The
    existing bounded-projector argument gives4pt/[N(1-p)] for an arbitrary
    initial CQ memory and4(1+p)t²/[N(1-p)]. A stationary PRODUCT removes
    the transient term; stationary memory alone does not remove correlations.
    """
    if type(steps) is not int or steps <= 0 or type(time) not in (int, float) or not math.isfinite(time) or time < 0:
        raise ValueError('positive steps and finite nonnegative model time required')
    if type(initial_memory_stationary) is not bool:
        raise ValueError('stationary product preparation must be explicit')
    p = float((architecture() if data is None else data)['p'])
    transient = 0 if initial_memory_stationary else 4*p*time/(steps*(1-p))
    return transient+4*(1+p)*time*time/(steps*(1-p))


def evolve_cq(blocks, steps, time=1, data=None):
    averaging_bound(steps, time, data=data)
    if (type(blocks) is not np.ndarray or blocks.shape != (12, 6, 6)
            or blocks.dtype != np.dtype('complex128') or not np.isfinite(blocks).all()):
        raise ValueError('all12 finite complex6-by6 CQ blocks are required')
    power = np.linalg.matrix_power(joint_channel(time/steps, data), steps)
    return (power @ blocks.reshape(432)).reshape(12, 6, 6)


def stationary_currents(data=None):
    data = architecture() if data is None else data
    flux = data['T'] * s.diag(*list(data['stationary']))
    current = flux-flux.T
    return {'flux': flux, 'current': current, 'detailed_balanced': current == s.zeros(12),
            'classical_current_is_coherent_U1_phase': False,
            'reversal_scope': 'this declared classical memory; no claim about the base carrier or full environment'}


def _clock_weights(steps):
    if (not isinstance(steps, (tuple, list, np.ndarray)) or len(steps) != 2
            or any(not math.isfinite(float(h)) or not 0 < float(h) < math.pi/2 for h in steps)):
        raise ValueError('two declared positive small apparatus steps are required')
    h, k = map(float, steps)
    denominator = math.cos(k)-math.cos(h)
    if abs(denominator) < 1e-8:
        raise ValueError('distinct well-separated clock cosines are required')
    return ((math.cos(k)-1)/denominator, -(math.cos(h)-1)/denominator)


def forward_records(*, burnin=BURNIN, steps=READOUT_STEPS):
    """Complete fixed-basis two-clock records after finite classical burn-in.

    Burn-in contracts the actual T from reset2 with no quantum feedback on
    the fresh probe work. The final pulse keeps old memory, all128 raw source
    labels and all three Pauli-Y pointer outcomes. No averaged Y exponential,
    eigbasis or coefficient enters a preparation or a detector.
    """
    if type(burnin) is not int or burnin < 0:
        raise ValueError('a nonnegative integer burn-in is required')
    _clock_weights(steps)
    data = architecture()
    memory = np.eye(12)[RESET]
    T = np.asarray(data['T'], float)
    memories = [memory.copy()]
    for _ in range(burnin):
        memory = T @ memory
        memories.append(memory.copy())
    # Aggregate identical new-memory pulses only in computing the factors;
    # no raw, old-memory or complement outcome is removed from the output.
    records = np.empty((2, 9, 1536, 3), float)
    for clock, step in enumerate(steps):
        U = feedback_unitaries(float(step), data)
        outcomes = np.empty((9, 12, 3), float)
        for f, rho in enumerate(fixed_readout.fixed_inputs()):
            for new in range(12):
                output = U[new] @ rho @ U[new].conj().T
                for g in range(3):
                    law = np.einsum('ij,oji->o', output, fixed_readout.fixed_effects(f, g))
                    if np.max(abs(law.imag)) > NUMERICAL_CONTRAST_ALLOWANCE:
                        raise ArithmeticError('pointer law is not real')
                    outcomes[3*f+g, new] = law.real
        for i, row in enumerate(data['labels']):
            records[clock, :, i] = memory[row[4]] * float(data['records'][i//12][-1]) * outcomes[:, row[5]]
    return {'probabilities': records, 'labels': np.asarray(data['labels'], np.int64),
            'burnin_memories': np.asarray(memories), 'burnin': burnin, 'steps': tuple(map(float, steps))}


def readout_error_budget(*, burnin=BURNIN, steps=READOUT_STEPS):
    """Exact two-clock rank-one inverse and separate finite-mixing error.

    For one P²=P pulse, contrast(f,g)=-2 sin(h)P[g,f]
    *[1+(cos(h)-1)P[f,f]]. Thus A_h=-contrast/(2sin(h))
    =Y+(cos(h)-1)Z for the SAME finite-memory distribution at both clocks.
    The declared weights cancel Z exactly, without a small-step truncation.
    The final distribution after burnin+1 transitions differs from stationary
    by at most2p^(burnin+1). A contrast error epsilon at each clock contributes
    sum|weight_h|epsilon/(2|sin(h)|) per entry, at most three times that in
    operator norm. Numerical allowances are stated assumptions, not outward
    source-error certificates or finite-shot confidence intervals.
    """
    if type(burnin) is not int or burnin < 0:
        raise ValueError('valid finite mixing is required')
    weights = _clock_weights(steps)
    mixing = 2 * float(architecture()['p'])**(burnin+1)
    amplification = sum(abs(w)/(2*math.sin(h)) for w, h in zip(weights, steps, strict=True))
    numerical = 3*NUMERICAL_CONTRAST_ALLOWANCE*amplification
    radius = mixing+numerical
    return {'operator_error_upper': math.nextafter(radius, math.inf),
            'finite_rank_one_operator_bias_upper': 0.,
            'mixing_l1_upper': mixing, 'numerical_operator_allowance': numerical,
            'contrast_to_entry_error_amplification': amplification,
            'two_clock_inverse_weights': list(weights),
            'numerical_contrast_allowance': NUMERICAL_CONTRAST_ALLOWANCE,
            'numerical_allowance_is_outward_or_physical_certificate': False,
            'burnin_transitions': burnin, 'final_record_transition_additional': 1,
            'dimensionless_coherent_steps': list(steps),
            'shot_noise_error_included': False}


def decode_records(probabilities, *, steps, operator_error):
    """Record-only inverse; no architecture, coefficient or eigbasis input."""
    if (type(probabilities) is not np.ndarray or probabilities.shape != (2, 9, 1536, 3)
            or probabilities.dtype != np.dtype('float64') or not np.isfinite(probabilities).all()
            or probabilities.min() < -NUMERICAL_CONTRAST_ALLOWANCE):
        raise ValueError('all raw-memory and pointer records are required')
    weights = _clock_weights(steps)
    if (type(operator_error) not in (int, float) or not math.isfinite(operator_error) or operator_error < 0):
        raise ValueError('declared apparatus step and operator allowance are required')
    if np.max(abs(probabilities.sum(axis=(2, 3))-1)) > NUMERICAL_CONTRAST_ALLOWANCE:
        raise ValueError('complete terminal records are required without renormalization')
    law = probabilities.sum(axis=2)
    raw = sum(w*(-(law[c, :, 0]-law[c, :, 1])/(2*math.sin(h))).reshape(3, 3).T
              for c, (w, h) in enumerate(zip(weights, steps, strict=True)))
    estimate = (raw+raw.T)/2
    values = np.linalg.eigvalsh(estimate)
    return {'estimate': estimate, 'operator_error_upper': operator_error,
            'eigenvalue_intervals': [[float(v-operator_error), float(v+operator_error)] for v in values],
            'antisymmetric_inversion_diagnostic': float(np.linalg.norm(raw-raw.T, 2)),
            'effective_coefficient_is_physical_pole_mass': False}


def validate_records(records):
    """Bind retained laws to the actual source, apparatus and finite burn-in.

    A normalized replacement law, even with a matching newly decoded kernel,
    is not accepted as the source execution. There is no hash-only admission.
    The blind inverse deliberately remains separate from this source replay.
    """
    expected = {'probabilities', 'labels', 'burnin_memories', 'burnin', 'steps'}
    if not isinstance(records, dict) or set(records) != expected:
        raise ValueError('complete record factors and apparatus settings are required')
    if records['burnin'] != BURNIN or tuple(records['steps']) != READOUT_STEPS:
        raise ValueError('the declared fixed burn-in or clock protocol changed')
    fresh = forward_records()
    for key in ('probabilities', 'burnin_memories', 'labels'):
        value, target = records[key], fresh[key]
        if (type(value) is not np.ndarray or value.shape != target.shape
                or value.dtype != target.dtype or not np.isfinite(value).all()):
            raise ValueError('invalid complete record factor '+key)
        if key == 'labels':
            if not np.array_equal(value, target):
                raise ValueError('raw or old/new memory labels changed')
        elif np.max(abs(value-target)) > NUMERICAL_CONTRAST_ALLOWANCE:
            raise ValueError('retained law differs from actual source replay: '+key)
    decode_records(records['probabilities'], steps=READOUT_STEPS,
                   operator_error=readout_error_budget()['operator_error_upper'])


def execute():
    """Run this candidate with no comparison data or persistent output writes."""
    data = architecture()
    records = forward_records()
    bounds = readout_error_budget()
    inferred = decode_records(records['probabilities'], steps=records['steps'],
                              operator_error=bounds['operator_error_upper'])
    return {'schema': 'oph.native_family_record_memory.v1', 'kernel': kernel_diagnostics(data),
            'stationary_exact': list(map(str, data['stationary'])),
            'event_probabilities_exact': {'C': str(data['p']), 'E': str(data['q']), 'joint': str(data['coincidence'])},
            'complete_raw_labels': 128, 'old_memory_inputs': 12, 'complete_instrument_branches': 1536,
            'error_budget': bounds, 'inferred_kernel': inferred,
            'persistence': {'axis_weight_ratios': '1:a:a*b/(1-a), a=p(1-q), b=pq',
                            'light_shortest_history': 'C then EC: two memory updates and three rare C/E events',
                            'power_p3_requires_q_proportional_to_p': True,
                            'stationary_probabilities_are_couplings_only_by_declared_feedback': True},
            'scope': {'complete_effective_CQ_instrument_executed': True,
                      'literal_native_pulse_history_compiled': False,
                      'fixed_Pauli_detectors_supplied': True,
                      'classical_irreversibility_derives_U1_phase': False,
                      'physical_pole_mass_map_derived': False,
                      'empirical_mass_inputs': [], 'measured_spectrum_claimed': False}}, records
