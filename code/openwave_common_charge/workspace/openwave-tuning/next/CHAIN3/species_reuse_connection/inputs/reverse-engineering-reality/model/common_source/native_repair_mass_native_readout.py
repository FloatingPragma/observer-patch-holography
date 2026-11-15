"""Fixed-basis recovery from complete native feedback and detector records.

Eight source-record memory updates precede each fresh computational probe.
One existing h=1/64 feedback step is followed by literal native decoders,
geometric inverse, fixed Pauli-Y interferometers and all128 pointer labels.
Every old/final syndrome remains in the defining source factors. The inverse
uses records and apparatus constants, never the kernel or its eigenvectors.
This recovers a specified finite feedback coefficient, not physical masses.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
import hashlib
import math

import numpy as np

from . import native_repair_mass_feedback as source
from . import native_code_transport as transport
from .mass_response import history
from ..observer.source_scalar import compiler

ROOT = source.ROOT
STEPS = 8
TIME = F(1, 64)
NORMALIZER = TIME-TIME**3/6+TIME**5/120
ADMISSION = 2e-8
ARRAY_ADMISSION = 2e-10
CONTRAST_NUMERICAL_ALLOWANCE = 1e-10
PROBES = tuple((f, g) for f in range(3) for g in range(3))
BASE_SCHEMA = dict(source.NATIVE_SCHEMA, **source.ARRAY_SCHEMA)
SCHEMA = {
    'native_readout_geometry_inverse': ((64, 64), 'complex128'),
    'native_readout_control_rotation': ((128, 128), 'complex128'),
    'native_readout_routes': ((9, 128, 128), 'complex128'),
    'native_readout_detectors': ((9, 128, 128), 'complex128'),
    'native_readout_pointer': ((2, 2, 2), 'complex128'),
    'native_readout_roundtrip_C': ((512, 8), 'complex128'),
    'native_readout_roundtrip_E': ((512, 8), 'complex128'),
    'native_readout_burnin_records': ((5, 128), 'float64'),
    'native_readout_burnin_transition': ((5, 5), 'float64'),
    'native_readout_burnin_memories': ((9, 5), 'float64'),
    'native_readout_records': ((9, 640, 128), 'float64'),
    'native_readout_kernel': ((3, 3), 'float64'),
    'native_readout_probe_labels': ((9, 2), 'int64'),
    'native_readout_pointer_labels': ((128,), 'int64')}


def _control_word():
    word = []
    for mask in range(64):
        label = ['Z' if mask >> j & 1 else 'I' for j in range(6)]+['X']
        word.extend(compiler.pauli_rotation(label, -math.pi/256))
    return word


def _route_word(f, g):
    a, b = 8*f, 64+8*g; word = []
    for bit in range(6):
        if a >> bit & 1:
            word += [compiler.pulse('I', -math.pi/2, [bit]), compiler.pulse('X', math.pi/2, [bit])]
    for bit in range(6):
        if (a ^ b) >> bit & 1: word += compiler.cnot(6, bit)
    return word


def _inverse_word():
    return compiler.inverse([event for stage in transport._preparation_stages() for event in stage])


def detector_sources():
    inverse_word = _inverse_word(); control_word = _control_word()
    inverse = compiler.execute(inverse_word, 6, np.eye(64, dtype=complex))
    control = compiler.execute(control_word, 7, np.eye(128, dtype=complex))
    routes = np.asarray([compiler.execute(_route_word(f, g), 7, np.eye(128, dtype=complex)) for f, g in PROBES])
    detectors = np.asarray([control@R@np.kron(np.eye(2), inverse) for R in routes])
    pointer = history.diagonal_read(np.eye(2, dtype=complex), [1, 0])
    arrays = {'native_readout_geometry_inverse': inverse,
        'native_readout_control_rotation': control, 'native_readout_routes': routes,
        'native_readout_detectors': detectors, 'native_readout_pointer': np.asarray(pointer.outputs)}
    resources = {'geometric_inverse': compiler.resources(inverse_word),
        'controlled_X_interferometer': compiler.resources(control_word),
        'routing_words': [compiler.resources(_route_word(f, g)) for f, g in PROBES],
        'local_binary_pointer': pointer.resources, 'pointer_repetitions': 7,
        'all_interferometer_input_columns': 128}
    return arrays, resources


def _local_channels(values):
    value = values.reshape(-1, 64, 8, 8)
    return np.einsum('rhia,rhjc->rijac', value, value.conj(), optimize=True)


def _pair_density(vector, C, E):
    rho = np.outer(vector, vector.conj()).reshape(8, 8, 8, 8)
    return np.einsum('ijac,abcd,klbd->ikjl', C, rho, E, optimize=True).reshape(64, 64)


def _roundtrip_density(rho, C, E):
    return np.einsum('ijac,nabmcd,klbd->nikmjl', C, rho.reshape(2, 8, 8, 2, 8, 8),
                     E, optimize=True).reshape(128, 128)


@lru_cache(maxsize=1)
def _pointer_path():
    shape = np.zeros([2]*14)
    effect = np.zeros((2, 2, 2))
    return np.einsum_path('abcdefgABCDEFG,hAa,iBb,jCc,kDd,lEe,mFf,nGg->hijklmn',
        shape, *([effect]*7), optimize='greedy')[0]


def pointer_probabilities(rho, pointer):
    effect = np.einsum('rha,rhb->rab', pointer.conj(), pointer)
    result = np.einsum('abcdefgABCDEFG,hAa,iBb,jCc,kDd,lEe,mFf,nGg->hijklmn',
        rho.reshape([2]*14), *([effect]*7), optimize=_pointer_path()).reshape(128)
    if np.max(abs(result.imag)) > CONTRAST_NUMERICAL_ALLOWANCE:
        raise ArithmeticError('complete pointer probabilities lost reality')
    return result.real


def _initial_densities(arrays, probe):
    C, E = (_local_channels(arrays['family_'+law]) for law in ('C', 'E'))
    vector = arrays['feedback_source_preparation'][:, probe]
    for c in range(8):
        for e in range(16):
            low = _pair_density(vector, C[c], E[e])
            rho = np.zeros((128, 128), complex); rho[:64, :64] = low
            yield c, e, rho


def burnin(arrays):
    # Every cycle receives fresh known probe work. The classical memory is
    # carried forward; all discarded work remains an explicit quantum output.
    C, E = (arrays['feedback_encoder_'+law] for law in ('C', 'E'))
    final_gram = np.kron(np.eye(2), np.kron(C.conj().T@C, E.conj().T@E))
    records = np.zeros((5, 128))
    for c, e, rho in _initial_densities(arrays, 0):
        cc, cr = arrays['family_C_labels'][c]
        for old in range(5):
            new = (old+1) % 5 if (cc, cr) == (0, 1) else 0
            U = arrays['feedback_source_unitaries'][new]
            records[old, 16*c+e] = np.trace(final_gram@U@rho@U.conj().T).real
    T = np.zeros((5, 5))
    for cc, cr, ec, er, old, new in arrays['feedback_branch_labels']:
        T[new, old] += records[old, 16*(4*cc+cr)+8*ec+er]
    memories = [np.eye(5)[0]]
    for _ in range(STEPS): memories.append(T@memories[-1])
    return records, T, np.asarray(memories)


def forward_records(arrays, memory):
    """All9 probes x640 native/memory labels x128 raw pointer labels."""
    post = [_local_channels(arrays['native_readout_roundtrip_'+law][None])[0] for law in ('C', 'E')]
    output = np.zeros((9, 640, 128))
    for f in range(3):
        for c, e, rho in _initial_densities(arrays, f):
            cc, cr = arrays['family_C_labels'][c]
            updates = [(old, (old+1) % 5 if (cc, cr) == (0, 1) else 0) for old in range(5)]
            by_new = {}
            for old, new in updates:
                if new not in by_new:
                    U = arrays['feedback_source_unitaries'][new]
                    work = _roundtrip_density(U@rho@U.conj().T, *post)
                    by_new[new] = [pointer_probabilities(D@work@D.conj().T, arrays['native_readout_pointer'])
                                  for D in arrays['native_readout_detectors'][3*f:3*f+3]]
                for g in range(3):
                    output[3*f+g, 128*old+16*c+e] = memory[old]*by_new[new][g]
    return output


def error_budget():
    """Reference-safe source, mixing and finite-step inversion allowances.

    The ideal reset law contracts zero-sum memories by p. Eight fresh-work
    cycles therefore contribute2*p^8 plus their CP-norm-amplified source
    error. One final prepared step and the complete detector contribute
    the remaining channel allowance. A contrast has operator norm one.
    Nine coefficient errors give spectral norm at most three times their
    entry bound. The real rank-one identity gives ||sym Z||<=1; alternating
    sine/cosine bounds supply the separate finite-step term.
    """
    a = F.from_float(ADMISSION); pair = a*(1+a)+a
    prepared = F.from_float(source._allowances()['prepared_six_input_dilation_allowance'])
    burn_step = (2+prepared)*prepared
    # Fresh preparation has operator norm<=1+a from its isometric target.
    # The two native branch factors together have norm<=1+a by their local
    # eight-input Gram bounds. The feedback unitary and final encoder pair
    # each have norm<=1+a from the same full-domain Gram checks. Squaring
    # this four-factor dilation norm bounds the CP norm of each retained
    # fresh-work step by(1+a)^8, including two preparation powers.
    burn_norm = (1+a)**8
    burn_error = burn_step*sum(burn_norm**j for j in range(STEPS))
    p = F(1, 12); mixing = 2*p**STEPS
    before_detector = (1+a)**2*a+(1+a)*pair+a
    pointer_error = (1+a)**7-1
    detector = (1+a)**8*pair+(1+a)**7*a+pointer_error
    dilation = (1+a)**9*before_detector+detector
    readout_channel = (2+dilation)*dilation
    total = readout_channel+(1+dilation)**2*(burn_error+mixing)
    finite_bias = TIME**2/2+TIME**7/(5040*NORMALIZER)
    radius = finite_bias+3*(total+F.from_float(CONTRAST_NUMERICAL_ALLOWANCE))/(2*NORMALIZER)
    values = {'burnin_source_channel_error': burn_error, 'burnin_ideal_mixing_l1_bound': mixing,
        'readout_source_channel_allowance': readout_channel, 'complete_contrast_error_allowance': total,
        'finite_step_symmetric_operator_bias_upper': finite_bias,
        'recovered_kernel_operator_error_upper': radius}
    return {key: math.nextafter(float(value), math.inf) for key, value in values.items()}


def decode_records(records, *, normalizer=float(NORMALIZER), operator_error=None):
    """Blind inverse: no architecture, supplied coefficients or eigbasis."""
    if type(records) is not np.ndarray or records.shape != (9, 640, 128) or not np.isfinite(records).all():
        raise ValueError('complete native pointer records are required')
    if type(normalizer) is not float or not math.isfinite(normalizer) or normalizer <= 0:
        raise ValueError('positive apparatus normalization required')
    if type(operator_error) is not float or not math.isfinite(operator_error) or operator_error < 0:
        raise ValueError('explicit nonnegative operator allowance required')
    probabilities = records.sum(axis=1)
    contrast = probabilities[:, 0]-probabilities[:, 64]
    raw = (-contrast/(2*normalizer)).reshape(3, 3).T
    kernel = (raw+raw.T)/2
    values = np.linalg.eigvalsh(kernel)
    intervals = [[float(x-operator_error), float(x+operator_error)] for x in values]
    ratios = ([[float(lo/intervals[0][1]), float(hi/intervals[0][0])] for lo, hi in intervals[1:]]
              if intervals[0][0] > 0 else None)
    return {'kernel': kernel, 'eigenvalue_intervals': intervals, 'ordered_ratio_intervals': ratios,
        'antisymmetric_inversion_diagnostic': float(np.linalg.norm(raw-raw.T, 2)),
        'normalization_defect_diagnostic': float(np.max(abs(probabilities.sum(axis=1)-1))),
        'matrix_eigenvectors_enter_probes_or_inverse': False}


def validate_arrays(arrays):
    if F.from_float(source.STEP) != TIME or source.SOURCE_ADMISSION != ADMISSION:
        raise ValueError('the fixed feedback clock or source admission changed')
    schema = dict(BASE_SCHEMA, **SCHEMA)
    if set(arrays) != set(schema): raise ValueError('complete native readout factors and records required')
    for key, (shape, dtype) in schema.items():
        value = arrays[key]
        if (type(value) is not np.ndarray or value.shape != shape or value.dtype != np.dtype(dtype)
                or not np.isfinite(value).all()): raise ValueError('invalid native readout array '+key)
    if (not np.array_equal(arrays['native_readout_probe_labels'], np.asarray(PROBES, np.int64))
            or not np.array_equal(arrays['native_readout_pointer_labels'], np.arange(128, dtype=np.int64))):
        raise ValueError('a fixed probe or raw pointer outcome disappeared')


def metrics(arrays):
    validate_arrays(arrays)
    source.metrics({key: arrays[key] for key in BASE_SCHEMA})
    target_inverse = transport.preparation_plan()[1].conj().T
    target_control = np.eye(128, dtype=complex)
    target_control[np.ix_([0, 64], [0, 64])] = np.asarray([[1, -1j], [-1j, 1]])/math.sqrt(2)
    errors = []
    for k, (f, g) in enumerate(PROBES):
        permutation = np.zeros((128, 128))
        for old in range(128):
            out = old ^ (8*f)
            if out >> 6 & 1: out ^= (8*f) ^ (8*g)
            permutation[out, old] = 1
        ideal = target_control@permutation@np.kron(np.eye(2), target_inverse)
        D = arrays['native_readout_detectors'][k]
        errors.append(float(np.linalg.norm(D-ideal)))
        composed = arrays['native_readout_control_rotation']@arrays['native_readout_routes'][k]@np.kron(np.eye(2), arrays['native_readout_geometry_inverse'])
        if np.max(abs(D-composed)) > ARRAY_ADMISSION: raise ValueError('detector composition differs')
        if np.linalg.norm(D.conj().T@D-np.eye(128), 2) > ADMISSION:
            raise ValueError('a full detector input/reference direction was lost')
    pointer = arrays['native_readout_pointer']
    target = np.asarray([[[1, 0], [0, 0]], [[0, 0], [0, 1]]])
    pointer_error = float(np.linalg.norm(pointer-target))
    if np.linalg.norm(np.einsum('rha,rhb->ab', pointer.conj(), pointer)-np.eye(2), 2) > ADMISSION:
        raise ValueError('local source pointer is incomplete')
    roundtrip_errors = []
    for law in ('C', 'E'):
        M = arrays['native_readout_roundtrip_'+law]; target = np.eye(512, 6, dtype=complex)
        roundtrip_errors.append(float(np.linalg.norm(M[:, :6]-target)))
        if np.linalg.norm(M.conj().T@M-np.eye(8), 2) > ADMISSION:
            raise ValueError('a native decoded padded input was lost')
    if max(errors+roundtrip_errors+[pointer_error]) > ADMISSION:
        raise ValueError('literal detector source exceeds its fixed component allowance')
    records, T, memories = burnin(arrays)
    for name, expected in [('native_readout_burnin_records', records), ('native_readout_burnin_transition', T),
                           ('native_readout_burnin_memories', memories)]:
        if np.max(abs(arrays[name]-expected)) > ARRAY_ADMISSION: raise ValueError('source burn-in differs: '+name)
    recorded = arrays['native_readout_records']
    computed = forward_records(arrays, memories[-1])
    record_l1_error = float(np.max(np.sum(abs(recorded-computed), axis=(1, 2))))
    if record_l1_error > CONTRAST_NUMERICAL_ALLOWANCE:
        raise ValueError('complete native records differ from their source factors in total variation')
    if recorded.min() < -CONTRAST_NUMERICAL_ALLOWANCE:
        raise ValueError('a native pointer probability is negative beyond numerical allowance')
    inferred = decode_records(recorded, operator_error=error_budget()['recovered_kernel_operator_error_upper'])
    if inferred['normalization_defect_diagnostic'] > error_budget()['complete_contrast_error_allowance']+CONTRAST_NUMERICAL_ALLOWANCE:
        raise ValueError('complete native detector probability normalization exceeds the source allowance')
    if np.max(abs(inferred['kernel']-arrays['native_readout_kernel'])) > ARRAY_ADMISSION:
        raise ValueError('recovered kernel differs from its own native detector records')
    inferred.pop('kernel')
    return {'maximum_full128_detector_error_diagnostic': max(errors),
        'native_roundtrip_errors_diagnostic': roundtrip_errors,
        'local_pointer_dilation_error_diagnostic': pointer_error,
        'complete_record_l1_reconstruction_error_diagnostic': record_l1_error,
        'burnin_probability_normalization_defect_diagnostic': float(np.max(abs(records.sum(axis=1)-1))),
        'all128_pointer_labels_and_native_syndrome_outputs_retained': True,
        'inferred_kernel': inferred}


def execute(base_arrays=None):
    base = source.execute()[1] if base_arrays is None else base_arrays
    arrays = {key: base[key] for key in BASE_SCHEMA}
    detectors, resources = detector_sources(); arrays.update(detectors)
    for law in ('C', 'E'):
        arrays['native_readout_roundtrip_'+law] = transport.PAIR_CODES[law].decoding(arrays['feedback_encoder_'+law])
    records, T, memories = burnin(arrays)
    arrays.update({'native_readout_burnin_records': records, 'native_readout_burnin_transition': T,
        'native_readout_burnin_memories': memories, 'native_readout_probe_labels': np.asarray(PROBES, np.int64),
        'native_readout_pointer_labels': np.arange(128, dtype=np.int64)})
    arrays['native_readout_records'] = forward_records(arrays, memories[-1])
    arrays['native_readout_kernel'] = decode_records(arrays['native_readout_records'],
        operator_error=error_budget()['recovered_kernel_operator_error_upper'])['kernel']
    report = {'classification': 'fixed-basis native-record recovery of a specified feedback coefficient',
        'dimensionless_model_step_exact': str(TIME), 'rational_contrast_normalizer_exact': str(NORMALIZER),
        'probe_pairs': [list(pair) for pair in PROBES], 'burnin_steps': STEPS, 'initial_memory': 0,
        'burnin_uses_fresh_known_native_work_before_each_science_probe': True,
        'stationary_weights_inserted_as_preparation_inputs': False,
        'complete_domains': {'native_memory_labels': 640, 'raw_pointer_labels': 128,
            'joint_records_per_probe': 81920, 'fixed_probe_settings': 9,
            'old_decoder_syndrome_dimension': 4096, 'final_decoder_syndrome_dimension': 4096,
            'local_pointer_quantum_output_dimension': 2, 'pointer_quantum_output_dimension': 128,
            'syndromes_are_classical_measured_labels': False,
            'all_native_quantum_rows_defined_by_retained_factors': True},
        'finite_step_identity': '-contrast/(2*sin(h))=Y+(cos(h)-1)*Z; Z=sum pi*v*(v^3)^T, norm(sym Z)<=1',
        'inverse_certificate': [
            'real rank-one P=v*v^T gives contrast_fg=-2*sin(h)*sum pi*P_gf*(1+(cos(h)-1)*P_ff)',
            'norm(v*(v^3)^T)<=1, so norm(sym Z)<=1 for every probability mixture',
            'sin(h)<=s0=h-h^3/6+h^5/120 and s0-sin(h)<=h^7/5040',
            'norm(recovered symmetric kernel-Y)<=h^2/2+h^7/(5040*s0)+3*contrast_error/(2*s0)',
            'fresh probe preparation after burn-in makes memory and probe a product without assuming stationary memory',
            'all old/final native syndromes and all pointer outcomes enter the source factors and probability contractions'],
        'source_and_inversion_budget': error_budget(), 'source_resources': resources,
        'execution': metrics(arrays),
        'scope': {'empirical_packet_enters_apparatus_or_inverse': False,
            'mass_eigenvectors_or_pole_values_enter_probes': False,
            'complete_detector_record_laws_are_experimental_finite_shot_samples': False,
            'target': 'real symmetric averaged coefficient of the existing five-memory rank-one law',
            'Nambu_control_and_pointer_are_declared_register_operations': True,
            'charge_conserving_physical_CAR_detector_derived': False,
            'physical_interacting_particle_poles_or_absolute_mass_scale_recovered': False,
            'numerical_admissions_are_outward_hardware_certificates': False,
            'fresh_blank_and_coherent_control_services_are_declared': True}}
    return report, arrays


def report(): return execute()[0]


def source_bindings():
    names = set(source.source_bindings()) | {
        'model/common_source/native_repair_mass_native_readout.py',
        'model/common_source/test_native_repair_mass_native_readout.py'}
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(names)}
