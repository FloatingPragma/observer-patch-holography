"""A fixed, uncalibrated Higgs feedback kernel from complete native records.

The corrected C/E instruments supply the exact event probability. A declared
five-state memory and coherent feedback turn its correlated history into a
positive family kernel. The family triplet is a gauge-trivial ancillary copy
of native geometric coordinates, not the colored current triplet. All memory
states, zero raw labels and the discarded old-memory record remain explicit.

This is a specified effective repair law, not a claim that unchanged C/E
selects the feedback. Its finite channel has noise. Only the proved small-step
limit gives the stationary averaged Hamiltonian; no mass packet enters it.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import expm
import sympy as s

from tools import strict_json
from . import native_gauge_object as objects
from . import native_code_transport as transport
from . import native_repair_mass_h2 as h2_control

ROOT = Path(__file__).resolve().parents[2]
COUNTS = (8, 32, 128, 512)
ARRAY_ADMISSION = 2e-10
DIAGNOSTIC_ADMISSION = 2e-9


def _clean(matrix):
    return objects._clean(matrix)


def _zero(matrix):
    return not _clean(matrix).todok()


def _entries(matrix):
    return [[str(x) for x in row] for row in _clean(matrix).tolist()]


@lru_cache(maxsize=1)
def native_records():
    """All128 actual C-left/E-right labels; no event is conditioned away."""
    result = []
    for cc in (0, 1):
        for cr in range(4):
            wc = transport.code.branch_weight('C', (2, 3), cc, cr)
            for ec in (0, 1):
                for er in range(8):
                    we = transport.PAIR_CODES['E'].branch_weight('E', (2, 10), ec, er)
                    weight = s.Rational(wc*we)
                    result.append((cc, cr, ec, er, cc == 0 and cr == 1, weight))
    if len(result) != 128 or sum(row[-1] for row in result) != 1:
        raise ArithmeticError('complete native record menu changed')
    p = sum(row[-1] for row in result if row[-2])
    if p != s.Rational(1, 12):
        raise ArithmeticError('positive-current C event changed')
    return tuple(result)


@lru_cache(maxsize=1)
def architecture():
    """Freeze native ordering before any empirical comparison."""
    data = objects.geometry(); B = data['frame']; metric = (B.T*B)[0, 0]
    candidates = []
    for g in sorted(data['rotations']):
        if objects.permutation_order(g) != 5:
            continue
        U = _clean(B.T*objects.permutation_matrix(g)*B/metric)
        if s.expand(s.trace(U)-(1+s.sqrt(5))/2) == 0:
            candidates.append((g, U))
    if len(candidates) != 12:
        raise ArithmeticError('native near-fivefold class changed')
    g, U = candidates[0]
    seed = next(i for i in range(12) if g[i] != i)
    v = B[seed, :].T; P = _clean(v*v.T/(v.T*v)[0])
    projectors = tuple(_clean(U**k*P*(U.T)**k) for k in range(5))
    if not _zero(U.T*U-s.eye(3)) or not _zero(U**5-s.eye(3)):
        raise ArithmeticError('native rotation lost its exact order')
    if any(not _zero(Q*Q-Q) or s.trace(Q) != 1 for Q in projectors):
        raise ArithmeticError('feedback coefficient is not a normalized projector')
    p = sum(row[-1] for row in native_records() if row[-2])
    stationary = s.ImmutableMatrix([(1-p)*p**k/(1-p**5) for k in range(5)])
    shift = s.zeros(5)
    for k in range(5): shift[(k+1) % 5, k] = 1
    reset = s.zeros(5); reset[0, :] = s.ones(1, 5)
    transition = p*shift+(1-p)*reset
    if transition*stationary != stationary or sum(stationary) != 1:
        raise ArithmeticError('record memory lost its stationary law')
    Y = _clean(sum((stationary[k]*projectors[k] for k in range(5)), s.zeros(3)))
    e2 = s.factor((1-s.trace(Y*Y))/2)
    determinant = s.factor(Y.det())
    expected_det = s.Rational(100224, 511709641)-s.Rational(418176, 12792741025)*s.sqrt(5)
    if (s.trace(Y) != 1 or e2 != s.Rational(1392, 22621)
            or s.simplify(determinant-expected_det) != 0):
        raise ArithmeticError('fixed source hierarchy polynomial changed')
    orbit = s.Matrix.hstack(v, U*v, U**2*v)
    if orbit.det() == 0:
        raise ArithmeticError('native seed orbit does not span the family triplet')
    return {'rotation': g, 'U': U, 'seed': seed, 'P': P,
            'projectors': projectors, 'p': p, 'stationary': stationary,
            'transition': s.ImmutableMatrix(transition), 'shift': shift,
            'Y': Y, 'e2': e2, 'determinant': s.expand(determinant)}


def averaging_bound(steps, time=1, *, initial_memory_stationary=False):
    """Diamond bound including arbitrary external reference and memory mixing.

    On CQ block trace norm, T eta=p S eta for sum-zero eta. Since T x is a
    positive trace-one state, eta'=(I-Pi)T x+(I-Pi)(Ad_U-I)T x gives
    ||eta'||<=p||eta||+4h. The marginal local error is <=4h^2+2hp||eta||.
    Summation gives 2hp||eta0||/(1-p)+4N h^2(1+p)/(1-p), with ||eta0||<=2.
    The stationary-start option requires x_k=pi_k sigma, not only a stationary
    classical marginal correlated with different quantum states.
    """
    if type(steps) is not int or steps <= 0 or type(time) not in (int, float) or not math.isfinite(time) or time < 0:
        raise ValueError('positive integer steps and finite nonnegative model time required')
    if type(initial_memory_stationary) is not bool:
        raise ValueError('stationary product preparation must be explicit')
    p = float(architecture()['p']); h = time/steps
    transient = 0 if initial_memory_stationary else 4*h*p/(1-p)
    return transient+4*steps*h*h*(1+p)/(1-p)


def feedback_unitaries(step):
    if type(step) not in (int, float) or not math.isfinite(step):
        raise ValueError('finite coherent model step required')
    X = np.asarray([[0, 1], [1, 0]], complex)
    result = []
    for P in architecture()['projectors']:
        P = np.asarray(P, complex)
        result.append(np.eye(6)+(math.cos(step)-1)*np.kron(np.eye(2), P)
                      -1j*math.sin(step)*np.kron(X, P))
    return np.asarray(result)


def raw_kraus(step):
    """All640 raw-record/old-memory Kraus matrices on the full30-state domain."""
    unitaries = feedback_unitaries(step)
    result = np.zeros((640, 30, 30), complex)
    labels = []
    index = 0
    for cc, cr, ec, er, advance, weight in native_records():
        for old in range(5):
            new = (old+1) % 5 if advance else 0
            result[index, 6*new:6*new+6, 6*old:6*old+6] = math.sqrt(float(weight))*unitaries[new]
            labels.append((cc, cr, ec, er, old, new))
            index += 1
    return result, np.asarray(labels, np.int64)


def joint_channel(step):
    """Exact factor contraction of the same complete recorded CQ instrument."""
    A = np.asarray(architecture()['transition'], float)
    U = feedback_unitaries(step)
    result = np.zeros((180, 180), complex)
    for new in range(5):
        pulse = np.kron(U[new], U[new].conj())
        for old in range(5):
            result[36*new:36*new+36, 36*old:36*old+36] = A[new, old]*pulse
    return result


def evolve_cq(blocks, steps, time=1):
    """Evolve all five classical-memory quantum blocks, without marginalizing."""
    averaging_bound(steps, time)
    blocks = np.asarray(blocks)
    if blocks.shape != (5, 6, 6) or blocks.dtype.kind != 'c' or not np.isfinite(blocks).all():
        raise ValueError('all five finite complex6-by6 memory blocks required')
    power = np.linalg.matrix_power(joint_channel(time/steps), steps)
    return (power@blocks.reshape(180)).reshape(5, 6, 6)


@lru_cache(maxsize=16)
def finite_channels(steps, time=1):
    """Five complete marginal channels, one per initial classical memory label."""
    averaging_bound(steps, time)
    partial = np.hstack([np.eye(36)]*5)
    power = np.linalg.matrix_power(joint_channel(time/steps), steps)
    result = np.stack([partial@power[:, 36*k:36*k+36] for k in range(5)])
    result.flags.writeable = False
    return result


def _choi(channel):
    return channel.reshape(6, 6, 6, 6).transpose(0, 2, 1, 3).reshape(36, 36)


def execute():
    data = architecture()
    Y = np.asarray(data['Y'], complex)
    H = np.kron(np.asarray([[0, 1], [1, 0]], complex), Y)
    U = expm(-1j*H)
    ideal = np.kron(U, U.conj())
    channels = [finite_channels(steps) for steps in COUNTS]
    records = native_records()
    labels = np.asarray([row[:4] for row in records], np.int64)
    weights = np.asarray([[int(row[-1].p), int(row[-1].q)] for row in records], np.int64)
    return {'channels': np.asarray(channels), 'averaged_unitary': U,
            'averaged_channel': ideal, 'kernel': Y,
            'projectors': np.asarray(data['projectors'], complex),
            'stationary': np.asarray(data['stationary'], float).reshape(5),
            'transition': np.asarray(data['transition'], float),
            'raw_labels': labels, 'raw_weights': weights,
            'raw_unitaries': feedback_unitaries(1/COUNTS[0])}


def _schema(arrays):
    expected = {'channels': ((4, 5, 36, 36), 'complex128'),
                'averaged_unitary': ((6, 6), 'complex128'),
                'averaged_channel': ((36, 36), 'complex128'),
                'kernel': ((3, 3), 'complex128'), 'projectors': ((5, 3, 3), 'complex128'),
                'stationary': ((5,), 'float64'), 'transition': ((5, 5), 'float64'),
                'raw_labels': ((128, 4), 'int64'), 'raw_weights': ((128, 2), 'int64'),
                'raw_unitaries': ((5, 6, 6), 'complex128')}
    from . import native_repair_mass_readout as readout
    from . import native_family_selection as family
    from . import native_mass_dressing as dressing
    from . import native_repair_mass_feedback as source_feedback
    from . import native_mass_record_cycles as cycles
    from . import native_mass_bath as bath
    from . import native_repair_mass_native_readout as native_readout
    for schema in (readout.ARRAY_SCHEMA, family.ARRAY_SCHEMA, dressing.SCHEMA,
                   source_feedback.ARRAY_SCHEMA, cycles.SCHEMA, bath.SCHEMA, native_readout.SCHEMA):
        if set(arrays) & set(schema): expected.update(schema)
    if set(arrays) != set(expected): raise ValueError('complete source array menu required')
    for key, (shape, dtype) in expected.items():
        value = arrays[key]
        if (type(value) is not np.ndarray or value.shape != shape
                or value.dtype != np.dtype(dtype) or not np.isfinite(value).all()):
            raise ValueError('invalid source array '+key)


def metrics(arrays):
    _schema(arrays)
    identity = np.eye(6).reshape(36)
    rows = []
    for steps, channels in zip(COUNTS, arrays['channels'], strict=True):
        trace_defect = max(float(np.linalg.norm(identity@C-identity)) for C in channels)
        min_choi = min(float(np.linalg.eigvalsh(_choi(C)).min()) for C in channels)
        errors = [6*float(np.linalg.norm(C-arrays['averaged_channel'])) for C in channels]
        bound = averaging_bound(steps)
        # This Frobenius certificate bounds the diamond error from the full Choi
        # matrices. It can be looser than the independent analytic bound.
        if trace_defect > ARRAY_ADMISSION or min_choi < -ARRAY_ADMISSION:
            raise ValueError('finite record channel is not complete positive trace preserving')
        actual = max(float(np.linalg.norm(C-arrays['averaged_channel'], 2)) for C in channels)
        if max(errors) > bound+ARRAY_ADMISSION:
            raise ValueError('complete Choi certificate exceeds correlated-memory averaging bound')
        rows.append({'steps': steps, 'step': 1/steps,
                     'analytic_diamond_bound_all_initial_classical_memories': bound,
                     'complete_Choi_Frobenius_diamond_upper_diagnostic': max(errors),
                     'superoperator_spectral_difference_diagnostic': actual,
                     'trace_preservation_defect_diagnostic': trace_defect,
                     'minimum_complete_Choi_eigenvalue_diagnostic': min_choi})
    return {'coherent_model_time': 1, 'all5_initial_classical_memories_checked': True,
            'independent_record_sampling_assumed': False,
            'complete_quantum_reference_domain_dimension': 6, 'finite_steps': rows}


def build():
    data = architecture(); arrays = execute(); results = metrics(arrays)
    from . import native_repair_mass_readout as readout
    from . import native_repair_mass_attachment as attachment
    from . import native_matter_theorems as matter_theorems
    from . import native_family_selection as family
    from . import native_mass_holonomy as holonomy
    from . import native_mass_dressing as dressing
    from . import native_repair_mass_feedback as source_feedback
    from . import native_family_potential as potential
    from . import native_mass_record_cycles as cycles
    from . import native_mass_bath as bath
    from . import native_repair_mass_native_readout as native_readout
    readout_report, readout_arrays = readout.execute()
    arrays.update(readout_arrays)
    family_report, family_arrays = family.execute()
    feedback_report, feedback_arrays = source_feedback.execute(family_arrays)
    native_readout_report, native_readout_arrays = native_readout.execute(feedback_arrays)
    dressing_report, dressing_arrays = dressing.execute()
    bath_report, bath_arrays = bath.execute()
    cycles_report, cycles_arrays = cycles.execute()
    arrays.update(family_arrays); arrays.update(feedback_arrays); arrays.update(dressing_arrays)
    arrays.update(bath_arrays); arrays.update(cycles_arrays)
    arrays.update(native_readout_arrays)
    _schema(arrays)
    spectrum = np.linalg.eigvalsh(arrays['kernel'])
    report = {'schema': 'oph.native_repair_mass.v1',
        'architecture': {'port_carrier': 'unchanged12-port incidence and proper rotations',
            'C_seam': [2, 3], 'E_seam': [2, 10],
            'rotation_selection': 'lexicographically first order5 native proper rotation with geometric trace(1+sqrt5)/2',
            'rotation': list(data['rotation']), 'seed_selection': 'first port moved by the selected rotation',
            'seed_port': data['seed'], 'family_action': 'gauge identity on the separately protected native geometric family band',
            'memory_states': 5, 'advance_event': 'C coin0, raw1, signed current+1',
            'all_other_raw_records': 'reset memory to0; identicallyzero labels retained',
            'ordering': 'correct native C/E, update classical memory, apply exp(-i step sigma_x tensor P_memory)',
            'old_memory_retained_as_reset_record': True,
            'new_coherent_Higgs_feedback_is_declared_architecture': True,
            'unchanged_C_E_alone_selects_the_feedback': False,
            'same_kernel_for_all_three_Higgs_channels': True},
        'exact_kernel': {'native_U': _entries(data['U']), 'native_seed_projector': _entries(data['P']),
            'p': str(data['p']), 'stationary_weights': list(map(str, data['stationary'])),
            'transition': _entries(data['transition']), 'Y': _entries(data['Y']),
            'trace': '1', 'second_elementary_symmetric': str(data['e2']),
            'determinant': str(data['determinant']),
            'formula': '(1-p)/(1-p^5) sum_k=0..4 p^k U^k P U^(-k)',
            'positive_definite_reason': 'positive stationary weights and span(v,Uv,U^2v)=C3',
            'eigenvalues_diagnostic': list(map(float, spectrum)),
            'ordered_ratios_diagnostic': list(map(float, spectrum/spectrum[0]))},
        'record_channel': {'raw_C_E_labels': 128, 'old_memory_labels': 5,
            'complete_instrument_branches': 640, 'joint_memory_Nambu_family_dimension': 30,
            'postselection': False, 'unread_memory_dephased_by_explicit_old_memory_record': True,
            'finite_step_is_unitary_on_marginal_quantum_system': False,
            'stationary_memory_product_required_for_zero_transient_bound': True},
        'averaging': {'time': 'dimensionless coherent model parameter t=N*step',
            'centered_memory_contraction': 'T eta=p S eta, norm<=p norm(eta)',
            'correlation_bound': 'norm(eta_n)<=p^n norm(eta0)+4 step/(1-p)',
            'all_initial_classical_memory_reference_states': True,
            'diamond_bound': '[4 p t+4(1+p)t^2]/[N(1-p)]',
            't1_bound': '56/(11*N)', 'stationary_product_t1_bound': '52/(11*N)',
            'independent_records_or_iid_memory_assumed': False},
        'execution': results, 'fixed_basis_readout': readout_report,
        'native_H2_domain_control': h2_control.report(),
        'same_current_Higgs_and_heldout_test': attachment.report(data['Y'], projectors=data['projectors']),
        'actual_source_finite_matter_theorems': matter_theorems.run(),
        'operational_family_band': family_report,
        'native_loop_holonomy': holonomy.report(),
        'finite_current_auxiliary_response': dressing_report,
        'native_feedback_source': feedback_report,
        'native_fixed_basis_readout': native_readout_report,
        'native_family_scalar_potential': potential.report(),
        'noncommuting_record_memory': cycles_report,
        'full_selector_incidence_bath': bath_report,
        'empirical_mass_or_coupling_inputs': [],
        'scope': {'fixed_repair_law_determines_its_dimensionless_family_kernel': True,
            'all_three_species_coefficients_separately_fitted': False,
            'absolute_mass_scale_derived': False, 'measured_mass_spectrum_claimed': False,
            'finite_native512_row_feedback_history_compiled_here': True,
            'complete_effective_CQ_instrument_executed': True,
            'physical_fermion_statistics_or_spin_derived_by_memory': False,
            'simulated_record_laws_are_finite_shot_experimental_data': False}}
    return report, arrays


def source_bindings():
    from . import native_repair_mass_readout as readout
    from . import native_repair_mass_attachment as attachment
    from . import native_matter_theorems as matter_theorems
    from . import native_family_selection as family
    from . import native_mass_holonomy as holonomy
    from . import native_mass_dressing as dressing
    from . import native_repair_mass_feedback as source_feedback
    from . import native_family_potential as potential
    from . import native_mass_record_cycles as cycles
    from . import native_mass_bath as bath
    from . import native_repair_mass_native_readout as native_readout
    paths = (set(transport.source_bindings()) | set(h2_control.source_bindings())
             | set(readout.source_bindings()) | set(attachment.source_bindings())
             | set(matter_theorems.source_bindings()) | set(family.source_bindings())
             | set(holonomy.source_bindings()) | set(dressing.source_bindings())
             | set(source_feedback.source_bindings()) | set(potential.source_bindings())
             | set(cycles.source_bindings()) | set(bath.source_bindings())
             | set(native_readout.source_bindings()))
    paths.update(('model/common_source/native_repair_mass.py',
                  'model/common_source/test_native_repair_mass.py'))
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


def _compare(a, b, path='receipt'):
    if type(a) is not type(b): raise ValueError(path+' type differs')
    if isinstance(a, dict):
        if set(a) != set(b): raise ValueError(path+' keys differ')
        for key in a: _compare(a[key], b[key], path+'.'+key)
    elif isinstance(a, list):
        if len(a) != len(b): raise ValueError(path+' length differs')
        for i, (x, y) in enumerate(zip(a, b, strict=True)): _compare(x, y, path+f'[{i}]')
    elif isinstance(a, float):
        derived_interval = any(name in path for name in ('eigenvalue_intervals', 'ordered_ratio_intervals'))
        resource_diagnostic = path.endswith(('.integrated_processor_norm_upper',
                                              '.maximum_numerical_valid_code_leakage'))
        tolerance = DIAGNOSTIC_ADMISSION if 'diagnostic' in path or derived_interval or resource_diagnostic else 0
        if not math.isfinite(b) or abs(a-b) > tolerance: raise ValueError(path+' numerical value differs')
    elif a != b: raise ValueError(path+' exact value differs')


def check(path):
    path = Path(path); report = strict_json.load(path); before = source_bindings()
    if report.get('source_bindings') != before: raise ValueError('source pins differ')
    artifact = path.with_suffix('.npz')
    if report.get('artifact') != {'name': artifact.name, 'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest()}:
        raise ValueError('binary custody differs')
    with np.load(artifact, allow_pickle=False) as archive:
        if len(archive.files) != len(set(archive.files)): raise ValueError('duplicate array key')
        arrays = {key: archive[key] for key in archive.files}
    from . import native_repair_mass_readout as readout
    from . import native_family_selection as family
    from . import native_mass_dressing as dressing
    from . import native_repair_mass_feedback as source_feedback
    from . import native_mass_record_cycles as cycles
    from . import native_mass_bath as bath
    from . import native_repair_mass_native_readout as native_readout
    _schema(arrays)
    for name, schema in (('readout', readout.ARRAY_SCHEMA), ('family', family.ARRAY_SCHEMA),
                         ('dressing', dressing.SCHEMA), ('feedback', source_feedback.ARRAY_SCHEMA),
                         ('cycles', cycles.SCHEMA), ('bath', bath.SCHEMA),
                         ('native_readout', native_readout.SCHEMA)):
        if not set(schema) <= set(arrays): raise ValueError('complete '+name+' arrays missing')
    readout_arrays = {key: arrays[key] for key in readout.ARRAY_SCHEMA}
    readout.validate_arrays(readout_arrays)
    inferred = readout.decode_records(arrays['readout_probabilities'], time=readout.TIME,
        entry_error_upper=readout.error_budget()['entry_error_upper'])
    inferred.pop('estimate')
    _compare(inferred, report['fixed_basis_readout']['inferred_kernel'], 'receipt.fixed_basis_readout.inferred_kernel')
    family_arrays = {key: arrays[key] for key in family.ARRAY_SCHEMA}
    _compare(family.native_metrics(family_arrays), report['operational_family_band']['actual_native_composition'],
             'receipt.operational_family_band.actual_native_composition')
    own_dressing = _bind_dressing(arrays, report['finite_current_auxiliary_response'])
    own_bath = _bind_bath(arrays, report['full_selector_incidence_bath'])
    own_native_inference = _bind_native_readout(arrays, report['native_fixed_basis_readout']['execution'])
    feedback_keys = set(source_feedback.ARRAY_SCHEMA) | set(source_feedback.NATIVE_SCHEMA)
    _compare(source_feedback.metrics({key: arrays[key] for key in feedback_keys}),
             report['native_feedback_source']['execution'], 'receipt.native_feedback_source.execution')
    _compare(cycles.metrics({key: arrays[key] for key in cycles.SCHEMA}),
             report['noncommuting_record_memory']['execution'], 'receipt.noncommuting_record_memory.execution')
    _compare(report['execution'], metrics(arrays))
    fresh, fresh_arrays = build(); _schema(fresh_arrays)
    # Derived readout intervals are verified above against their own retained
    # complete records. Fresh forward records and inferred arrays are compared
    # below at the stated numerical gate, avoiding an ulp/time amplification
    # in printed ratio endpoints becoming a platform-specific source change.
    fresh['fixed_basis_readout']['inferred_kernel'] = inferred
    # Native inference is bound to its complete retained pointer record laws,
    # including their total-variation reconstruction gate. Only the derived
    # intervals share that custody; fixed source/inversion budgets stay exact.
    fresh['native_fixed_basis_readout']['execution']['inferred_kernel'] = own_native_inference
    # The complete pole and residue summaries are reconstructed above from
    # their own retained factors. All law, domain and analytic-bound fields
    # are checked too. Fresh factors, poles and residues must pass the same
    # full-array gate below before these summaries can share that custody.
    fresh['finite_current_auxiliary_response'] = own_dressing
    fresh['full_selector_incidence_bath'] = own_bath
    _compare(fresh, {k: v for k, v in report.items() if k not in ('source_bindings', 'artifact')})
    if set(arrays) != set(fresh_arrays): raise ValueError('source array menu differs')
    for key in arrays:
        if arrays[key].dtype.kind in 'iu':
            if not np.array_equal(arrays[key], fresh_arrays[key]): raise ValueError('source integer array differs: '+key)
        elif np.max(abs(arrays[key]-fresh_arrays[key])) > ARRAY_ADMISSION:
            raise ValueError('source numerical array differs: '+key)
    if source_bindings() != before: raise ValueError('source changed during replay')


def _bind_dressing(arrays, retained_report):
    from . import native_mass_dressing as dressing
    own = dressing.report_from_arrays({key: arrays[key] for key in dressing.SCHEMA})
    _compare(own, retained_report, 'receipt.finite_current_auxiliary_response')
    return own


def _bind_bath(arrays, retained_report):
    from . import native_mass_bath as bath
    own = bath.report_from_arrays({key: arrays[key] for key in bath.SCHEMA})
    _compare(own, retained_report, 'receipt.full_selector_incidence_bath')
    return own


def _bind_native_readout(arrays, retained_execution):
    from . import native_repair_mass_native_readout as native_readout
    keys = set(native_readout.BASE_SCHEMA) | set(native_readout.SCHEMA)
    own = native_readout.metrics({key: arrays[key] for key in keys})
    _compare(own, retained_execution, 'receipt.native_fixed_basis_readout.execution')
    return own['inferred_kernel']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path); parser.add_argument('--check', type=Path)
    args = parser.parse_args()
    if bool(args.output) == bool(args.check): parser.error('choose --output or --check')
    if args.check:
        check(args.check); print('Fixed native record-memory Higgs kernel: PASS'); return
    before = source_bindings(); report, arrays = build()
    if source_bindings() != before: raise ValueError('source changed during generation')
    artifact = args.output.with_suffix('.npz'); np.savez_compressed(artifact, **arrays)
    report['source_bindings'] = before
    report['artifact'] = {'name': artifact.name, 'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest()}
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False)+'\n', encoding='utf-8')


if __name__ == '__main__': main()
