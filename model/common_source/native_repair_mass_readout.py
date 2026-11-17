"""Fixed-basis records of the complete finite-memory Higgs feedback channel.

Three computational inputs, each used in three detector settings, infer all
nine entries of the real family kernel. Neither preparation nor detector
uses its eigenvectors, eigenvalues or a measured mass. The forward law keeps
all five classical memory outputs and executes the correlated finite-step
channel. Its averaged Hamiltonian is used only in the independent error
bound connecting these finite records to the inferred kernel.
"""
from __future__ import annotations

import hashlib
import math
from pathlib import Path

import numpy as np

from . import native_repair_mass as source
from . import native_code_transport as transport
from . import native_gauge_object as objects

ROOT = Path(__file__).resolve().parents[2]
TIME = 2.**-14
STEPS = 64
NUMERICAL_CONTRAST_ALLOWANCE = 1e-12
ARRAY_ADMISSION = 2e-10
ARRAY_SCHEMA = {'readout_terminal_blocks': ((3, 5, 6, 6), 'complex128'),
                'readout_probabilities': ((9, 5, 3), 'float64'),
                'readout_estimate': ((3, 3), 'float64')}


def fixed_inputs():
    """Input f is |Nambu0,f>; no action coefficients enter this basis."""
    result = np.zeros((3, 6, 6), complex)
    for f in range(3): result[f, f, f] = 1
    return result


def fixed_effects(f, g):
    """Pauli-Y plus/minus and the full four-dimensional complement."""
    if type(f) is not int or type(g) is not int or not 0 <= f < 3 or not 0 <= g < 3:
        raise ValueError('fixed family coordinate indices0,1,2 required')
    a, b = np.eye(6)[:, f], np.eye(6)[:, 3+g]
    plus = (a+1j*b)/math.sqrt(2)
    minus = (a-1j*b)/math.sqrt(2)
    return np.asarray([np.outer(plus, plus.conj()), np.outer(minus, minus.conj()),
                       np.eye(6)-np.outer(a, a)-np.outer(b, b)])


def error_budget():
    """Record-to-entry bound; the coherent parameter is not physical time.

    For ||Y||<=1, real symmetric Y, and H=X tensor Y, the fixed contrast is
    -2 cos(tY)[f,f] sin(tY)[g,f]. Its linear inverse differs from Y[g,f] by
    at most(2/3)t^2. The complete finite-channel diamond error contributes
    delta/(2t); the separately stated numerical contrast allowance does too.
    """
    diamond = source.averaging_bound(STEPS, TIME, initial_memory_stationary=True)
    finite_time = (2/3)*TIME*TIME
    channel = diamond/(2*TIME)
    numerical = NUMERICAL_CONTRAST_ALLOWANCE/(2*TIME)
    return {'time': TIME, 'steps': STEPS, 'kernel_operator_norm_upper': 1.,
            'stationary_product_initialization_required': True,
            'finite_channel_diamond_bound': diamond,
            'finite_time_entry_remainder': finite_time,
            'finite_channel_entry_allowance': channel,
            'numerical_contrast_allowance': NUMERICAL_CONTRAST_ALLOWANCE,
            'numerical_entry_allowance': numerical,
            'entry_error_upper': finite_time+channel+numerical,
            'numerical_allowance_is_outward_or_physical_certificate': False}


def decode_records(probabilities, *, time, entry_error_upper):
    """Blind inverse: terminal records and declared apparatus error only."""
    if (type(probabilities) is not np.ndarray or probabilities.shape != (9, 5, 3)
            or probabilities.dtype != np.dtype('float64') or not np.isfinite(probabilities).all()
            or np.min(probabilities) < -NUMERICAL_CONTRAST_ALLOWANCE):
        raise ValueError('all finite memory and plus/minus/complement record laws required')
    if (type(time) not in (int, float) or not math.isfinite(time) or time <= 0
            or type(entry_error_upper) not in (int, float) or not math.isfinite(entry_error_upper)
            or entry_error_upper < 0):
        raise ValueError('positive apparatus time and nonnegative entry error required')
    if np.max(abs(probabilities.sum(axis=(1, 2))-1)) > NUMERICAL_CONTRAST_ALLOWANCE:
        raise ValueError('terminal records must be complete without renormalization')
    inferred = np.empty((3, 3), float)
    for f in range(3):
        for g in range(3):
            law = probabilities[3*f+g].sum(axis=0)
            inferred[g, f] = -(law[0]-law[1])/(2*time)
    # Real symmetry is part of the declared generator class. Averaging the
    # two independent estimates does not enlarge their entrywise error.
    estimate = (inferred+inferred.T)/2
    radius = 3*entry_error_upper  # Frobenius norm bounds operator norm.
    values = np.linalg.eigvalsh(estimate)
    intervals = [[float(v-radius), float(v+radius)] for v in values]
    positive = intervals[0][0] > 0
    ratios = ([[1., 1.]]+[[lo/intervals[0][1], hi/intervals[0][0]]
                           for lo, hi in intervals[1:]]) if positive else None
    return {'estimate': estimate, 'entry_error_upper': entry_error_upper,
            'symmetric_estimator_disagreement_diagnostic': float(np.max(abs(inferred-inferred.T))),
            'kernel_operator_error_upper': radius, 'eigenvalue_intervals': intervals,
            'strictly_positive_kernel_resolved': positive, 'ordered_ratio_intervals': ratios,
            'particle_names_or_calibrated_mass_inputs': None}


def probabilities_from_blocks(blocks):
    if (type(blocks) is not np.ndarray or blocks.shape != (3, 5, 6, 6)
            or blocks.dtype != np.dtype('complex128') or not np.isfinite(blocks).all()):
        raise ValueError('all three inputs and five final quantum memory blocks required')
    result = np.empty((9, 5, 3), float)
    for f in range(3):
        for g in range(3):
            effects = fixed_effects(f, g)
            for k in range(5):
                law = np.einsum('ij,oji->o', blocks[f, k], effects)
                if np.max(abs(law.imag)) > NUMERICAL_CONTRAST_ALLOWANCE:
                    raise ValueError('terminal probability acquired an imaginary part')
                result[3*f+g, k] = law.real
    return result


def execute():
    """Run only finite correlated-memory source dynamics for the forward law."""
    stationary = np.asarray(source.architecture()['stationary'], float).reshape(5)
    terminal = np.asarray([source.evolve_cq(stationary[:, None, None]*rho, STEPS, TIME)
                           for rho in fixed_inputs()])
    probabilities = probabilities_from_blocks(terminal)
    bounds = error_budget()
    inferred = decode_records(probabilities, time=TIME, entry_error_upper=bounds['entry_error_upper'])
    arrays = {'readout_terminal_blocks': terminal, 'readout_probabilities': probabilities,
              'readout_estimate': inferred.pop('estimate')}
    report = {'fixed_computational_preparations': 3, 'detector_settings': 9,
              'settings_order': '3*f+g; input|0,f>, PauliY on|0,f> and|1,g>',
              'memory_outputs_per_setting': 5, 'pointer_outcomes_per_memory': 3,
              'pointer_order': ['plus', 'minus', 'orthogonal_complement'],
              'finite_feedback_steps_per_setting': STEPS,
              'one_execution_each_setting_feedback_steps': 9*STEPS,
              'initial_memory': 'exact stationary product pi tensor|0,f><0,f|',
              'error_budget': bounds, 'inferred_kernel': inferred,
              'scope': {'source_kernel_or_poles_in_preparation_or_detector': False,
                        'forward_evolution_uses_averaged_kernel_exponential': False,
                        'all_final_memory_and_pointer_outcomes_retained': True,
                        'ideal_Pauli_detector_and_stationary_preparation_supplied': True,
                        'analytic_record_laws_not_a_physical_shot_campaign': True,
                        'effective_kernel_eigenvalues_are_physical_pole_masses': False}}
    validate_arrays(arrays)
    return report, arrays


def validate_arrays(arrays):
    if set(arrays) != set(ARRAY_SCHEMA): raise ValueError('complete readout array menu required')
    for name, (shape, dtype) in ARRAY_SCHEMA.items():
        value = arrays[name]
        if (type(value) is not np.ndarray or value.shape != shape or value.dtype != np.dtype(dtype)
                or not np.isfinite(value).all()):
            raise ValueError('invalid readout array '+name)
    blocks = arrays['readout_terminal_blocks']
    for output in blocks:
        if np.max(abs(output-output.conj().transpose(0, 2, 1))) > NUMERICAL_CONTRAST_ALLOWANCE:
            raise ValueError('terminal memory blocks are not Hermitian')
        if min(np.linalg.eigvalsh(block).min() for block in output) < -NUMERICAL_CONTRAST_ALLOWANCE:
            raise ValueError('terminal memory block is not positive')
    law = probabilities_from_blocks(blocks)
    if np.max(abs(law-arrays['readout_probabilities'])) > ARRAY_ADMISSION:
        raise ValueError('records differ from retained final memory blocks')
    bounds = error_budget()
    inferred = decode_records(law, time=TIME, entry_error_upper=bounds['entry_error_upper'])
    if np.max(abs(inferred['estimate']-arrays['readout_estimate'])) > ARRAY_ADMISSION:
        raise ValueError('inferred kernel differs from retained records')


def report():
    return execute()[0]


def source_bindings():
    paths = set(transport.source_bindings()) | set(objects.source_bindings()) | {
        'model/common_source/native_repair_mass.py', 'model/common_source/test_native_repair_mass.py',
        'model/common_source/native_repair_mass_readout.py',
        'model/common_source/test_native_repair_mass_readout.py'}
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(paths)}
