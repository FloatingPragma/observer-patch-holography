"""Native preparation of the fixed-number eta reference and local phase kick.

The eight-bit amplitude tree acts on the even physical modes while the
odd modes are initialized blank. Compact evaluation is an exact relabelling
of these spectator blanks, not a charge projection. The subsequent copy
and phase word is executed on all 65536 physical occupation rows.
"""
from functools import lru_cache
from fractions import Fraction as F
from types import MappingProxyType
import math

import numpy as np

from ..observer.source_scalar import compiler
from ..reference_recovery.recovery import CompleteRun, _failure_vectors, execute_complete
from ..reference_observation.observation import frozen
from ..certified_reference import certificate as cert
from . import cube, repair


ADMISSION = 3e-9


def reference_values():
    return [((-1)**sum(sum(cube.GEOMETRY.sites[j]) for j in range(8) if mask >> j & 1)
             if mask.bit_count() == 3 else 0) for mask in range(256)]


def tree_word():
    return compiler.compile_real_state(np.array(reference_values())/math.sqrt(56))


def copy_word():
    return [event for site in range(8) for event in compiler.cnot(2*site, 2*site+1)]


def phase_word():
    # exp(i pi/2 n_0)=exp(i pi/2) exp(-i pi/4 Z_0) exp(-i pi/4 Z_1).
    return [compiler.pulse('I', math.pi/2, [0]),
            compiler.pulse('Z', -math.pi/4, [0]),
            compiler.pulse('Z', -math.pi/4, [1])]


def complete16(word, initial):
    """Short declared-domain source word with all first-failure outputs."""
    if np.ma.isMaskedArray(initial): raise ValueError('masked initial states forbidden')
    a = np.asarray(initial)
    if a.ndim == 1: a = a[:, None]
    if (a.ndim != 2 or a.shape[0] != 65536 or not 1 <= a.shape[1] <= 4
            or a.dtype.kind not in 'fiuc' or not np.isfinite(a).all()):
        raise ValueError('all sixteen physical modes and one through four reference columns required')
    a = np.asarray(a, complex)
    validated = compiler._validate(word, 16)
    if len(validated) > 32: raise ValueError('the declared full-cube short source word has at most32 pulses')
    tape = tuple((e['kind'], e['angle'], tuple(e['targets'])) for e in validated)
    matrices = {}
    for kind, angle, _ in tape:
        if (kind, angle) not in matrices:
            full = compiler.native_pulse(kind, angle)
            if full.shape != (6, 6) or np.linalg.norm(full.conj().T@full-np.eye(6)) > 2e-11:
                raise ArithmeticError('literal full-source primitive lost unitarity admission')
            matrices[kind, angle] = frozen(full)
    state = a.copy(); weights = np.zeros((len(tape), 4)); gram = np.zeros((a.shape[1], a.shape[1]), complex)
    for event, (kind, angle, targets) in enumerate(tape):
        full = matrices[kind, angle]
        for slot, branch in enumerate(_failure_vectors(state, full, targets, 16)):
            gram += branch.conj().T@branch
            weights[event, slot] = np.vdot(branch, branch).real
        size = 1 << len(targets)
        state = compiler._apply(state, full[:size, :size], targets, 16)
    residual = np.linalg.norm(state.conj().T@state+gram-a.conj().T@a)
    if not math.isfinite(residual) or residual > ADMISSION*max(1., np.linalg.norm(a)**2):
        raise ArithmeticError('complete cube source output/reference accounting failed')
    return CompleteRun(tape, 16, frozen(a), MappingProxyType(matrices), frozen(state), frozen(weights), frozen(gram))


def expand_blanks(compact):
    """Exact spectator-blank isometry; all compact rows are retained."""
    a = np.asarray(compact)
    if a.ndim == 1: a = a[:, None]
    if a.ndim != 2 or a.shape[0] != 256: raise ValueError('complete eight-mode input required')
    out = np.zeros((65536, a.shape[1]), complex)
    for mask in range(256):
        out[sum(1 << (2*j) for j in range(8) if mask >> j & 1)] = a[mask]
    return out


def ideal_reference():
    out = np.zeros(65536, complex)
    for mask, value in enumerate(reference_values()):
        out[sum(3 << (2*j) for j in range(8) if mask >> j & 1)] = value/math.sqrt(56)
    return out


@lru_cache(maxsize=1)
def runs():
    initial = np.zeros(256); initial[0] = 1
    tree = execute_complete(tree_word(), 8, initial)
    copied = complete16(copy_word(), expand_blanks(tree.success))
    phased = complete16(phase_word(), copied.success)
    return tree, copied, phased


def record(name, run):
    return {'stage': name, 'pulses': len(run.word),
            'pair_code_event_indices': [j for j, (_, _, t) in enumerate(run.word) if len(t) == 2],
            'source_abort_slots': sum(6-(1 << len(t)) for _, _, t in run.word),
            'all_first_export_failure_weights': run.failure_weights.tolist(),
            'first_export_failure_probability': run.failure_probability}


def certificate():
    tree, copy, phase = tree_word(), copy_word(), phase_word()
    a = cert.amplitude_certificate(tree, reference_values(), bits=8)
    if a['offset'] != len(tree): raise ArithmeticError('unconsumed native amplitude-tree event')
    copy_template = [event for site in range(8) for event in cert.cnot(2*site, 2*site+1)]
    phase_template = [('I', (0,), cert._pi(1, 2)), ('Z', (0,), cert._pi(-1, 4)),
                      ('Z', (1,), cert._pi(-1, 4))]
    errors = []
    for word, template in ((copy, copy_template), (phase, phase_template)):
        offset, error, variables = cert.consume(word, 0, template)
        if offset != len(word) or variables: raise ArithmeticError('unconsumed short cube source event')
        errors.append(error)
    bound = a['tree_error']+a['fixed_angle_error']+sum(errors, F())+cert.controls_bound(tree+copy+phase, {})
    return {'signed_amplitude_tree_nodes_checked': a['checked_conditional_nodes'],
            'tree_pulses': len(tree), 'copy_pulses': len(copy), 'phase_pulses': len(phase),
            'compact_tree_word_sha256': cert.word_digest(tree),
            'physical_tree_targets': 'compact bit j means physical mode2j; all odd modes remain blank',
            'copy_word_sha256': cert.word_digest(copy), 'phase_word_sha256': cert.word_digest(phase),
            'complete_source_state_error_bound': str(bound),
            'complete_source_state_error_bound_upper': math.nextafter(float(bound), math.inf),
            'complete_65536_basis_state_word_executed': '24-pulse copy and3-pulse local phase',
            'source_abort_outputs_retained': True, 'charge_postselection_used': False,
            'native_field_wait_or_full_detector_executed_here': False,
            'physical_hardware_precision_proved': False}


def report():
    tree, copied, phased = runs()
    target = ideal_reference()
    occupancy = np.array([(mask & 3).bit_count() for mask in range(65536)])
    kicked = np.exp(1j*math.pi*occupancy/2)*target
    error = math.sqrt(float(np.linalg.norm(phased.success[:, 0]-kicked)**2)+
                      sum(run.failure_probability for run in (tree, copied, phased)))
    continuation = phased.success
    records = []
    for law in ('C', 'E'):
        for seam in repair.SEAMS:
            outcomes = repair.apply(continuation, law, seam)
            for label, state in outcomes:
                records.append({'law': law, 'seam': [repair.PORTS[j] for j in seam],
                                'record': list(label), 'probability': float(np.vdot(state, state).real)})
            continuation = outcomes[-1][1]
    native_mass = math.fsum(row['probability'] for row in records if row['record'][0] == 'native')
    total = float(np.vdot(continuation, continuation).real)+native_mass+sum(run.failure_probability for run in (tree, copied, phased))
    if error > ADMISSION or abs(total-1) > ADMISSION:
        raise ArithmeticError('native cube preparation/phase differs from complete ideal output')
    return {'source_records': [record(name, run) for name, run in
                               zip(('signed_tree', 'pair_copy', 'local_phase'), (tree, copied, phased), strict=True)],
            'native_preparation_error': float(np.linalg.norm(copied.success[:, 0]-target)),
            'complete_preparation_and_phase_error': error,
            'all130_following_repair_records': records,
            'following_repair_native_probability': native_mass,
            'complete_output_probability': total, 'host_numerical_admission': ADMISSION,
            'native_source_preparation_and_local_phase_executed': True,
            'spectrum_detection_is_analytic_not_empirically_resolved': True}
