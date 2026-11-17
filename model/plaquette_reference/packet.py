"""Retained same-action plaquette preparation and regional modular comparison."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from tools import strict_json
from ..field_thermal.thermal import car_restriction, restrict
from ..full_fock_field import field
from ..full_fock_field.local_simulation import split_bound
from ..observer.scalar_response.experiment import canonical, compare_receipt
from ..reference_recovery.recovery import execute_complete
from ..regional_modular_expansion import expansion
from . import preparation as p


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
SOURCE_PATHS = (
    'model/plaquette_reference/preparation.py',
    'model/plaquette_reference/packet.py',
    'model/plaquette_reference/test_preparation.py',
    'model/plaquette_reference/README.md',
    'model/plaquette_reference/contract.json',
    'model/observer/source_scalar/compiler.py',
    'model/reference_recovery/recovery.py',
    'model/field_thermal/thermal.py',
    'model/field_symmetry/symmetry.py',
    'model/field_stress/energy.py',
    'model/full_fock_field/field.py',
    'model/full_fock_field/local_simulation.py',
    'model/regional_modular_expansion/expansion.py',
    'code/source_selection_model/response.json',
    'code/m1_source_realization/receipt.json',
)


def negative_log(matrix):
    values, vectors = np.linalg.eigh(matrix)
    if not np.isfinite(values).all() or values[0] <= 0:
        raise ArithmeticError('the complete reference must remain faithful')
    return -(vectors*np.log(values))@vectors.conj().T, float(values[0])


def resource_report(prepared):
    runs = (prepared.probability_run,)+prepared.conditional_runs
    counts = [len(run.word) for run in runs]
    pairs = [sum(kind == 'J' for kind, _, _ in run.word) for run in runs]
    failures = sum(4*(n-j)+2*j for n, j in zip(counts, pairs, strict=True))
    labels = p.DIMENSION+failures+4*prepared.readout.events
    return {
        'field_work_bits': p.MODES,
        'initialized_field_bits': p.MODES,
        'retained_classical_record_bits': p.MODES,
        'all_occupation_inputs_and_records': p.DIMENSION,
        'probability_preparation_pulses': counts[0],
        'conditional_pulses_min': min(counts[1:]),
        'conditional_pulses_max': max(counts[1:]),
        'trajectory_pulses_max': counts[0]+max(counts[1:]),
        'all_enumerated_native_pulses': sum(counts),
        'all_enumerated_load_and_export_events': 2*sum(counts),
        'first_native_export_failure_Kraus_slots': failures,
        'success_and_first_failure_labels_enumerated': labels,
        'bits_to_encode_a_complete_output_label': (labels-1).bit_length(),
        'abort_spectator_quantum_bits_at_most': p.MODES-1,
        'readout_adaptive_events_enumerated': prepared.readout.events,
        'readout_load_and_export_events_per_trajectory': 2*p.MODES,
        'readout_first_failure_Kraus_slots': 4*prepared.readout.events,
        'every_zero_angle_native_request_retained': True,
        'record_conditioned_blanking_is_reversible_and_source_compiled': True,
        'no_postselection_or_success_renormalization': True,
        'physical_initialization_readout_precision_and_storage_costs_calibrated': False,
    }


def word_digest(run):
    """Bind the immutable, actually executed native tape, including zeros."""
    digest = hashlib.sha256()
    for kind, angle, targets in run.word:
        digest.update((kind+'|'+angle.hex()+'|'+','.join(map(str, targets))+';').encode())
    return digest.hexdigest()


def regional_comparison(prepared):
    data = p.design()
    theorem = p.error_theorem()
    h = data.hamiltonian
    raw = expm(-float(p.BETA)*(h-np.eye(p.DIMENSION)))
    gibbs = raw/np.trace(raw)
    kraus = car_restriction(p.MODES, expansion.KEPT)
    actual = restrict(prepared.field_density, kraus)
    target = restrict(gibbs, kraus)
    ideal_columns = data.target_columns()
    ideal = restrict(ideal_columns@ideal_columns.conj().T, kraus)
    k_actual, minimum_actual = negative_log(actual)
    k_target, _ = negative_log(target)
    polynomial = np.asarray(expansion.polynomial(expansion.plaquette(), p.BETA), complex)
    # E=Tr_B/16; the exact Gibbs normalization scalar remains in K_A.
    normalization_scalar = float(np.log(np.trace(restrict(raw, kraus)/16).real))
    normalized_polynomial = polynomial+normalization_scalar*np.eye(16)
    tail = expansion.logarithm_tail(p.BETA, F(3))
    admission = p.EXECUTION_ADMISSION
    admitted_density = (2+admission)*admission
    admitted_minimum = theorem['regional_minimum_weight']-admitted_density
    if admitted_minimum <= 0:
        raise ArithmeticError('the declared complete-output admission does not retain faithfulness')
    admitted_bridge = (theorem['density_trace_error']+admitted_density)/admitted_minimum
    # (2+eps)*eps also covers the small floating departure from unit trace;
    # the exact ideal complete instrument is trace preserving.
    execution_density = (2+prepared.complete_output_error)*prepared.complete_output_error
    total_density = float(theorem['density_trace_error'])+execution_density
    minimum_bound = float(theorem['regional_minimum_weight'])-execution_density
    numerical_bridge = total_density/minimum_bound
    combined_bound = float(tail)+numerical_bridge
    residual = float(np.linalg.norm(k_actual-normalized_polynomial, 2))
    if residual > combined_bound+1e-13 or minimum_actual < minimum_bound-1e-13:
        raise ArithmeticError('prepared reference lost the regional modular error budget')
    centered = lambda a: a-np.trace(a)/16*np.eye(16)
    return {
        'source_success_trace_without_renormalization': float(np.trace(actual).real),
        'source_regional_minimum_weight': minimum_actual,
        'ideal_polynomial_vs_Gibbs_density_trace_diagnostic': float(np.linalg.norm(ideal-target, 'nuc')),
        'executed_vs_Gibbs_regional_modular_operator_diagnostic': float(np.linalg.norm(k_actual-k_target, 2)),
        'exact_Gibbs_normalization_scalar_retained': normalization_scalar,
        'executed_normalized_P8_operator_error': residual,
        'executed_trace_zero_tangent_operator_error': float(np.linalg.norm(centered(k_actual-polynomial), 2)),
        'exact_P8_remainder': str(tail),
        'conditional_admission_density_error_exact': str(admitted_density),
        'conditional_admission_regional_minimum_exact': str(admitted_minimum),
        'conditional_admission_modular_bridge_exact': str(admitted_bridge),
        'conditional_admission_combined_P8_bound_exact': str(tail+admitted_bridge),
        'conditional_bound_premise': 'complete-output error <= 2e-11; this is an admitted implementation tolerance, not a theorem of hardware precision',
        'execution_diagnostic_bridge_bound': numerical_bridge,
        'combined_P8_and_execution_diagnostic_bound': combined_bound,
        'all_Hermitian_trace_zero_tangents': '|DS(delta)-Tr(delta P8)| <= combined_bound ||delta||_1',
        'floating_source_execution_is_an_admission_check_not_an_interval_proof': True,
        'exact_ideal_source_word_theorem_has_only_the_rational_design_error': True,
        'executed_regional_object': 'unnormalized success field block; every abort complement is a separate retained output',
        'ideal_reference_is_the_unconditioned_field_marginal': True,
    }


def continuation(prepared):
    """Same original action, source-executed on every retained ensemble record."""
    parameter, steps = F(1, 2), 2
    duration = -p.BETA*parameter
    word = field.source_word(p.GEOMETRY, float(duration), steps)
    runs = tuple(execute_complete(word, p.MODES, prepared.columns[:, start:start+16])
                 for start in range(0, p.DIMENSION, 16))
    actual = np.column_stack([run.success for run in runs])
    target = expm(-1j*float(duration)*p.design().hamiltonian)@prepared.columns
    failures = math.fsum(run.failure_probability for run in runs)
    error = math.sqrt(float(np.linalg.norm(actual-target))**2+failures)
    bound = split_bound(p.GEOMETRY.sides, abs(duration), steps)
    if error > float(bound)+1e-10:
        raise ArithmeticError('continued prepared reference exceeds the source action splitting bound')
    return {
        'modular_parameter': str(parameter), 'Hamiltonian_model_duration': str(duration),
        'source_pulses_per_record': len(word), 'retained_classical_records': p.DIMENSION,
        'all_record_columns_executed_in_batches_of': 16,
        'split_operator_error_bound_exact': str(bound),
        'ideal_prepared_modular_generator_error_up_to_scalar_exact': str(p.error_theorem()['full_modular_operator_error']),
        'ideal_prepared_modular_word_error_up_to_scalar_exact': str(
            bound+abs(parameter)*p.error_theorem()['full_modular_operator_error']),
        'complete_output_error_diagnostic': error,
        'first_export_failure_probability': failures,
        'field_stationarity_error_diagnostic': float(np.linalg.norm(actual@actual.conj().T-prepared.field_density, 'nuc')),
        'ideal_polynomial_reference_commutes_with_original_H_exactly': True,
        'all_native_exports_and_preparation_records_retained': True,
        'physical_clock_calibration_established': False,
    }


def run():
    prepared = p.prepare()
    theorem = p.error_theorem()
    exact = {key: str(value) for key, value in theorem.items()}
    exact['ideal_P8_entropy_tangent_bound'] = str(
        expansion.logarithm_tail(p.BETA, F(3))+theorem['regional_modular_operator_error'])
    return {
        'schema': 'oph.plaquette_reference.v1', 'issues': [1138, 1145],
        'gates': ['CP-02', 'CP-09'], 'all_pass': True,
        'domain': {
            'spatial_sites': [2, 2, 1], 'full_CAR_modes': p.MODES,
            'full_Fock_dimension': p.DIMENSION,
            'number_sector_dimensions': [1, 8, 28, 56, 70, 56, 28, 8, 1],
            'speed': '1/8', 'interaction': '1', 'beta': str(p.BETA),
            'action': 'original full-Fock Weyl/contact H; scalar centering and exact onsite frame only',
            'reference': 'normalized Q8^2, Q8=sum_{k=0}^8[-beta(H-I)/2]^k/k!, with rigorous error to Gibbs',
            'extension': 'all 256 classical column records retained; not a pure thermofield purification',
        },
        'exact_design': p.design().certificate(), 'exact_error_theorem': exact,
        'resources': resource_report(prepared),
        'execution': {
            'complete_output_error_diagnostic': prepared.complete_output_error,
            'numerical_admission': str(p.EXECUTION_ADMISSION),
            'first_export_failure_probability': prepared.failure_probability,
            'all_branch_probabilities': [float(np.vdot(x, x).real) for x in prepared.columns.T],
            'total_success_plus_abort': float(np.vdot(prepared.columns, prepared.columns).real+prepared.failure_probability),
            'readout_probability_is_source_prepared': True,
            'exact_valid_code_export_failure_probability': 0,
            'probability_native_word_sha256': word_digest(prepared.probability_run),
            'conditional_native_word_sha256': [word_digest(run) for run in prepared.conditional_runs],
        },
        'regional_modular_comparison': regional_comparison(prepared),
        'continued_modular_action': continuation(prepared),
        'scope': {
            'closed': 'a complete same-action finite source ensemble reference, with quantitative Gibbs and regional-modular errors and continued original action',
            'retained_inputs': ['declared action and temperature', 'initialized work bits',
                                'coherent load/export', 'complete source readout and conditional scheduler',
                                'analogue angle precision and memory services'],
            'remaining': ['thermodynamic/scaling or interacting continuum theorem',
                          'physical clock, cap geometry and Einstein normalization',
                          'a pure full-plaquette thermofield source word, if separately required'],
            'no_physical_particle_or_Einstein_gate_accepted': True,
            'classical_record_conditioning_is_not_a_claim_that_each_conditional_field_is_Gibbs': True,
        },
        'source_sha256': {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in SOURCE_PATHS},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--write', action='store_true')
    modes.add_argument('--check', nargs='?', const=OUTPUT, type=Path)
    args = parser.parse_args()
    result = run()
    if args.write:
        OUTPUT.write_text(canonical(result))
    if args.check:
        compare_receipt(strict_json.loads(args.check.read_text()), result)
    summary = {key: result[key] for key in ('all_pass', 'resources', 'regional_modular_comparison', 'continued_modular_action')}
    summary['execution'] = {key: result['execution'][key] for key in (
        'complete_output_error_diagnostic', 'numerical_admission',
        'first_export_failure_probability', 'total_success_plus_abort')}
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
