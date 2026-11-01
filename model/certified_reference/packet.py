"""Replay the exact full-source error and regional modular certificate."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from tools import strict_json
from ..observer.scalar_response.experiment import canonical, compare_receipt
from ..plaquette_reference import preparation as p
from ..regional_modular_expansion import expansion
from . import certificate, continuation, design_check


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
SOURCE_PATHS = (
    'model/certified_reference/interval.py', 'model/certified_reference/source.py',
    'model/certified_reference/certificate.py', 'model/certified_reference/continuation.py',
    'model/certified_reference/design_check.py', 'model/certified_reference/packet.py',
    'model/certified_reference/test_certificate.py', 'model/certified_reference/README.md',
    'model/certified_reference/contract.json',
    'model/plaquette_reference/preparation.py', 'model/plaquette_reference/results.json',
    'model/observer/source_scalar/compiler.py', 'model/reference_recovery/recovery.py',
    'model/field_symmetry/symmetry.py', 'model/field_stress/energy.py',
    'model/full_fock_field/field.py', 'model/full_fock_field/local_simulation.py',
    'model/common_source/chiral_candidate/gaussian_source.py',
    'model/interacting_field/source.py', 'model/regional_modular_expansion/expansion.py',
    'code/source_selection_model/response.json', 'code/m1_source_realization/receipt.json',
)


def encode(value):
    if type(value) is F:
        return str(value)
    if isinstance(value, dict):
        return {k: encode(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(v) for v in value]
    return value


def modular_bounds(error, continuation_error):
    if type(error) is not F or not 0 <= error < F(1, 10**6):
        raise ValueError('proved small rational complete-output error required')
    theorem = p.error_theorem()
    execution_density = (2+error)*error
    minimum = theorem['regional_minimum_weight']-execution_density
    full_minimum = theorem['full_minimum_weight']-execution_density
    if min(minimum, full_minimum) <= 0:
        raise ArithmeticError('literal-source certificate did not retain faithfulness')
    density = theorem['density_trace_error']+execution_density
    bridge = density/minimum
    tail = expansion.logarithm_tail(p.BETA, F(3))
    return {'complete_output_vector_error': error,
            'all_abort_probability_upper': error**2,
            'source_success_density_trace_error': execution_density,
            'source_plus_design_Gibbs_trace_error': density,
            'actual_regional_success_minimum_weight': minimum,
            'actual_full_success_minimum_weight': full_minimum,
            'regional_modular_operator_error': bridge,
            'regional_P8_remainder': tail,
            'regional_P8_total_operator_and_trace_zero_tangent_bound': tail+bridge,
            'full_modular_generator_error': density/full_minimum,
            'continued_actual_success_modular_word_error_up_to_Gibbs_scalar':
                continuation_error+F(1, 2)*density/full_minimum,
            'regional_log_is_of_unnormalized_success_block': True,
            'Gibbs_normalization_scalar_retained_in_absolute_P8_comparison': True,
            'all_abort_outputs_remain_separate_without_success_renormalization': True}


def run():
    design = design_check.certify()
    prepared = certificate.certify_preparation()
    continued = continuation.certify()
    prior = strict_json.loads((ROOT/'model/plaquette_reference/results.json').read_text())
    if prepared['probability_native_word_sha256'] != prior['execution']['probability_native_word_sha256']:
        raise ArithmeticError('this is not the previously executed probability word')
    if [r['native_word_sha256'] for r in prepared['conditional_records']] != prior['execution']['conditional_native_word_sha256']:
        raise ArithmeticError('a retained conditional preparation word changed')
    bounds = modular_bounds(prepared['complete_preparation_output_error'],
                            continued['full_action_complete_operator_error'])
    return encode({
        'schema': 'oph.certified_reference.v1', 'issues': [1138, 1145],
        'gates': ['CP-02', 'CP-09'], 'all_pass': True,
        'domain': {'geometry': [2, 2, 1], 'field_modes': 8, 'full_Fock_dimension': 256,
                   'beta': F(1, 32), 'speed': F(1, 8), 'interaction': F(1),
                   'prepared_target': 'Q8^2/Tr(Q8^2), Q8=sum_{k=0}^8[-(H-I)/64]^k/k!',
                   'source_law': 'exact retained Q(sqrt(5)) M6 response with literal binary64 integrated port values',
                   'all_256_classical_records_and_every_export_abort_retained': True},
        'independent_design_binding': design, 'preparation_certificate': prepared,
        'continuation_certificate': continued, 'modular_certificate': bounds,
        'previous_complete_output_admission': p.EXECUTION_ADMISSION,
        'previous_admission_now_proved_for_literal_source_specification':
            prepared['complete_preparation_output_error'] < p.EXECUTION_ADMISSION,
        'numerical_expm_residual_used_in_proof': False,
        'source_target_angle_evaluation_trusted': False,
        'outward_kernel': {'dyadic_precision_bits': 160, 'trigonometric_Taylor_terms': 64,
                           'pi_Machin_atan_terms': 80, 'integer_square_roots': True},
        'scope': {
            'closed': 'literal finite source-control preparation and modular error, including represented-angle and represented-port rounding; former complete-output tolerance is no longer an unproved premise for this mathematical source',
            'declared_services': ['exact retained six-level response law', 'coherent proper-code load/export',
                                  'complete readout and conditional scheduler', 'initialized work register',
                                  'memory and all retained record/abort outputs'],
            'not_claimed': ['correct rounding of host expm or all floating simulation matrices',
                            'hardware agreement with the declared response and service laws',
                            'scalable local/thermodynamic Gibbs preparation',
                            'physical clock calibration, continuum fields or Einstein dynamics'],
            'no_physical_gate_or_canonical_premise_promoted': True,
        },
        'source_sha256': {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in SOURCE_PATHS},
    })


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
    print(json.dumps({key: result[key] for key in ('all_pass', 'modular_certificate',
                     'previous_admission_now_proved_for_literal_source_specification', 'scope')}, indent=2))


if __name__ == '__main__':
    main()
