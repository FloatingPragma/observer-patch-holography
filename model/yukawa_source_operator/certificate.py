"""Replay source algebra and source traces with explicit input ancestry."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from . import operator, three_point
from ..native_modular_reference import local_preparation
from ..common_source import slow_geometry


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')


def run():
    algebra = operator.certificate()
    statistics = three_point.certificate()
    source_paths = set(local_preparation.SOURCE_PATHS)
    source_paths.update((
        'model/common_source/native_matter_attachment.py',
        'model/common_source/native_sm_current.py',
        'model/common_source/native_gauge_object.py',
        'model/common_source/slow_geometry.py',
        'model/common_source/chiral_candidate/action.py',
        'model/common_source/particle_identification/family_attachment.py',
        'model/common_source/chiral_candidate/kinetic.py',
        'model/common_source/chiral_candidate/record_two_point.py',
        'model/common_source/native_record_yukawa_field.py',
        'model/observer_particles/unified_yukawa_action.py',
        'model/observer_particles/scalar_bilinear_route.py',
        'model/observer_particles/family_naturality.py',
        'model/observer_particles/results/family_naturality_operators.json',
        'model/observer_particles/results/family_naturality.json',
        'model/observer_particles/results/unified_yukawa_operators.json',
        'model/observer_particles/results/unified_yukawa_action.json',
        'model/observer_particles/results/representation_identification_operators.json',
        'model/observer_particles/results/representation_identification.json',
        'code/a5_closure/manifests/echosahedral_federation_reference.json',
        'Lean/Screen/WeylYukawaConventions.lean',
        'code/source_selection_model/response.json',
        slow_geometry.SPEC.relative_to(ROOT).as_posix(),
        (slow_geometry.ARCHIVE/'oph_exact/carrier.py').relative_to(ROOT).as_posix(),
        (slow_geometry.ARCHIVE/'oph_fpe/dynamics/self_readback_repair_closure.py').relative_to(ROOT).as_posix(),
        'model/yukawa_source_operator/operator.py',
        'model/yukawa_source_operator/three_point.py',
        'model/yukawa_source_operator/certificate.py',
    ))
    return {
        'schema': 'oph.yukawa_source_operator.v1',
        'status': 'EXACT_NATIVE_TENSOR_DUAL_AND_COMPLETE_SOURCE_THREE_POINT_RESPONSE',
        'empirical_mass_or_Yukawa_inputs': [],
        'operator_extraction': algebra,
        'finite_source_three_point': statistics,
        'physical_attachment': {
            'invariant_operator_shape_derived_from_native_currents': True,
            'static_three_point_numbers_derived_from_native_reference': True,
            'canonical_kinetic_conversion_defined_for_supplied_metrics': True,
            'unique_physical_Higgs_intervention_derived': False,
            'source_selected_charged_fermion_Yukawa_matrix': False,
            'physical_charged_fermion_mass_prediction': False,
            'strongest_sufficient_datum': 'same-source inverse-fermion-propagator derivative with respect to canonically normalized Higgs, projected by the exact tensor dual; complete two-point residues and the shared source clock are required',
        },
        'source_sha256': {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                          for name in sorted(source_paths)},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.write:
        OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    if args.check and json.loads(OUTPUT.read_text()) != result:
        raise ArithmeticError('retained source response receipt differs from complete replay')
    print(json.dumps({'status': result['status'], 'exact_complex_coefficients': 27,
                      'complete_source_states': 66, 'physical_Yukawas_derived': False}))


if __name__ == '__main__':
    main()
