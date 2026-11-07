"""Recompute the exact native tree/string sector and complete record receipt."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from ..native_quantum_gauss import packet as owning
from . import strings


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results.json'


def source_bindings():
    result = owning.source_bindings()
    for path in HERE.iterdir():
        if path.suffix in ('.py', '.md'):
            result[path.relative_to(ROOT).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
    return dict(sorted(result.items()))


def produce():
    bindings = source_bindings()
    tree = strings.Tree()
    result = {
        'schema': 'oph.native_gauge_strings.v1',
        'representations': strings.representation_certificate(),
        'same_owning_Hamiltonian': strings.one_edge_identity(),
        'complete_block_census': strings.central_sector_census(tree),
        'complete_tree_sector': strings.exact_tree_certificate(tree),
        'complete_position_records': [
            {'central_eigenvalue': charge, 'dimension': len(indices),
             'records': strings.complete_history(tree, charge)}
            for charge, indices in strings.native.native_blocks()],
        'continuous_Hamiltonian_volume_control': strings.tail_bound(12, 0, 3, 6, strings.s.Rational(1, 8)),
        'scope': {
            'all_one_fermion_Gauss_states_on_the_declared_tree_classified': True,
            'all_five_matter_blocks_and_sixty_link_states_in_the_ambient_carrier': True,
            'root_position_included': True,
            'same_existing_local_Hamiltonian_on_every_edge': True,
            'literal_raising_link_normalization_retained': True,
            'static_dual_reference_and_initial_root_singlet_are_inputs': True,
            'graph_coherent_coupling_electric_coefficient_and_clock_are_inputs': True,
            'same_CE_repair_law_or_compact_defect_band_identified': False,
            'cycles_or_additional_matter_included_in_kernel_classification': False,
            'relativistic_or_chiral_continuum_derived': False,
            'physical_or_asymptotic_particle_mass_identified': False,
            'detector_reservoir_or_reset_cost_derived': False},
        'source_bindings': bindings,
    }
    if bindings != source_bindings():
        raise ArithmeticError('source bindings changed during tree certificate generation')
    return result


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def verify(path=OUTPUT):
    if json.loads(Path(path).read_text()) != produce():
        raise ValueError('native gauge-string receipt differs from exact recomputation')
    return True


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_mutually_exclusive_group(required=True)
    actions.add_argument('--write', action='store_true')
    actions.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.write:
        OUTPUT.write_text(canonical(produce()))
        print(OUTPUT)
    else:
        verify()
        print('native gauge-string receipt OK')
