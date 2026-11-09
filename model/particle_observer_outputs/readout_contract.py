"""Complete same-source finite instruments and inherited compiler budgets.

The finite Hamiltonian, preparation service and operational tick are declared
model data. This module composes their observer outputs without selecting an
architecture or identifying any finite energy with a laboratory particle.
"""
from __future__ import annotations

from fractions import Fraction as F

import sympy as s

from model.spin_operational_position.bounds import generic_instrument_plan


def _zero(matrix):
    return all(s.simplify(value) == 0 for value in matrix.todok().values())


def verify_projectors(projectors):
    """All-input PVM identities, including coherent degenerate eigenspaces."""
    projectors = dict(projectors)
    if not projectors:
        raise ValueError('at least one retained outcome is required')
    size = next(iter(projectors.values())).rows
    if size < 1:
        raise ValueError('a nonempty finite field is required')
    for label, value in projectors.items():
        if not isinstance(label, str) or not label or value.shape != (size, size):
            raise ValueError('labeled square projectors on one field space are required')
        if any(v.has(s.Float) for v in value.todok().values()):
            raise ValueError('exact operator entries are required')
        if not _zero(value-value.H) or not _zero(value*value-value):
            raise ValueError('every readout branch must be its complete orthogonal projection')
    values = list(projectors.values())
    if not _zero(sum(values, s.zeros(size))-s.eye(size)):
        raise ValueError('the readout omits an input subspace')
    for i, left in enumerate(values):
        for right in values[i+1:]:
            if not _zero(left*right):
                raise ValueError('projective readout outcomes overlap')
    return {'field_dimension': size,
            'retained_outcome_ranks': {label: int(s.trace(p)) for label, p in projectors.items()},
            'all_input_columns_included': True,
            'whole_degenerate_subspaces_preserved': True,
            'spectral_or_physical_interpretation_supplied_by_owner': True}


def energy_menu(data):
    """Use the exact owning Higgs Hamiltonian and its three spectral blocks."""
    menu = {'lower_energy': data['ground_projector'],
            'upper_energy': data['excited_projector'],
            'dark_energy': data['dark_projector']}
    verify_projectors(menu)
    energies = {'lower_energy': data['ground_energy'],
                'upper_energy': data['excited_energy'],
                'dark_energy': data['dark_energy']}
    h = data['H']
    if any(not _zero(h*p-energies[label]*p) for label, p in menu.items()):
        raise ValueError('spectral readout does not belong to this Hamiltonian')
    if not _zero(h-sum((energies[label]*p for label, p in menu.items()), s.zeros(h.rows))):
        raise ValueError('the complete spectral reconstruction differs from the source')
    return menu, energies


def scalar_number_menu(data):
    """The two complete scalar-number outcomes on the same conserved sector."""
    empty = data['vacuum_projector']
    menu = {'empty_scalar': empty, 'occupied_scalar': s.eye(empty.rows)-empty}
    verify_projectors(menu)
    return menu


def pad_menu(projectors):
    """Encode all field inputs and retain every unused binary input as outside."""
    verified = verify_projectors(projectors)
    size = verified['field_dimension']
    bits = max(1, (size-1).bit_length())
    padded_size = 1 << bits
    result = {label: s.ImmutableSparseMatrix(padded_size, padded_size, p.todok())
              for label, p in projectors.items()}
    if padded_size != size:
        result['outside_field_encoding'] = s.ImmutableSparseMatrix(
            padded_size, padded_size, {(j, j): 1 for j in range(size, padded_size)})
    verify_projectors(result)
    return result


def pointer_blocks(projectors):
    """Exact unitary controlled addition, represented by pointer matrix blocks.

    U=sum_r Pi_r tensor shift^r, with zero padded record labels. The returned
    blocks use pointer-major ordering, the corresponding tensor permutation.
    Acting on blank pointer0 gives Kraus Pi_r. The algebra proves unitarity
    without allocating a large dense apparatus matrix.
    """
    verified = verify_projectors(projectors)
    values = tuple(projectors.values())
    size = verified['field_dimension']
    pointer_size = 1 << (len(values)-1).bit_length()
    zero = s.ImmutableSparseMatrix(size, size, {})
    padded = values+(zero,)*(pointer_size-len(values))
    blocks = tuple(tuple(padded[(a-b) % pointer_size] for b in range(pointer_size))
                   for a in range(pointer_size))
    return blocks


def compiled_menu_plan(projectors, error=F(1, 10**6)):
    """Consume the existing finite-instrument theorem for this exact menu."""
    padded = pad_menu(projectors)
    size = next(iter(padded.values())).rows
    plan = generic_instrument_plan((size-1).bit_length(), len(padded), error)
    plan['field_dimension_before_encoding'] = next(iter(projectors.values())).rows
    plan['retained_outcome_ranks'] = {label: int(s.trace(p)) for label, p in padded.items()}
    plan['unused_pointer_labels_retained_as_zero'] = (1 << plan['coherent_record_bits'])-len(padded)
    plan['source_primitive_route_also_admissible'] = True
    plan['compiled_tape_executed_by_this_packet'] = False
    return plan


def composition_error(*, ticks, detector_uses, preparation_error=F(0),
                      step_error=F(0), detector_error=F(0), initial_state_error=F(0)):
    """Half-trace/TV bound for complete adaptive histories, with references.

    One identical finite macro-step is repeated at every counted tick.
    Preparation/step/detector errors are complete half-diamond bounds; the
    initial-state term is a complete half trace-distance bound.
    """
    if type(ticks) is not int or ticks < 0 or type(detector_uses) is not int or detector_uses < 0:
        raise ValueError('nonnegative exact use counts are required')
    errors = (preparation_error, step_error, detector_error, initial_state_error)
    if any(type(value) is not F or not 0 <= value <= 1 for value in errors):
        raise ValueError('exact complete-channel error bounds in [0,1] are required')
    return min(F(1), initial_state_error+preparation_error+ticks*step_error+detector_uses*detector_error)


def macro_clock_plan(field_dimension=307, error=F(1, 10**6)):
    """One fixed compiled unitary word repeated; no physical-duration claim."""
    if type(field_dimension) is not int or field_dimension < 2:
        raise ValueError('a finite field dimension at least two is required')
    bits = (field_dimension-1).bit_length()
    plan = generic_instrument_plan(bits, 1, error)
    return {'declared_fixture_tick': 'delta_t=pi/45 in the chosen Hamiltonian unit',
            'n_tick_ideal_operation': '(exp(-i delta_t H))^n=exp(-i n delta_t H)',
            'identical_compiled_word_repeated': True,
            'unitary_compile_plan': plan,
            'primitive_clock_means_complete_macro_word_count': True,
            'preparation_and_detection_are_separate_declared_services': True,
            'autonomous_H_paused_or_refocused_during_these_services': True,
            'first_failure_and_all_pointer_records_retained': True,
            'physical_seconds_or_GeV_assigned': False}


def certificate():
    # The owning module certifies its construction directly against the full
    # finite-Higgs vertex. No separate substitute Hamiltonian is introduced.
    from . import higgs_records
    data = higgs_records.sector_data(s.diag(1, 2, 3), 1)
    energy, values = energy_menu(data)
    number = scalar_number_menu(data)
    if data['H'].rows != 307 or data['gap'] != 15:
        raise ArithmeticError('the shared exact Higgs fixture changed')
    unitary_plan = macro_clock_plan()
    prepare = higgs_records.preparation_unitary(data)
    if not _zero(prepare.H*prepare-s.eye(307)) or not _zero(
            prepare*data['vacuum_projector']*prepare.H-data['ground_projector']):
        raise ArithmeticError('the common-source preparation fails its complete unitary contract')
    return {'schema': 'oph.finite_observer_readout_contract.v1',
            'source_fixture': {'Y': 'diag(1,2,3)', 'pair_cost': '1', 'scalar_tags': [0, 1, 2],
                               'conserved_sector': 'N_F-2N_H=0', 'dimension': 307},
            'spectral_instrument': verify_projectors(energy),
            'spectral_values': {label: str(value) for label, value in values.items()},
            'ground_referenced_spectral_values': {
                label: str(value-data['ground_energy']) for label, value in values.items()},
            'scalar_number_instrument': verify_projectors(number),
            'same_source_ground_preparation': {
                'service': 'phase-shifted finite-Higgs vertex pulse with declared number-term pause',
                'all_307_input_columns_unitary': True,
                'empty_input_is_sent_to_exact_ground': True,
                'same_inherited_unitary_compile_plan_as_one_tick': True},
            'spectral_compiler': compiled_menu_plan(energy),
            'number_compiler': compiled_menu_plan(number),
            'macro_clock': unitary_plan,
            'two_tick_two_detector_complete_error': str(composition_error(
                ticks=2, detector_uses=2, preparation_error=F(1, 10**6),
                step_error=F(1, 10**6), detector_error=F(1, 10**6))),
            'original_full_CAR_scalar_dimension': 3*(1 << 18),
            'outside_conserved_sector_dimension': 3*(1 << 18)-307,
            'outside_sector_route': 'retain I-P_sector as its own outcome when the ambient full field is used',
            'entire_ambient_energy_spectrum_computed': False,
            'quantum_semantics_and_coherent_control_interface_declared': True,
            'architecture_uniqueness_required': False,
            'bare_occupation_repair_alone_compiles_detector': False,
            'native_compiled_tape_execution_claimed': False,
            'laboratory_particle_identity_or_mass_assigned': False}
