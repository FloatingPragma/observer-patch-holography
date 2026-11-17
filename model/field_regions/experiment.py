"""The same interacting CAR history read on nested spatial regions."""
from __future__ import annotations
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy import sparse
from ..interacting_field.lattice import TwoFermions, WeylLattice
from ..interacting_field import source, bounds
from ..interacting_field.experiment import protocol, pure_cq_distance, apparatus_resources
from .regions import (
    regional_basis, restrict_pure, restrict_regional, lift_number_observable,
    conditional_expectation, cube_mode_map, two_particle_embedding,
)

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
TIME = F(1, 16)
STEPS = 8


def density_cq_distance(first, second):
    return sum(float(np.linalg.eigvalsh(a-b).__abs__().sum())/2
               for a, b in zip(first, second, strict=True))


def sector_weights(density, kept, reference_dimension=1):
    basis = regional_basis(kept)
    return [float(sum(density[j*reference_dimension+r, j*reference_dimension+r].real
                      for j, occupied in enumerate(basis) if len(occupied) == n
                      for r in range(reference_dimension))) for n in range(3)]


def ramsey_effect(kept, p, q):
    """Local full-Fock orbital occupation restricted to N_local<=2."""
    m = len(kept)
    h = np.zeros((m, m), complex)
    ip, iq = kept.index(p), kept.index(q)
    h[ip, ip] = h[iq, iq] = .5
    h[ip, iq] = -.5j; h[iq, ip] = .5j
    result = np.zeros((1+m+m*(m-1)//2,)*2, complex)
    result[1:1+m, 1:1+m] = h
    result[1+m:, 1+m:] = TwoFermions(m).lift_one_body(h).toarray()
    return result


def history_level(radius):
    lattice = WeylLattice(radius)
    probabilities, outputs = protocol(lattice, steps=STEPS)
    exact_probabilities, exact_outputs = protocol(lattice, exact=True)
    kept = (lattice.origin, lattice.origin+1, lattice.neighbor, lattice.neighbor+1)
    reduced = [restrict_pure(lattice, state, kept) for state in outputs]
    exact_reduced = [restrict_pure(lattice, state, kept) for state in exact_outputs]
    single_site = kept[:2]
    nested_error = max(float(np.linalg.norm(restrict_regional(state, kept, single_site)
                    - restrict_pure(lattice, vector, single_site), 2))
                    for state, vector in zip(reduced, outputs, strict=True))
    record_error = max(abs(float(np.trace(state).real)-probability)
                       for state, probability in zip(reduced, probabilities, strict=True))
    temporal = F(bounds.history(((TIME, STEPS), (TIME, STEPS)))['complete_history_half_trace_upper'])
    spatial = F(bounds.boundary(radius, 2*TIME)['complete_history_half_trace_upper'])
    row = {
        'radius': radius, 'sites': lattice.volume, 'N2_dimension': lattice.dimension,
        'region_sites': [[0, 0, 0], [1, 0, 0]], 'regional_dimension': len(regional_basis(kept)),
        'model_time_intervals': [str(TIME), str(TIME)], 'split_steps_per_interval': STEPS,
        'complete_record_probabilities': probabilities.tolist(),
        'regional_sector_weights_by_record_N0_N1_N2': [sector_weights(state, kept) for state in reduced],
        'record_probability_restriction_residual': float(record_error),
        'nested_restriction_operator_residual': nested_error,
        'regional_source_to_sparse_oracle_cq_distance': density_cq_distance(reduced, exact_reduced),
        'full_source_to_sparse_oracle_cq_distance': pure_cq_distance(outputs, exact_outputs),
        'ideal_temporal_half_trace_bound': str(temporal),
        'ideal_fixed_spacing_infinite_history_half_trace_bound': str(min(F(1), temporal+spatial)),
        'retained_history_normalization_residual': abs(float(sum(probabilities))-1),
        'number_of_complete_records': len(outputs), 'local_sector_postselection': False,
        'floating_residuals_are_outward_certified': False,
        'apparatus_resources': apparatus_resources(lattice),
        'evolution_resources': source.resources(lattice, float(TIME), STEPS),
        'evolution_resources_include_both_intervals': True,
    }
    return row, lattice, outputs, reduced, exact_outputs


def reference_history():
    """An admitted joint pure input, with no reference reset between two reads."""
    lattice = WeylLattice(1)
    p, q = lattice.origin, lattice.neighbor
    kept = (p, p+1, q, q+1)
    state = np.column_stack((lattice.basis_state(p, p+1), 1j*lattice.basis_state(p, q+1)))/np.sqrt(2)
    input_reference = state.conj().T@state
    state = source.evolve(lattice, state, float(TIME), STEPS)
    first = source.read_ramsey(lattice, state, p+1, q+1)
    local = restrict_pure(lattice, state, kept)
    effect = ramsey_effect(kept, p+1, q+1)
    read_square_error = 0.
    for bit, branch in enumerate(first):
        k = effect if bit else np.eye(len(effect))-effect
        kr = np.kron(k, np.eye(2))
        read_square_error = max(read_square_error, float(np.linalg.norm(
            restrict_pure(lattice, branch, kept)-kr@local@kr.conj().T, 2)))
    outputs = []
    for bit, branch in enumerate(first):
        if bit:
            branch = source.phase_mode(lattice, branch, q+1, np.pi/2)
        branch = source.evolve(lattice, branch, float(TIME), STEPS)
        outputs.extend(source.read_ramsey(lattice, branch, p+1, q+1))
    reduced = [restrict_pure(lattice, branch, kept) for branch in outputs]
    final_reference = sum(restrict_regional(rho, kept, (), 2) for rho in reduced)
    probabilities = [float(np.trace(rho).real) for rho in reduced]
    return {
        'radius': 1, 'reference_dimension': 2, 'joint_input_preparation': 'supplied entangled N2/reference state',
        'reference_is_not_reset': True, 'all_four_record_blocks_retained': True,
        'subsequent_word_resources_bounded_by_radius1_level_counts': True,
        'source_preparation_for_entangled_input_executed': False,
        'complete_record_probabilities': probabilities,
        'local_instrument_restriction_square_residual': read_square_error,
        'reference_marginal_residual': float(np.linalg.norm(final_reference-input_reference.T, 2)),
        'retained_history_normalization_residual': abs(sum(probabilities)-1),
        'regional_sector_weights_by_record_N0_N1_N2': [sector_weights(rho, kept, 2) for rho in reduced],
    }


def algebra_controls():
    space = TwoFermions(4)
    h = np.zeros((4, 4)); h[2, 3] = h[3, 2] = 1
    x = space.lift_one_body(h).toarray()
    ex = conditional_expectation(space, x, (0, 1))
    defect = conditional_expectation(space, x@x, (0, 1))-ex@ex
    identity_error = float(np.linalg.norm(conditional_expectation(space, np.eye(6), (0, 1))-np.eye(6), 2))
    return {
        'fixed_N2_ambient_dimension': 6, 'retained_modes': [0, 1],
        'expectation_of_exterior_hopping_norm': float(np.linalg.norm(ex, 2)),
        'multiplicativity_defect_operator_norm': float(np.linalg.norm(defect, 2)),
        'multiplicativity_defect_trace': float(np.trace(defect).real),
        'unitality_residual': identity_error,
        'conditional_expectation_is_star_homomorphism': False,
        'full_Fock_CAR_algebra_is_identified_with_N2_compression': False,
        'homomorphic_retraction_obstruction_is_existing_theorem': 'OPH.QFT.no_scalar_hom_of_matrix',
        'actual_surviving_interface': 'OPH.QFT.CPRegionalNet',
        'separated_marginals_determine_joint_state': False,
        'opposite_phase_states_identical_separated_marginals': all(np.array_equal(
            restrict_pure(space, (space.basis_state(0, 2)+space.basis_state(1, 3))/np.sqrt(2), region),
            restrict_pure(space, (space.basis_state(0, 2)-space.basis_state(1, 3))/np.sqrt(2), region))
            for region in ((0, 1), (2, 3))),
        'opposite_phase_states_joint_half_trace_distance': 1,
    }


def cube_comparison(first, second):
    row1, small, out1, reduced1, exact1 = first
    row2, large, out2, reduced2, exact2 = second
    j = two_particle_embedding(small, large, cube_mode_map(small, large))
    compression = j.conj().T@large.hamiltonian@j-small.hamiltonian
    compression.eliminate_zeros()
    crossing = large.hamiltonian@j-j@small.hamiltonian
    # Frobenius norm is a diagnostic of non-intertwining, not a propagation rate.
    crossing_norm = float(sparse.linalg.norm(crossing))
    full_distance = pure_cq_distance([j@state for state in out1], out2)
    exact_full = pure_cq_distance([j@state for state in exact1], exact2)
    spatial = F(bounds.boundary(small.radius, 2*TIME)['complete_history_half_trace_upper'])
    temporal = F(row1['ideal_temporal_half_trace_bound'])+F(row2['ideal_temporal_half_trace_bound'])
    return {
        'radii': [small.radius, large.radius],
        'isometry_defect_nonzero_entries': (j.conj().T@j-sparse.eye(small.dimension)).nnz,
        'Hamiltonian_compression_defect_nonzero_entries': compression.nnz,
        'Hamiltonian_intertwining_defect_Frobenius_norm': crossing_norm,
        'compressed_Hamiltonian_identity_is_not_dynamic_intertwining': True,
        'source_complete_history_cq_distance_after_embedding': full_distance,
        'source_regional_history_cq_distance': density_cq_distance(reduced1, reduced2),
        'sparse_oracle_complete_history_cq_distance': exact_full,
        'ideal_exact_history_boundary_bound': str(spatial),
        'ideal_compiled_complete_history_bound': str(min(F(1), spatial+temporal)),
        'comparison_is_spatial_continuum': False,
    }


def source_bindings():
    from ..interacting_field.experiment import source_bindings as field_bindings
    paths = set(field_bindings()) | {
        'model/field_regions/regions.py', 'model/field_regions/experiment.py',
        'model/field_regions/__main__.py',
        'model/field_regions/README.md', 'model/field_regions/campaign.json',
        'Lean/QFT/FiniteCausalObserverNet.lean', 'Lean/QFT/CPRestrictionNet.lean',
        'Lean/QFT/NoncommutativeWitness.lean', 'Lean/QFT/TwoSlotCPNetWitness.lean',
    }
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(paths)}


def run():
    first, second = history_level(1), history_level(2)
    levels = [first[0], second[0]]
    reference = reference_history()
    comparison = cube_comparison(first, second)
    controls = algebra_controls()
    checks = {
        'normalized_complete_histories': all(row['retained_history_normalization_residual'] < 1e-10 for row in levels+[reference]),
        'regional_restriction_preserves_records': all(row['record_probability_restriction_residual'] < 1e-12 for row in levels),
        'nested_regional_maps_compose': all(row['nested_restriction_operator_residual'] < 1e-12 for row in levels),
        'reference_and_local_instrument_square': max(reference['reference_marginal_residual'], reference['local_instrument_restriction_square_residual']) < 1e-10,
        'compression_identity': comparison['Hamiltonian_compression_defect_nonzero_entries'] == 0,
        'isometry_identity': comparison['isometry_defect_nonzero_entries'] == 0,
        'exact_dynamic_intertwining_not_claimed': comparison['Hamiltonian_intertwining_defect_Frobenius_norm'] > 0,
        'regional_distance_contracts': comparison['source_regional_history_cq_distance'] <= comparison['source_complete_history_cq_distance_after_embedding']+1e-12,
        'fixed_spacing_history_within_ideal_bound_numerical_check': comparison['source_complete_history_cq_distance_after_embedding'] < float(F(comparison['ideal_compiled_complete_history_bound'])),
        'homomorphic_restriction_control_detected': controls['multiplicativity_defect_operator_norm'] == 1,
        'all_local_particle_sectors_retained': all(all(len(weights) == 3 for weights in row['regional_sector_weights_by_record_N0_N1_N2']) for row in levels),
    }
    return {'schema': 'oph.field-regions.v1', 'owner_gates': ['CP-08', 'CP-13'],
            'model_action': 'dGamma[-(1/8) sigma.sin(k)] + sum_x n_x0 n_x1',
            'time_parameter': 'supplied Hamiltonian model coordinate',
            'levels': levels, 'reference_history': reference, 'nested_cube_comparison': comparison,
            'algebra_controls': controls, 'checks': checks, 'all_pass': all(checks.values()),
            'physical_Lorentzian_or_continuum_attachment_claimed': False,
            'premise_register_dispositions_changed': False,
            'source_bindings': source_bindings()}
