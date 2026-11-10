"""Replay the finite joint interaction and complete operational countercontrols."""
from fractions import Fraction as F
import hashlib
from pathlib import Path
import sympy as s

from model.particle_common_contract.provenance import import_closure
from .model import PROGRAM_ID, JointObserver, SparseState, combine
from . import preparation, gauss

ROOT = Path(__file__).resolve().parents[2]


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def input_hashes():
    paths = {p for p in Path(__file__).parent.glob('*.py') if not p.name.startswith('test_')}
    paths.update(ROOT/p for p in (
        'model/requirements.txt', 'model/particle_common_contract/provenance.py',
        'evidence/native_geometric_source_20260925/spec.json',
        'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
        'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py'))
    paths = {ROOT/name for name in import_closure({str(p.relative_to(ROOT)) for p in paths})}
    paths = {p for p in paths if p.name in {q.name for q in p.parent.iterdir()}}
    return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}


def field_distribution(observer, state):
    branches = observer.field_measurement(state)
    require(sum(branch.norm_squared() for _, branch in branches) == state.norm_squared(),
            'complete link measurement')
    return {str(record['outcome']): str(branch.norm_squared()) for record, branch in branches}


def interaction_certificate():
    observer = JointObserver()
    g0 = observer.initial_word
    g1 = g0 ^ (1 << observer.law.edges.index((0, 5)))
    # Freeze inputs before reading outputs; inverse-wave preparation is a declared service.
    pure_outputs = []
    pure_face_outputs = []
    for port in (5, 6):
        prepared = observer.inverse_wave(observer.embedded_port(g0, port))
        output = observer.step(prepared)
        require(output.norm_squared() == prepared.norm_squared() == 1, 'unit norm')
        require(observer.inverse_step(output) == prepared, 'full joint inverse')
        pure_outputs.append(field_distribution(observer, output))
        face_branches = observer.region_face_measurement(output)
        require(sum(branch.norm_squared() for _, branch in face_branches) == 1, 'complete invariant flux detector')
        pure_face_outputs.append({str(record['face_index']): str(branch.norm_squared())
                                  for record, branch in face_branches})
    require(pure_outputs[0]['4194307'] == '1' and pure_outputs[1]['1'] == '1',
            'wave preparation controls the full-tick defect output')

    target = observer.embedded_ports(g0, (5, 6))
    prepared = observer.inverse_wave(target)
    entangled = observer.step(prepared)
    dist = field_distribution(observer, entangled)
    require(dist['1'] == dist['4194307'] == '1/2', 'two equal Schmidt weights')
    # Removing only the controlled transport coherently restores the superposition.
    recombined = observer.controls(entangled, inverse=True)
    coherent = observer.wave_preparation_measurement(recombined, (5, 6))
    coherent_probabilities = [branch.norm_squared() for _, branch in coherent]
    require(coherent_probabilities == [1, 0], 'coherent inverse-control return')
    decohered_probabilities = [s.Integer(0), s.Integer(0)]
    complete_history = []
    for record, branch in observer.region_face_measurement(entangled):
        inverse = observer.controls(branch, inverse=True)
        for i, (detector, output) in enumerate(observer.wave_preparation_measurement(inverse, (5, 6))):
            probability = output.norm_squared()
            decohered_probabilities[i] += probability
            complete_history.append({'face_index': record['face_index'], 'return': detector['outcome'],
                                     'probability': str(probability)})
    require(decohered_probabilities == [F(1, 2), F(1, 2)], 'retaining a flux-face readout changes interference')
    require(sum(F(row['probability']) for row in complete_history) == 1, 'all10 history outcomes')

    address = observer.addresses.index((0, 5))
    wave_results = []
    for g in (g0, g1):
        output = observer.step(observer.embedded_ports(g, (0, 5)))
        branches = observer.wave_basis_measurement(output, 0, address=address)
        probabilities = [branch.norm_squared() for _, branch in branches]
        require(sum(probabilities) == 1, 'wave detector complement retained')
        wave_results.append(list(map(str, probabilities)))
    require(wave_results == [['1/120', '119/120'], ['0', '1']], 'full-tick reciprocal interaction')
    # The triplet sees identity transport on every central link word.
    triplet_probabilities = []
    for g in (g0, g1):
        output = observer.step(observer.embedded_ports(g, (0, 5), internal=0))
        triplet_probabilities.append(observer.wave_basis_measurement(output, 0, internal=0, address=address)[0][1].norm_squared())
    require(triplet_probabilities[0] == triplet_probabilities[1], 'central-triplet countercontrol')
    return {
        'wave_to_defect': {
            'initial_link_word': g0,
            'preparations': ['W_g* A|5,internal3>', 'W_g* A|6,internal3>'],
            'complete_link_distributions': pure_outputs,
            'deterministic_final_link_words': [4194307, 1],
            'complete_gauge_invariant_face_distributions': pure_face_outputs,
            'deterministic_final_face_indices': [17, 0],
            'actual_full_tick_executed': True,
        },
        'defect_to_wave': {
            'initial_link_words': [g0, g1],
            'same_wave_preparation': 'A(|0,internal3>+|5,internal3>)/sqrt(2)',
            'wave_detector': {'port': 0, 'internal': 3, 'parity': 0, 'oriented_address': [0, 5]},
            'yes_and_complement_probabilities': wave_results,
            'triplet_comparison_yes_probabilities': list(map(str, triplet_probabilities)),
            'actual_full_tick_executed': True,
        },
        'coherence': {
            'domain': 'ambient register factorization; the preparation-return projector is not a constrained gauge-invariant detector',
            'input': '|g=1> tensor W_g* A(|5,internal3>+|6,internal3>)/sqrt(2)',
            'schmidt_weights': ['1/2', '1/2'],
            'final_link_distribution': dist,
            'coherent_inverse_controls_return_and_complement': list(map(str, coherent_probabilities)),
            'after_complete_face_readout_return_and_complement': list(map(str, decohered_probabilities)),
            'complete_face_then_return_history': complete_history,
            'recorded_idle_move_flag': False,
        },
    }


def constraint_execution_certificate():
    observer = JointObserver()
    wave_probabilities, face_probabilities = [], []
    def gamma(state):
        return SparseState(gauss.gamma_amplitudes(state.as_dict()), state.norm_denominator)
    def plus(state):
        return SparseState(gauss.project_amplitudes(state.as_dict()), state.norm_denominator)
    for word in (1, 3):
        state = observer.embedded_ports(word, (0, 5))
        require(observer.step(gamma(state)) == gamma(observer.step(state)), 'executed gauge covariance')
        output = observer.step(plus(state))
        require(gamma(output) == output, 'preserved positive Gauss sector')
        branches = observer.wave_basis_measurement(output, 0, address=observer.addresses.index((0, 5)))
        require(all(gamma(branch) == branch for _, branch in branches), 'invariant complete wave detector')
        wave_probabilities.append(str(branches[0][1].norm_squared()/output.norm_squared()))
    for port in (5, 6):
        state = plus(observer.inverse_wave(observer.embedded_port(1, port)))
        output = observer.step(state)
        require(gamma(output) == output and observer.inverse_step(output) == state, 'physical inverse')
        branches = observer.region_face_measurement(output)
        require(all(gamma(branch) == branch for _, branch in branches), 'invariant complete face detector')
        require(sum(branch.norm_squared() for _, branch in branches) == output.norm_squared(), 'physical detector completeness')
        face = next(branch for record, branch in branches if record['face_index'] == 17)
        face_probabilities.append(str(face.norm_squared()/output.norm_squared()))
    require(wave_probabilities == ['1/120', '0'] and face_probabilities == ['1', '0'], 'constrained reciprocal controls')
    return {'positive_Gauss_subspace_used': True,
            'conditional_wave_yes_probabilities': wave_probabilities,
            'conditional_target_flux_face_probabilities': face_probabilities,
            'all_detector_branches_preserve_constraint': True,
            'complete_inverse_checked': True}


def invariant_coherence_certificate():
    observer = JointObserver()
    raw = observer.inverse_wave(observer.embedded_ports(1, (5, 6)))
    # P+ has norm squared1/2 on this raw preparation. Multiply by sqrt2
    # through the common normalization to retain exact rational coefficients.
    projected = gauss.project_amplitudes(raw.as_dict())
    prepared = SparseState(projected, raw.norm_denominator//2)
    def gamma(state):
        return SparseState(gauss.gamma_amplitudes(state.as_dict()), state.norm_denominator)
    require(prepared.norm_squared() == 1 and gamma(prepared) == prepared, 'normalized constrained preparation')
    output = observer.step(prepared)
    recombined = observer.controls(output, inverse=True)
    direct = gauss.covariant_preparation_return(observer, recombined, (5, 6))
    require(all(gamma(branch) == branch for _, branch in direct), 'invariant return branches')
    direct_probabilities = [branch.norm_squared() for _, branch in direct]
    require(direct_probabilities == [1, 0], 'constrained coherent return')
    history = []
    totals = [s.Integer(0), s.Integer(0)]
    for face_record, face_branch in observer.region_face_measurement(output):
        require(gamma(face_branch) == face_branch, 'invariant recorded face branch')
        inverse = observer.controls(face_branch, inverse=True)
        for index, (record, branch) in enumerate(gauss.covariant_preparation_return(observer, inverse, (5, 6))):
            require(gamma(branch) == branch, 'both post-record return branches invariant')
            probability = branch.norm_squared()
            totals[index] += probability
            history.append({'face_index': face_record['face_index'], 'return': record['outcome'],
                            'probability': str(probability)})
    require(totals == [F(1, 2), F(1, 2)] and sum(totals) == 1, 'constrained complete recording contrast')
    # Test the operator intertwining on a non-invariant input, separately
    # from branch preservation on a positive-sector input.
    left = gauss.covariant_preparation_return(observer, gamma(raw), (5, 6))
    right = gauss.covariant_preparation_return(observer, raw, (5, 6))
    require(all(a == gamma(b) for (_, a), (_, b) in zip(left, right)), 'return projector commutes with Gamma')
    require(sum(branch.norm_squared() for _, branch in right) == raw.norm_squared(), 'complete invariant apparatus')
    return {
        'domain': 'positive interior Gauss sector with declared conditional-frame return apparatus',
        'normalized_preparation': 'sqrt(2) P_plus (|1> tensor W_1* A(|5,3>+|6,3>)/sqrt(2))',
        'coherent_return_and_complement': list(map(str, direct_probabilities)),
        'after_complete_face_readout_return_and_complement': list(map(str, totals)),
        'complete_face_then_return_history': history,
        'return_projector_commutes_with_Gamma': True,
        'all_recorded_branches_preserve_positive_constraint': True,
        'idle_joint_U_changed': False,
    }


def structure_certificate():
    observer = JointObserver()
    law = observer.law
    initial = observer.initial_word
    outside = law.syndrome(initial) & ~law.region_mask
    permutations = {}
    fixed_counts = []
    for edge in law.interior:
        permutation = [law.transport(g, edge)[0] for g in observer.link_words]
        require(set(permutation) == set(observer.link_words), 'closed transport permutation')
        require(all(law.transport(h, edge)[0] == g for g, h in zip(observer.link_words, permutation)), 'involution')
        permutations[str(edge)] = permutation
    for g in observer.link_words:
        require(g & law.boundary_mask == initial & law.boundary_mask, 'fixed boundary links')
        require(law.region_charge(g) == 1, 'protected odd regional flux')
        require(law.syndrome(g) & ~law.region_mask == outside, 'fixed exterior flux')
        require((law.syndrome(g) & law.region_mask).bit_count() == 1, 'one regional defect')
    for port in range(12):
        fixed = 0
        for g in observer.link_words:
            out = observer.controls(observer.basis_state(g, port))
            fixed += next(iter(out.as_dict()))[0] == g
        fixed_counts.append(fixed)
    # Sum the literal diagonal coefficient of QR on each wave coordinate.
    # R acts only on seam address; Q preserves it, so no off-diagonal term
    # contributes to (QR)_ii. Central signs occur off the Q diagonal.
    per_port = []
    for port in range(12):
        per_port.append(sum((F(-29, 30)*(F(1, 2) if port in ends else 1)
                             for ends in observer.addresses for _ in range(10)), F(0)))
    trace = sum(n*t for n, t in zip(fixed_counts, per_port))
    require(trace == F(preparation.certificate()['joint_unitary_trace']), 'same literal-U trace')
    return {
        'dimension': observer.dimension,
        'wave_dimension_per_link_word': 7200,
        'link_words': list(observer.link_words),
        'interior_edge_indices': list(law.interior),
        'interior_edges': [list(law.edges[e]) for e in law.interior],
        'transport_permutations': permutations,
        'control_fixed_word_counts': fixed_counts,
        'unitary_trace': str(trace),
        'boundary_charge': 1,
        'interior_flux_count': 1,
        'internal_Q5_diagonal': [2, 2, 2, -3, -3],
        'all_30_link_signs_used_in_wave_repair': True,
        'complete_joint_clock': 'U=C4 C3 C2 C1 C0 direct_sum_g(Q_g R)',
        'forward_event_order': ['wave'] + [list(law.edges[e]) for e in law.interior],
    }


def certificate():
    return {
        'schema': 'oph.coherent_defect_wave.v1', 'program_id': PROGRAM_ID,
        'input_sha256': input_hashes(),
        'source': structure_certificate(),
        'reference': preparation.certificate(),
        'gauss': dict(gauss.certificate(), executed_constraint_readouts=constraint_execution_certificate()),
        'interaction': dict(interaction_certificate(), constrained_coherence=invariant_coherence_certificate()),
        'ancestry': {
            'reused': ['pinned icosahedral carrier', 'reversible seam reflection', 'fixed-boundary transport grammar'],
            'declared': ['complex amplitude framework', 'coherent link register', 'one-defect boundary sector',
                         'occupation-controlled transport', 'fixed joint tick order', 'complete filter and detectors'],
            'particle_data_calibrations': [],
        },
        'scope': {
            'finite_bidirectional_interaction': True,
            'protected_regional_flux': True,
            'same_joint_clock_reference': True,
            'one_interior_gauss_constraint': True,
            'continuous_dynamical_gauge_theory': False,
            'complete_boundary_gauss_constraints': False,
            'complete_joint_effect_tomography': False,
            'interacting_3plus1_continuum': False,
            'physical_particle_species': False,
            'physical_energy_or_GeV_mass': False,
        },
    }
