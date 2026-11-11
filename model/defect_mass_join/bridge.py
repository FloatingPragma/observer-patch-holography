"""Source embeddings and the complete quantum-scalar contact tensor."""
from __future__ import annotations

from functools import lru_cache
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import linalg as la
import sympy as s

from .source import (CHANNELS, LEPTONS, Source, compact_source, local_interactions,
                     one_cell, one_cell_polynomial, lepton_mask, map_lepton_column_back,
                     car, higgs_memory, representation)
from ..yukawa_defect_vertex import vertex as lepton_readout
from ..defect_mode_overlap import overlap
from ..defect_scalar_modes import modes as scalar_modes
from ..particle_common_contract.provenance import import_closure
from . import encounter

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent


@lru_cache(maxsize=1)
def local_operator_checks():
    """Full thirty-mode polynomial identities, never a quark-free shortcut."""
    h = one_cell_polynomial()
    parts, couplings = car.quantum_yukawa()
    independent = car.add(*parts.values(), car.monomial(matrix=higgs_memory()[2]))
    independent = car.clean({key: value.subs({x: 1 for x in couplings.values()})
                             for key, value in independent.items()})
    if h != independent:
        raise ArithmeticError('the common one-cell member is not the literal full native law')
    if h != car.adjoint(h):
        raise ArithmeticError('the full one-cell member is not Hermitian')
    for index in range(12):
        if car.commutator(car.gauge_current(index), h):
            raise ArithmeticError('a native current fails on the full F30 law')
    expected = {}
    epsilon = s.Matrix([[0, 1], [-1, 0]])
    fields = higgs_memory()[0]
    for name, vertex in representation()[3].items():
        for k, tensor in enumerate(vertex['components']):
            coefficient = -s.kronecker_product(epsilon, tensor.conjugate())
            # The equation of motion uses the adjoint of the scalar multiplying
            # the annihilation pair: phi for down/e, phi^dagger for up.
            field = fields[k] if vertex['conjugate_Higgs'] else fields[k].H
            for (i, j), value in s.ImmutableSparseMatrix(coefficient).todok().items():
                term = car.scale(car.monomial(matrix=field), value)
                expected[i, j] = car.add(expected.get((i, j), {}), term)
    operators = [car.monomial(((i, False),)) for i in range(30)]
    nonzero = []
    for i, ci in enumerate(operators):
        derivative = car.commutator(ci, h)
        for j, cj in enumerate(operators):
            actual = car.anticommutator(derivative, cj)
            if actual != expected.get((i, j), {}):
                raise ArithmeticError(('full native contact identity failed', i, j))
            if actual:
                nonzero.append((i, j))
    return {'full_F30_normal_ordered_H_identity': True,
            'all_native_global_current_commutators': 12,
            'all_fermion_pair_contacts_checked': 900,
            'nonzero_oriented_native_contacts': len(nonzero),
            'equation': '{[c_(x,s,a),H],c_(y,t,b)}=-delta_xy D_x epsilon_st sum_channel,k conjugate(C_channel,k[a,b]) adjoint(Phi_channel,x,k)',
            'full_thirty_mode_physical_ground_control': 'only the texture-off one-cell member, on the native-current singlet domain, matches local_reference; no spatial ground transfer',
            'scalar_number_or_stiffness_term_used_as_Higgs_insertion': False}


@lru_cache(maxsize=1)
def lepton_control_embedding():
    """Every column of the exact y_u=y_d=0 one-cell reducing control."""
    reduced = one_cell(('charged_lepton',))
    full = one_cell()
    target = lepton_readout.source_data()['H']
    for mask in range(64):
        for scalar in range(5):
            actual = map_lepton_column_back(reduced.column(lepton_mask(mask), scalar))
            expected = {row: complex(value) for (row, col), value in target.todok().items()
                        if col == 5*mask+scalar}
            if actual != expected:
                raise ArithmeticError('a complete 320-state control column differs')
    all_unit_column = full.column(0, 0)
    quark_mask = ((1 << 30)-1) ^ sum(1 << i for i in LEPTONS)
    leakage = sum(abs(value)**2 for (mask, _), value in all_unit_column.items() if mask & quark_mask)
    if leakage <= 0:
        raise ArithmeticError('the full unit-channel non-reduction control vanished')
    return {'exact_control_parameters': {'texture': False, 'sites': [[0, 0, 0]],
                'up': 0, 'down': 0, 'charged_lepton': 1, 'scalar_frequency': 1},
            'complete_intertwined_source_columns': 320,
            'six_to_thirty_CAR_mode_map': list(LEPTONS),
            'scalar_embedding': 'the same five scalar coordinates',
            'all_unit_channel_quark_creation_norm_squared_on_blank': leakage,
            'empty_quarks_reduce_all_three_unit_channel_member': False,
            'control_ground_or_Z_H_transferred_to_different_members': False}


@lru_cache(maxsize=1)
def spatial_contact():
    """Every site and every quantum scalar mode, with all exterior leakage."""
    joined = compact_source()
    band = overlap.build()
    scalar = scalar_modes.build_scalar_source()
    if (joined.points != tuple(tuple(p) for p in scalar['points'])
            or joined.points != tuple(band['raw']['points'])):
        raise ArithmeticError('the scalar and fermion spatial sources differ')
    if not np.array_equal(joined.K, scalar['K']):
        raise ArithmeticError('the scalar stiffness is not the same operator')
    if la.norm(joined.h-band['h']) > 1e-11 or la.norm(joined.Omega-scalar['Omega']) > 1e-9:
        raise ArithmeticError('a source coefficient was silently substituted')
    w, c = band['W'], band['c']
    # The local spin epsilon with the fixed stagger phase is exactly C.
    literal = np.zeros_like(c)
    epsilon = np.array([[0., 1.], [-1., 0.]])
    for site, phase in enumerate(joined.stagger):
        literal[2*site:2*site+2, 2*site:2*site+2] = phase*epsilon
    if la.norm(literal-c) != 0:
        raise ArithmeticError('the microscopic pair phase does not reproduce the defect source')
    profiles = np.asarray(scalar['modes'])
    raw = np.einsum('xr,xij->rij', profiles, band['kernel'])
    canonical = np.sqrt(2*np.asarray(scalar['frequencies']))[:, None, None]*raw
    rows = []
    for mode, (frequency, profile) in enumerate(zip(scalar['frequencies'], profiles.T, strict=True)):
        response = overlap.profile_response(profile, value=band)
        if la.norm(response['vertex']-raw[mode]) > 1e-11:
            raise ArithmeticError('a scalar mode contact differs from the full local operator')
        rows.append({'mode': mode, 'free_scalar_frequency': float(frequency),
            'raw_spatial_singular_values': la.svdvals(raw[mode]).tolist(),
            'free_vacuum_canonical_singular_values': la.svdvals(canonical[mode]).tolist(),
            'complete_exterior_leakage': response['leakage'],
            'free_vacuum_canonical_exterior_leakage': float(np.sqrt(2*frequency)*response['leakage'])})
    return {'source': joined.descriptor(), 'W': w, 'site_kernel': band['kernel'],
            'scalar_modes': profiles, 'raw_mode_kernel': raw, 'canonical_mode_kernel': canonical,
            'complete_mode_records': rows,
            'fermion_band_metric_error': float(la.norm(w.conj().T@w-np.eye(4))),
            'full_native_band_annihilators': 60, 'full_scalar_modes_per_component': joined.sites,
            'operator_identity': 'K_band=-sum_x Y_x tensor conjugate(C_channel,k) adjoint(Phi_channel,x,k)',
            'raw_scalar_mode_identity': 'Y_j=sum_x u_xj Y_x',
            'free_vacuum_canonical_identity': 'Y_j,canonical=sqrt(2 omega_j) Y_j',
            'interacting_ground_scalar_residue_substituted_from_vacuum': False,
            'complete_band_is_reducing_under_local_quantum_Higgs_interaction': False}


def reference_and_records():
    """One source/state/reference/clock, with complete finite preparation."""
    source = compact_source()
    return {'source': 'the complete compact Source member, all three primitive channels one',
        'source_singlet_projector': 'P=1_{0}(sum_A Q_A^2), the twelve GLOBAL native currents',
        'local_gauge_field_or_local_Gauss_constraint_supplied': False,
        'ground_definition': 'P0 is the lowest-energy spectral projector of P H P on ran P; rho=P0/rank(P0)',
        'ground_executed_for_complete_spatial_member': False,
        'reference': 'full conjugate source carrier, idle Hamiltonian zero, retained in every output',
        'initial_preparation': 'vec(P)/sqrt(rank(P)); declared finite singlet/parity-preserving preparation service',
        'complete_preparation_energy_probabilities': 'p(E)=rank(P P_E)/rank(P); retain EVERY energy outcome',
        'ground_branch_original_probability': 'rank(P0)/rank(P), never replaced by one',
        'ground_reference_vector': 'vec(P0)/sqrt(rank(P0))',
        'clock': 'U(t)=exp(-itH), hbar=1 in this one source clock',
        'Hamiltonian_entire_carrier_norm_bound': str(source.norm_bound()),
        'charged_fields': 'a_mu=sum_x,s conjugate(W_xs,alpha)c_x,s,internal, all60 band coordinates',
        'complete_charged_isometry': '(sum_mu a_mu^dagger tensor |mu dual,+> + a_mu tensor |mu,->)/sqrt(60)',
        'complete_charged_Gram': 'V_F^dagger V_F=I on the ENTIRE CAR/scalar carrier',
        'two_point_energy_matrices': 'C+_mu,nu(E)=Tr rho a_mu P_E a_nu^dagger; C-_mu,nu(E)=Tr rho a_nu^dagger P_E a_mu',
        'full_energy_resolution': 'every P_E of the same full H, including complement-mediated and zero-weight outcomes',
        'fermion_CAR_normalization': 'sum_E(C++C-)=I60',
        'charged_excitation_energies_presumed_positive': False,
        'scalar_isometry': 'phi_xk/sqrt(4N), phi_xk^dagger/sqrt(4N), K_noop=sqrt(1-1/(4N))*P_scalar_occupied',
        'scalar_pointer_completeness': 'sum K_r^dagger K_r=I; both component reference spaces and no-operation retained',
        'actual_scalar_pulse': 'X=sum_x,k(phi_xk^dagger tensor |x,k dual><0| + h.c.)',
        'three_point_readout': 'Tr rho {[a_mu,H],a_nu} phi_xk^dagger or its conjugate orientation, with all required charge-reference factors',
        'contact_reconstruction': 'T=Gamma C_scalar on the full field-coordinate span; any null Gram directions remain explicitly unresolved',
        'interacting_scalar_normalization': 'C0=<[phi,phi^dagger]>; Z1=<[phi,[H,phi^dagger]]>, both computed on this rho before choosing an amputation',
        'one_cell_lepton_normalization_reused': False,
        'free_scalar_vacuum_normalization_reused_for_interacting_ground': False}


def bosonic_source_content():
    """The scalar spectrum actually supplied, with its reference kept explicit."""
    source = compact_source()
    modes = scalar_modes.build_scalar_source()
    if la.norm(source.Omega-modes['Omega']) > 1e-9:
        raise ArithmeticError('the scalar spectrum does not use the joined source frequency')
    # Reconstruct the complete free retarded function and first moment from
    # every residue; the coupled many-body ground is not this reference.
    z = 1+2j
    u, omega = modes['modes'], modes['frequencies']
    retarded = (u*(2*omega/(z*z-omega*omega)))@u.T
    literal = 2*source.Omega@la.inv(z*z*np.eye(source.sites)-source.K)
    first_moment = (u*(2*omega))@u.T
    errors = {'complete_free_response': float(la.norm(retarded-literal, 2)),
              'complete_free_first_moment': float(la.norm(first_moment-2*source.Omega, 2))}
    if max(errors.values()) > 2e-8:
        raise ArithmeticError('the complete scalar response or moment differs')
    return {'source_bosonic_carrier': 'one scalar blank plus particle and antiparticle Higgs doublets at every site; total occupation<=1',
        'scalar_complex_field_coordinates': 2*source.sites,
        'scalar_carrier_dimension': source.scalar_dimension,
        'scalar_suboperator_blank_energy': 0,
        'scalar_suboperator_blank_rank': 1,
        'all_scalar_suboperator_frequencies': list(map(float, omega)),
        'scalar_suboperator_rank_per_spatial_mode': 4,
        'complete_scalar_suboperator_rank': 1+4*len(omega),
        'free_reference': 'scalar blank for H_H alone, tensor an idle conjugate reference',
        'free_reference_is_joint_interacting_ground': False,
        'retarded_free_response': 'D_phi(z)=2 Omega (z^2 I-K)^(-1), for each doublet component',
        'all_positive_and_negative_frequency_poles_retained': True,
        'equal_time_free_commutator_expectation': 'zero',
        'free_first_energy_moment': '2 Omega',
        'free_canonical_field': 'q=(2 Omega)^(-1/2) phi; D_q(z)=(z^2 I-K)^(-1)',
        'same_operator_numeric_checks': errors,
        'global_native_current_count': 12,
        'dynamical_gauge_link_or_field_operator_count': 0,
        'global_current_matrices_are_gauge_boson_oscillators': False,
        'comparison_quantum_link_source': 'chiral_joint finite_program has 60-state link factors, covariant path operators, link energy and tensor-product local five-state scalars; those operators are absent from this member',
        'identical_clock_reference_embedding_of_quantum_link_source_proved': False,
        'one_cell_common_native_polynomial_does_not_transfer_gauge_link_spectrum': True,
        'interacting_scalar_poles_or_residues_executed': False,
        'W_Z_photon_gluon_mass_records_supplied': False,
        'scalar_suboperator_frequencies_identified_as_physical_rest_masses': False}


def serial(value):
    if isinstance(value, s.Basic):
        return str(value)
    if isinstance(value, np.ndarray):
        raise TypeError('large tensor arrays belong in the numeric operator packet')
    if isinstance(value, dict):
        return {k: serial(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [serial(v) for v in value]
    return value


def bindings():
    paths = ['model/defect_mass_join/source.py', 'model/defect_mass_join/bridge.py',
             'model/defect_mass_join/encounter.py', 'model/defect_mass_join/test_join.py',
             'model/defect_mass_join/__init__.py', 'model/defect_mass_join/__main__.py',
             'model/defect_mass_join/README.md', 'model/requirements.txt',
             'model/defect_mass_join/PROOF.md', 'model/defect_scalar_modes/modes.py',
             'model/compact_defect_pair/spectral.py', 'model/compact_defect_pair/topology.py',
             'model/coherent_defect_wave/chiral_joint/higgs_composites/model.py',
             'model/coherent_defect_wave/chiral_higgs/model.py',
             'model/coherent_defect_wave/physical_attachment/chiral.py',
             'model/common_source/chiral_candidate/action.py',
             'model/observer_particles/results/representation_identification_operators.json',
             'model/observer_particles/results/unified_yukawa_operators.json']
    for name in paths:
        if not (ROOT/name).is_file():
            raise FileNotFoundError('missing declared joined-source input: '+name)
    paths = import_closure(paths)
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


def certificate():
    pins = bindings()
    contact = spatial_contact()
    result = serial({'schema': 'oph.defect_mass_join.v1', 'source_bindings': pins,
        'local_full_domain_checks': local_operator_checks(),
        'literal_lepton_control': lepton_control_embedding(),
        'joined_source': contact['source'],
        'complete_spatial_scalar_contacts': {key: value for key, value in contact.items()
            if key not in ('source', 'W', 'site_kernel', 'scalar_modes', 'raw_mode_kernel', 'canonical_mode_kernel')},
        'state_reference_clock_and_complete_records': reference_and_records(),
        'three_site_complete_source_propagation_control': encounter.execute(),
        'bosonic_source_content_and_complete_free_scalar_control': bosonic_source_content(),
        'scope': {'one_explicit_full_finite_Hamiltonian_family': True,
            'all_scalar_modes_and_CAR_occupations_retained': True,
            'full_Fock_anomalous_contact_vertex_derived': True,
            'finite_time_scalar_to_distant_fermion_pair_record_executed': True,
            'localized_family_matrix_supplied': False, 'individual_mass_fitting': False,
            'complete_spatial_interacting_ground_and_poles_executed': False,
            'amputated_frequency_dependent_1PI_vertex_executed': False,
            'physical_rest_mass_or_Lorentz_limit': False}})
    if bindings() != pins:
        raise ArithmeticError('the joined source changed during calculation')
    return result


def write():
    result = certificate()
    (HERE/'results.json').write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    contact = spatial_contact()
    np.savez_compressed(HERE/'operators.npz', **{name: contact[name] for name in
        ('W', 'site_kernel', 'scalar_modes', 'raw_mode_kernel', 'canonical_mode_kernel')})
    return result
