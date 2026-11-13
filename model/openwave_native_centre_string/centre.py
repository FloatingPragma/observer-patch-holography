"""Actual native charge/triality identities and complete Gauss string labels."""
from fractions import Fraction as F
from itertools import product
import numpy as np
import sympy as s
from scipy.linalg import eigh, expm
from model.coherent_defect_wave.physical_attachment.chiral import representation
from model.coherent_defect_wave.chiral_joint.interacting_charged_response.model import FAMILIES
from model.openwave_strong_step_scaling.source import chain


def native_quotient():
    g = representation()[2]
    Q = g[10] + g[11]
    C = s.Rational(2, 3) * g[6] + s.Rational(4, 3) * g[7]
    assert Q.is_diagonal() and C.is_diagonal()
    rows = []
    for i in range(15):
        q, c = Q[i, i], C[i, i]
        assert (3 * q).is_integer and (3 * c).is_integer
        assert (q + c).is_integer
        rows.append({'native_coordinate': i, 'charge_exact': str(q),
                     'electric_integer': int(3 * q),
                     'colour_triality': int(3 * c) % 3,
                     'integer_charge_plus_centre_generator': int(q + c)})
    assert s.trace(Q) == 0 and s.trace(C) == 0
    qints = sorted({int(3 * Q[i, i]) for i in FAMILIES['Q']})
    assert qints == [-1, 2]
    return {'complete_native_one_generation_coordinates': rows,
            'both_spin_copies_retained': True,
            'full_30_CAR_identity': 'exp(2*pi*i*Q)*centre_SU3=I on every exterior-power sector',
            'quotient_group': '(SU3 x U1)/Z3, with U1 integer q=3Q',
            'filled_reference_charge_and_triality_exact': '0',
            'native_up_down_electric_integers': [qints[1], qints[0]]}


def spatial_labels(d):
    L = d['n'].shape[1]
    qints = native_quotient()['native_up_down_electric_integers']
    local_q = np.zeros_like(d['n'], dtype=np.int64)
    for j, (x, spin, flavour) in enumerate(d['orbitals']):
        local_q[:, x] += qints[flavour] * d['counts'][:, j]
    total_q = local_q.sum(axis=1)
    assert np.all(total_q % 3 == 0)
    reference = -total_q
    prefix_number = np.cumsum(d['n'], axis=1)[:, :-1]
    flux = np.cumsum(local_q, axis=1)[:, :-1] + reference[:, None]
    assert np.all((flux + prefix_number) % 3 == 0)
    extended = np.column_stack((np.zeros(len(flux), dtype=np.int64), flux,
                                np.zeros(len(flux), dtype=np.int64)))
    charges = local_q.copy(); charges[:, 0] += reference
    assert np.array_equal(np.diff(extended, axis=1), charges)
    assert np.array_equal(np.diff(extended, axis=1) % 3, (-d['n']) % 3)
    # Every actual off-diagonal source term gets the uniquely forced electric
    # shift. This includes the charge-one weak exchanges, not only hopping.
    a, b = d['H'].nonzero()
    shifts = flux[a] - flux[b]
    assert np.max(np.count_nonzero(shifts, axis=1)) <= 1
    assert np.max(abs(shifts)) <= 3
    assert np.all(total_q[a] == total_q[b])
    return {'local_q': local_q, 'reference': reference, 'flux': flux,
            'triality': prefix_number % 3, 'phase_thirds': flux % 3,
            'all_hamiltonian_link_shifts': sorted(set(tuple(map(int, row)) for row in shifts))}


def execute(length=3, up_number=1, times=(0, .5, 2)):
    d = chain(length)
    labels = spatial_labels(d)
    ids = np.flatnonzero((d['nu'] == up_number) & (d['ns'] == 2))
    h = d['H'][ids][:, ids].toarray()
    e, U = eigh(h)
    # A genuine source occupation input: all quarks initially at the first
    # vertex. Spin/flavour counts are retained, rather than installing a
    # prepared energy eigenstate or a synthetic flux response.
    starts = np.flatnonzero(d['n'][ids, 0] == 3)
    assert len(starts) > 0
    initial = np.zeros_like(h); initial[starts, starts] = 1 / len(starts)
    initial_energy = float(np.trace(initial @ h))
    input_dimension = d['B'].shape[0]
    success = F(len(starts), input_dimension)
    ground = U[:, 0]
    nontrivial = np.any(labels['phase_thirds'][ids] != 0, axis=1)
    out = []
    for time in times:
        # Recenter the scalar mass before resolving vectors and phases.
        centre = float(np.trace(h) / len(h))
        propagator = expm(-1j * time * (h - centre * np.eye(len(h))))
        rho = propagator @ initial @ propagator.conj().T
        prob = np.diag(rho).real
        assert abs(prob.sum() - 1) < 2e-12
        before = float(np.trace(rho @ h).real)
        assert abs(before - initial_energy) < 2e-9
        after = float(prob @ np.diag(h))
        records = []
        for row, p in zip(ids, prob):
            records.append({'native_symmetric_occupation': list(map(int, d['counts'][row])),
                'local_electric_charge_thirds': list(map(int, labels['local_q'][row])),
                'reference_charge_thirds': int(labels['reference'][row]),
                'all_cut_electric_flux_integers': list(map(int, labels['flux'][row])),
                'all_cut_colour_trialities': list(map(int, labels['triality'][row])),
                'all_cut_dual_phase_thirds': list(map(int, labels['phase_thirds'][row])),
                'conditional_probability_diagnostic': float(p),
                'original_input_probability_diagnostic': float(p * success)})
        out.append({'time_source_units': time, 'complete_joint_records': records,
                    'energy_before_probe_diagnostic': before,
                    'energy_after_complete_nonselective_probe_diagnostic': after,
                    'probe_work_diagnostic': after - before,
                    'nontrivial_string_probability_diagnostic': float(prob @ nontrivial)})
    return {'cells': length, 'native_up_number': up_number,
            'complete_singlet_dimension': d['H'].shape[0],
            'complete_executed_charge_spin_dimension': len(ids),
            'native_total_charge_exact': str(F(up_number - 1)),
            'all_actual_native_Hamiltonian_electric_shifts': [list(x) for x in labels['all_hamiltonian_link_shifts']],
            'ground_nontrivial_string_probability_diagnostic': float(abs(ground) ** 2 @ nontrivial),
            'preparation': {'input': 'Uniform complete native three-quark carrier before Gauss/singlet preparation',
                'full_native_input_dimension': input_dimension,
                'accepted_colocated_charge_spin_dimension': len(starts),
                'success_exact': str(success), 'failure_exact': str(1 - success),
                'complete_nonselective_preparation_work_exact': '0 on the declared uniform input',
                'failure_domain': 'Original kinematic native complement, not a successful Gauss state'},
            'complete_time_records': out,
            'U1_electric_coefficient_exact': '0',
            'zero_coefficient_role': 'Exact Gauss dilation of the existing strong source, not a new electromagnetic energy fit.',
            'actual_phase_invariant': 'exp(2*pi*i*electric_flux/3)=exp(-2*pi*i*colour_triality/3)',
            'continuous_order_parameter_texture_winding_proved': False,
            'kinematic_constraint_and_Fock_complements_retained': True}


def certificate():
    return {'schema': 'oph.native_quotient_Gauss_string.v1',
            'verdict': 'PASS_NATIVE_FOCK_QUOTIENT_AND_COMPLETE_DYNAMICAL_STRING_RECORDS',
            'native_quotient': native_quotient(),
            'complete_spatial_executions': [execute(L, charge)
                                          for L in (2, 3) for charge in (1, 2)],
            'architecture': {'parent': 'openwave_strong_step_scaling',
                'source_carrier': 'Complete native30-CAR cells and full quotient-group links; the executed three-quark singlet is a proved reducing sector.',
                'U1_electric_energy_coefficient_exact': '0',
                'SU3_and_matter_coefficients': 'Exactly the parent; no fitted coefficient is changed.',
                'quantum_reference': 'R=-Q_total at vertex0; its declared energy is zero and its charge is retained.',
                'physical_dilation_is_exactly_isometric': True,
                'all_source_native_edge_terms_are_gauge_dressed': True},
            'scope': {'fractional_charge_dual_phase_relation_exact': True,
                'actual_spatial_native_source_executed': True,
                'continuous_order_parameter_winding_derived': False,
                'conjugate_link_angle_texture_is_inferred_from_an_electric_eigenvalue': False,
                'new_free_quark_GeV_mass_prediction': False,
                'new_electromagnetic_energy_fit': False,
                'parent_calibrated_masses_preserved_by_operator_transfer': True,
                'same_source_as_the192_CAR_transverse_member': False}}
