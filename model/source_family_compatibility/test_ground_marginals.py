"""Source algebra, complete records and marginal-bound negative controls."""
from fractions import Fraction as F
import json

import numpy as np
import pytest
import sympy as s

from . import ground_marginals as g


def test_absolute_ground_projector_includes_degenerate_and_other_charge_states():
    matrix = s.diag(2, -1, -1, 3)
    value, rho, rank = g.algebraic_ground_projector(matrix)
    assert value == -1 and rank == 2
    assert rho == s.diag(0, s.Rational(1, 2), s.Rational(1, 2), 0)
    assert s.trace(rho) == 1


def test_exact_dimer_projector_matches_independent_spectral_source():
    from .local_spectral import dimer_source
    h, expected, _, _ = dimer_source()
    value, rho, rank = g.algebraic_ground_projector(h)
    assert rank == 1 and rho == expected
    assert s.simplify(value-(4-s.sqrt(17))/8) == 0
    assert g.finite_source_controls()['all_Fock_sectors_used'] is True


def test_whole_source_bond_compression_and_trace_factor():
    row = g.finite_source_controls()
    assert row['whole_bond_norm'] == F(1, 8)
    assert row['one_particle_trace_squared'] == F(1, 64)
    assert row['single_occupation_hopping_compression_exactly_zero']
    # Counting only one spin gives the wrong kinetic trace-square factor.
    assert row['one_particle_trace_squared'] != F(1, 128)


def test_finite_source_controls_return_detached_rationals():
    before = repr(g.finite_source_controls())
    changed = g.finite_source_controls()
    changed['whole_bond_norm']._numerator = 999
    assert repr(g.finite_source_controls()) == before


def test_surface_bound_counts_all_missing_bonds_and_translates():
    row = g.energy_density_interval(10, F(1, 8), (F(-20), F(-19)),
                                    field=F(1, 4), observable_norm=F(2),
                                    spans=(2, 1, 0))
    assert row['missing_bond_norm_per_site'] == F(3, 80)
    assert row['omitted_perturbation_translate_fraction'] == F(7, 25)
    assert row['interval'] == (F(-20, 1000)-F(3, 80)-F(7, 50),
                               F(-19, 1000)+F(7, 50))
    assert row['finite_ground_substituted_for_infinite_marginal'] is False


def test_onsite_perturbation_has_no_removed_boundary_anchors():
    plus = g.energy_density_interval(10, F(1, 8), (F(-20), F(-19)),
                                     field=F(10**6))
    base = g.energy_density_interval(10, F(1, 8), (F(-20), F(-19)))
    assert plus == base


def test_pressure_secants_detect_nonunique_atomic_source_marginal():
    # Actual c=0 centered source: each site's two singly occupied ground
    # vectors give arbitrary spin-Z expectation. The ±h spin field has
    # full-source ground energy -|h|V for h<U/2.
    h, volume = F(1, 8), 8
    row = g.expectation_interval(2, F(), h, (F(), F()),
                                 (-h*volume, -h*volume),
                                 (-h*volume, -h*volume))
    assert row['expectation_interval'] == (-1, 1)
    assert row['unique_marginal_or_convergence_modulus_inferred'] is False


def test_source_surface_error_cannot_be_removed_by_tiny_field():
    # A single exact energy triple may yield a narrow algebraic secant;
    # its volume surface error must still be divided by the field step.
    energies = (F(), F())
    row = g.expectation_interval(2, F(1, 8), F(1, 1000),
                                 energies, energies, energies)
    assert row['expectation_interval'] == (-1, 1)


def test_pauli_trace_bound_against_actual_noncommuting_states():
    x = np.array([[0, 1], [1, 0]], complex)
    z = np.diag([1, -1])
    rho = (np.eye(2)+F(3, 5)*x.astype(object))/2
    sigma = (np.eye(2)+F(4, 5)*z.astype(object))/2
    trace_distance = np.linalg.eigvalsh(np.asarray(rho-sigma, complex))
    actual = np.abs(trace_distance).sum()/2
    reference = (F(1), F(), F(), F(4, 5))
    boxes = ((F(1), F(1)), (F(3, 5), F(3, 5)), (F(), F()), (F(), F()))
    row = g.pauli_marginal_bound(1, reference, boxes)
    assert row['half_trace_distance_squared_upper'] == F(1, 4)
    assert abs(actual**2-F(1, 4)) < 1e-14


def test_complete_record_error_keeps_original_small_probability():
    assert g.complete_record_interval((F(1, 100), F(1, 50)),
                                      F(1, 200), F(1, 500), F(1, 1000)) == (F(1, 500), F(7, 250))
    assert g.complete_record_interval((F(1, 10000), F(1, 10000)),
                                      F(1, 1000), F(), F()) == (0, F(11, 10000))


def test_default_uniform_faithfulness_is_source_rational_bound():
    row = g.default_faithfulness()
    assert row['doublon_density_strict_interval'] == (F(1, 90000), F(233, 1600))
    assert row['least_one_site_eigenvalue_strict_lower'] == F(1, 90000)
    assert F(1, 2)-row['doublon_density_strict_interval'][1] > F(1, 90000)
    assert row['thermodynamic_gap_assumed'] is False
    assert row['larger_region_faithfulness_inferred'] is False


@pytest.fixture(scope='module')
def packet():
    return g.certificate()


def test_existing_native_ground_and_six_pulses_are_distinguished(packet):
    row = packet['executed_native_local_records']
    assert row['existing_source_bindings_authenticated'] == 566
    assert row['existing_preparation_tape_reexecuted'] is False
    assert row['additional_native_tagging_pulses_executed'] == 6
    assert row['fresh_blank_tag_registers_consumed'] == 2
    assert row['pointer_measurement_assumed_to_reset_the_register'] is False
    assert row['full_field_rows_retained'] == 65536
    assert row['native_record_probabilities_renormalized'] is False


def test_all_raw_records_and_classical_failure_prefixes_retained(packet):
    row = packet['executed_native_local_records']
    assert len(row['complete_pointer_records']) == 36
    failures = row['complete_additional_first_failure_records']
    assert len(failures) == 12+72
    assert sum(x['stage'] == 'first_tag' for x in failures) == 12
    assert {x['first_pointer_record'] for x in failures if x['stage'] == 'second_tag'} == set(range(6))
    assert row['earlier_preparation_failure_archive_retained']
    proper = [x for x in row['complete_pointer_records']
              if all(j < 2 for j in x['raw_pointer_records'])]
    assert len(proper) == 4
    assert all(F(x['original_probability_interval'][0]) > F(3, 200) for x in proper)
    assert sum(F(x['original_probability_interval'][0]) for x in proper) < 1
    assert sum(F(x['original_probability_interval'][1]) for x in proper) > 1


def test_actual_cube_energy_enclosure_retains_large_surface_error(packet):
    row = packet['actual_cube_energy_density_enclosure']
    assert row['full_Fock_rows'] == 65536 and row['scalar_shift_K_to_Q'] == '4'
    assert row['missing_bond_norm_per_site'] == '3/16'
    lo, hi = map(F, row['interval'])
    assert hi-lo > F(3, 16)
    assert packet['thermodynamic_local_marginal_convergence_modulus_proved'] is False


def test_prepared_one_site_reference_has_coarse_actual_uniform_bound(packet):
    row = packet['prepared_one_site_thermodynamic_comparison']
    assert F(row['ideal_cube_to_symmetric_ground_halftrace_upper']) == F(12967, 50000)
    assert row['complete_prepared_one_site_halftrace_strict_upper'] == '13/50'
    assert row['every_selected_symmetric_centered_ground_covered']
    assert row['full_source_or_multisite_state_bound'] is False


def test_packet_results_are_detached(packet):
    frozen = json.dumps(packet, sort_keys=True)
    value = g.certificate()
    value['executed_native_local_records']['complete_pointer_records'][0]['original_probability_interval'][0] = '1'
    assert json.dumps(g.certificate(), sort_keys=True) == frozen


@pytest.mark.parametrize('call', [
    lambda: g.energy_density_interval(3, F(1), (F(), F())),
    lambda: g.energy_density_interval(2, F(-1), (F(), F())),
    lambda: g.energy_density_interval(2, F(1), (F(), F()), spans=(2, 0, 0)),
    lambda: g.expectation_interval(2, F(1), F(), (F(), F()), (F(), F()), (F(), F())),
    lambda: g.pauli_marginal_bound(1, (F(1),), ((F(1), F(1)),)),
    lambda: g.algebraic_ground_projector(s.Matrix([[0, 1], [0, 0]])),
])
def test_invalid_admissions_rejected(call):
    with pytest.raises(ValueError):
        call()
