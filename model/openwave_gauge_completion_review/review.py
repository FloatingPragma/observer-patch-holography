"""Independent CAR, group-moment and infinite-rotor comparison calculations."""
import json
import ast
import hashlib
from fractions import Fraction as F
import numpy as np
import sympy as s
from scipy.linalg import eigh, expm
from model.openwave_pure_su2_source.custody import ROOT


def car_preparation():
    annihilators = []
    for bit in range(6):
        a = np.zeros((64, 64), dtype=np.int64)
        for state in range(64):
            if (state >> bit) & 1:
                a[state^(1 << bit), state] = (-1)**((state & ((1 << bit)-1)).bit_count())
        annihilators.append(a)
    B = sum((annihilators[3+i].T@annihilators[i] for i in range(3)), np.zeros((64, 64), dtype=np.int64))
    P = np.diag([int(x < 8) for x in range(64)])
    K = B@P; G = K.T@K
    failure = 3*P-G
    assert np.linalg.eigvalsh(failure).min() == 0
    assert np.array_equal(G+failure+3*(np.eye(64, dtype=np.int64)-P), 3*np.eye(64, dtype=np.int64))
    assert max(np.diag(G)) == 3
    u = np.zeros(64, dtype=np.int64); u[7] = 1
    w = B@u
    assert int(np.dot(B@w, B@w)) == 4*int(np.dot(w, w))
    return {'full_Fock_dimension': 64, 'exact_complete_instrument_Gram': True,
            'unprojected_pump_squared_gain_exact': '4/3',
            'reference_test_is_necessary': True}


def group_moments(receipt):
    # Schur orthogonality and 3xbar3=1+8 give <chi8>=2/3 in this class state.
    # Conjugation invariance fixes the two-index second moment from unitarity
    # and <|chi3|²>=1+<chi8>=5/3, independently of entry-wise Weingarten sums.
    A, B = s.symbols('A B')
    solution = s.solve([A+3*B-1, 9*A+3*B-s.Rational(5, 3)], (A, B))
    assert solution == {A: s.Rational(1, 12), B: s.Rational(11, 36)}
    z = -s.Rational(1, 2)+s.I*s.sqrt(3)/2
    for row in receipt['exact_SU3_Haar_moments']:
        r = row['centre_index']
        for i in range(3):
            for j in range(3):
                assert s.simplify(s.sympify(row['all_nine_U_means_exact'][3*i+j])-(z**r/3 if i == j else 0)) == 0
                assert s.sympify(row['all_nine_U_absolute_square_means_exact'][3*i+j]) == (solution[A] if i == j else 0)+solution[B]
        assert row['TrU_absolute_square_mean_exact'] == '5/3'
    for row in receipt['complete_spatial_records']:
        probabilities = [s.sympify(x['probability_exact']) for x in row['all_spatial_flux_colour_pointer_records']]
        assert s.simplify(sum(probabilities)) == 1
        assert all(x.is_nonnegative is True for x in probabilities)
    return {'independent_class_function_Schur_solution_exact': {str(k): str(v) for k, v in solution.items()},
            'all_colour_flux_pointer_probabilities_rechecked': True,
            'singlet_postselection_not_used': True}


def rotor_matrix(cap, q):
    states = [(x, n) for n in range(-cap, cap+1) for x in range(4)]
    where = {v: i for i, v in enumerate(states)}
    H = np.zeros((len(states), len(states)))
    for col, (x, n) in enumerate(states):
        flux = [n-q if edge < x else n for edge in range(4)]
        H[col, col] = sum(v*v for v in flux)/32+1/4
        for direction in (-1, 1):
            y = (x+direction)%4
            k = n+q if x == 0 and direction == -1 else n-q if x == 3 and direction == 1 else n
            if (y, k) in where:
                H[where[y, k], col] -= 1/16
            if (x, n+direction) in where:
                H[where[x, n+direction], col] -= 1/16
    assert np.array_equal(H, H.T)
    return states, where, H


def rotor_comparison(receipt):
    states, where, H = rotor_matrix(36, -3)
    lowest = eigh(H, eigvals_only=True)[0]
    n = np.arange(-80, 81)
    vacuum = np.diag((4*n*n+4)/32)+np.diag(np.full(160, -1/16), 1)+np.diag(np.full(160, -1/16), -1)
    vacuum_lowest = eigh(vacuum, eigvals_only=True)[0]
    energy = 1+lowest-vacuum_lowest
    bounds = list(map(F, receipt['infinite_loop_charged_rest_certificate']['charged_rest_minus_full_field_vacuum_electron_units_interval_exact']))
    assert float(bounds[0]) < energy < float(bounds[1])
    small_states, small_where, small = rotor_matrix(6, -3)
    initial = np.zeros(len(small)); initial[small_where[0, 0]] = 1
    large_initial = np.zeros(len(H)); large_initial[where[0, 0]] = 1
    time = 1.5
    actual_small = expm(-1j*time*small)@initial
    actual_large = expm(-1j*time*H)@large_initial
    embedded = np.zeros_like(actual_large)
    for k, state in enumerate(small_states):
        embedded[where[state]] = actual_small[k]
    error = float(np.linalg.norm(embedded-actual_large))
    z = F(3, 8); r = 3
    from math import factorial
    bound = 2*z**r/factorial(r)/(1-z/(r+1))
    assert error < float(bound)
    for row in receipt['complete_native_executions']:
        cutoff = F(row['infinite_source_exact_Galerkin_original_record_TV_truncation_upper_exact'])
        numerical = F(row['finite_position_flux_numerical_certificate']['original_input_position_flux_table_L1_half_error_upper_exact'])
        total = F(row['infinite_source_original_position_flux_table_L1_half_error_upper_exact'])
        assert total == cutoff+numerical and numerical > 0
    return {'independent_complete_large_cap_dimension': len(H),
            'charged_energy_above_field_vacuum_electron_units': energy,
            'lies_in_full_infinite_rotor_certificate': True,
            'independent_expm_cutoff_vector_difference': error,
            'Dyson_cutoff_upper_exact': str(bound),
            'physical_isolated_electron_or_Yang_Mills_gap_inferred': False}


def pi_certificate():
    source = ast.parse((ROOT/'model/openwave_transverse_completion/source.py').read_text())
    phase = next(x for x in source.body if isinstance(x, ast.FunctionDef) and x.name == 'exact_phase')
    assignment = next(x for x in phase.body if isinstance(x, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'pi' for t in x.targets))
    pi = F(ast.literal_eval(assignment.value.args[0]))
    def atan_bounds(q, terms):
        x = F(1, q)
        partial = sum(((-1)**k*x**(2*k+1)/F(2*k+1) for k in range(terms)), F(0))
        following = (-1)**terms*x**(2*terms+1)/F(2*terms+1)
        return min(partial, partial+following), max(partial, partial+following)
    a, b = atan_bounds(5, 40), atan_bounds(239, 12)
    lo, hi = 16*a[0]-4*b[1], 16*a[1]-4*b[0]
    assert pi < lo < hi < pi+F(1, 10**50)
    return {'identity': 'pi=16 atan(1/5)-4 atan(1/239)',
            'alternating_series_terms': [40, 12],
            'actual_owner_constant_exact': str(pi),
            'declared_constant_error_upper_exact': str(F(1, 10**50)),
            'independent_rational_enclosure_verifies_constant': True}


def complete_work_account(receipt):
    """Rebuild final field energy and detector Gauss charge from all records."""
    preparation = receipt['complete_preparation_and_detector']
    bare_field = F(preparation['input_field_energy_electron_units_exact'])
    assert bare_field == F(1, 8)
    checked = 0
    tolerance = F(1, 500000000000)
    for row in receipt['complete_native_executions']:
        charge = 3*F(row['charge_exact'])
        assert charge.denominator == 1
        charge = charge.numerator
        electric = F(0)
        for (position, loop), probability in zip(row['all_position_flux_labels'],
                                                row['conditional_position_flux_probabilities'], strict=True):
            flux = [loop-charge*(edge < position) for edge in range(4)]
            # After absorption, the local detector retains q and the original
            # compensating reference retains -q. No field record is reset.
            assert all(flux[v]-flux[(v-1)%4] == charge*((v == position)-(v == 0))
                       for v in range(4))
            electric += F(probability)*sum(e*e for e in flux)/32
        assert F(row['preparation_failure_probability_exact']) == F(1, 2)
        assert F(row['preparation_input_and_failure_field_energy_electron_units_exact']) == bare_field
        work = F(row['preparation_mean_source_work_electron_units_diagnostic'])
        work += F(row['destructive_readout_mean_source_work_electron_units_diagnostic'])
        assert abs(work-electric/2) < tolerance
        account = row['complete_energy_account']
        assert F(account['initial_mean_electron_units']) == bare_field
        assert abs(F(account['final_mean_with_failure_electron_units'])-bare_field-electric/2) < tolerance
        checked += 1
    assert checked == 30
    return {'complete_executions_checked': checked,
            'bare_input_and_failure_energy_exact': str(bare_field),
            'final_energy_rebuilt_from_every_position_flux_probability': True,
            'local_detector_and_original_reference_Gauss_charge_preserved': True,
            'displayed_diagnostic_energy_balance_tolerance_exact': str(tolerance),
            'true_field_vacuum_is_not_assumed_to_equal_bare_input': True}


def run():
    flux = json.loads((ROOT/'model/openwave_quark_flux_records/results.json').read_text())
    rotor = json.loads((ROOT/'model/openwave_transverse_completion/results.json').read_text())
    for owner in (flux, rotor):
        assert owner['verdict'].startswith('PASS')
        for path, digest in owner['source_sha256'].items():
            assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == digest, path
    assert flux['preparation_and_failure']['native_pump']['exact_complete_Gram_sum'] == 'I_64'
    return {'schema': 'oph.independent_compact_gauge_completion_review.v1',
            'verdict': 'PASS_INDEPENDENT_COMPACT_GAUGE_CARRIER_AND_RECORD_REVIEW',
            'complete_native_preparation': car_preparation(),
            'independent_SU3_class_moments': group_moments(flux),
            'independent_compact_rotor': rotor_comparison(rotor),
            'independent_phase_range_reduction_constant': pi_certificate(),
            'independent_complete_preparation_and_absorption_work': complete_work_account(rotor),
            'domains': {'full_Gauss_constraints_and_reference_retained': True,
                'fractional_flux_interference_is_a_self_sustained_texture': False,
                'finite_plaquette_dressing_is_physical_electron_pole': False,
                'finite_cutoff_tail_is_entire_numerical_output_error': False}}
