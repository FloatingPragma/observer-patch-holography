"""Exact weak-response dilution by complete independent geometry registers."""
from fractions import Fraction as F
from math import comb
import json
from pathlib import Path
from model.finite_area_response.response import Interval
from model.openwave_replicated_geometry.source import power, sum_weights
from model.openwave_universal_energy_geometry.source import QSIZE as Q, OMEGA, sinhcosh, profile
from model.openwave_gravity_unit_constraint.source import pi_interval

ROOT = Path(__file__).resolve().parents[2]


def validate_n(n):
    if type(n) is not int or n < 1:
        raise ValueError('a positive integer register count is required')


def susceptibility(n, f):
    """Complete binomial spectral sum, used only for small executable controls."""
    validate_n(n); f = F(f)
    if not 0 <= f <= 1:
        raise ValueError('profile must be in [0,1]')
    sa, ca = sinhcosh(Q/n); sb, cb = sinhcosh(Q*f/n)
    total, rows = sum_weights(n, sa/ca*sb/cb)
    return -2*power(ca,n)*total/OMEGA, rows


def coefficient_bounds(n):
    """All N, including enormous counts, without allocating their carrier.

    tanh(t)>=t-t^3/3, log(cosh(t))<=t^2/2 and exp(z)<=1/(1-z)
    give exact rational bounds. This is a theorem evaluation, not a 2^N replay.
    """
    validate_n(n)
    stiffness = profile()['stiffness_exact']
    plo, phi = pi_interval()
    prefactor = Q*Q*stiffness/(2*OMEGA*n)
    lower = 1-Q*Q/(3*n*n)
    upper = 1/(1-Q*Q/(2*n))
    return prefactor*lower/phi, prefactor*upper/plo


def inherited_units():
    old = json.loads((ROOT/'model/openwave_gravity_unit_constraint/results.json').read_text())
    units = old['external_reference_inputs']
    strong = json.loads((ROOT/'model/openwave_string_positive_cone/results.json').read_text())['common_length_attachment']
    return F(strong['electron_unit_GeV']), F(units['hbar_c_GeV_fm']), F(units['Newton_G_natural_GeV_minus2']), units, F(strong['a_fm'])


def calibrated_count(u, target):
    """One declared empirical hardware calibration, never an independent G prediction."""
    u, target = F(u), F(target)
    if min(u, target) <= 0:
        raise ValueError('positive energy anchor and target required')
    plo, phi = pi_interval()
    prefactor = Q*Q*profile()['stiffness_exact']/(2*OMEGA)
    lo, hi = prefactor/(phi*target*u*u), prefactor/(plo*target*u*u)
    nearest = lambda v: (v+F(1,2)).numerator // (v+F(1,2)).denominator
    assert nearest(lo) == nearest(hi), 'the inherited pi interval does not identify one integer'
    n = nearest(lo); validate_n(n)
    return n, (lo, hi)


def run():
    u, hc, target, references, strong_length = inherited_units()
    n, ideal = calibrated_count(u, target)
    glo, ghi = coefficient_bounds(n)
    Glo, Ghi = glo/u**2, ghi/u**2
    relative = max(abs(Glo/target-1), abs(Ghi/target-1))
    assert relative < F(1,10**40)
    mean_z = Q*Q/(2*n); variance_z = Q*Q/n
    examples = []
    for count in (1,2,4,8,16):
        slope, weights = susceptibility(count, F(1,2))
        probabilities = [F(comb(count,k),2**count) for k in range(count+1)]
        assert sum(probabilities) == 1
        assert sum(comb(count,r) for r in range(count+1)) == 2**count
        examples.append({'registers': count, 'full_geometry_dimension': 2**count,
            'complete_free_spectrum': [{'energy_exact': str(OMEGA*r), 'multiplicity': comb(count,r)} for r in range(count+1)],
            'complete_coordinate_probabilities_exact': [str(p) for p in probabilities],
            'log_lapse_derivative_at_profile_one_half_interval': slope.encoded(),
            'complete_excitation_response_weights': weights,
            'Newton_coefficient_interval_exact': [str(x) for x in coefficient_bounds(count)]})
    return {'schema': 'oph.collective_calibrated_gravity_strength.v1',
        'verdict': 'PASS_CALIBRATED_COMPLETE_SOURCE_WEAK_STRENGTH_AND_UNIT_COMPATIBILITY',
        'architecture': {'Q_N': '(q/N) sum_i Z_i', 'H_geometry': '(omega/2) sum_i(I-X_i)',
            'H_complete': 'H_geometry+exp(Q_N)K_source+exp(F Q_N)K_probe',
            'q_exact': str(Q), 'per_register_frequency_exact': str(OMEGA),
            'native_matter_factors': 'all eight source and all eight probe levels of the retained positive native energy matrix',
            'positive_terms_and_complete_carrier': True,
            'difference_from_previous_member': 'per-register frequency is fixed; the earlier frequency normalization preserving large gravity strength is not used',
            'free_geometry_gap_exact': str(OMEGA)},
        'all_N_theorem': {'complete_geometry_dimension': '2^N',
            'complete_free_spectrum': 'energy r*omega with multiplicity binomial(N,r), every r=0,...,N',
            'coordinate_readout': 'q*(2k-N)/N with probability binomial(N,k)/2^N, every k=0,...,N',
            'log_lapse_derivative': '-(2/omega)*cosh(q/N)^N*sum_{r=1}^N binomial(N,r)*(tanh(q/N)*tanh(qF/N))^r/r',
            'exact_far_field_coefficient': 'g_N=S*q*cosh(q/N)^N*tanh(q/N)/(2*pi*omega)',
            'rational_enclosure': 'g0/N*(1-q^2/(3*N^2)) <= g_N <= g0/N/(1-q^2/(2*N)); g0=S*q^2/(2*pi*omega)',
            'all_profile_relative_error': '1-q^2*F^2/(3*N^2) <= chi_N(F)/(chi_N_prime(0)*F) <= exp(q^2*F/N)',
            'analytic_complete_sector_certificate': True, 'huge_carrier_enumerated_or_simulated': False},
        'calibrated_member': {'registers': n, 'ideal_leading_order_count_interval_exact': [str(x) for x in ideal],
            'calibration_input': 'the inherited observed Newton coefficient chooses N by the stated nearest-integer rule',
            'electron_energy_unit_GeV_exact': str(u), 'electron_anchor_changed': False,
            'free_symbol_speed_matched_length_fm_exact': str(hc/u),
            'free_symbol_speed_ratio_exact': '1',
            'dimensionless_Newton_interval_exact': [str(glo),str(ghi)],
            'conditional_Newton_GeV_minus2_interval_exact': [str(Glo),str(Ghi)],
            'target_Newton_GeV_minus2_exact': str(target),
            'relative_target_discrepancy_upper_exact': str(relative),
            'relative_target_discrepancy_less_than_exact': '1/10^40',
            'vacuum_lapse_mean_bias_upper_exact': str(mean_z/(1-mean_z)),
            'vacuum_relative_lapse_variance_upper_exact': str(variance_z/(1-variance_z)),
            'nonselective_coordinate_readout_work_source_units_exact': str(n*OMEGA/2),
            'nonselective_coordinate_readout_work_GeV_exact': str(n*OMEGA*u/2),
            'supplied_blank_record_capacity': 'enough memory to retain all N+1 coordinate outcomes; no free erasure claim',
            'diagnostic_float': {'registers': float(n), 'conditional_G_GeV_minus2': float((Glo+Ghi)/2),
                'speed_matched_length_fm': float(hc/u), 'relative_target_discrepancy_upper': float(relative),
                'readout_work_GeV': float(n*OMEGA*u/2)}},
        'unchanged_strong_length_control': {'length_fm_exact': str(strong_length),
            'free_symbol_speed_ratio_exact': str(u*strong_length/hc),
            'matches_physical_light_speed': False,
            'different_spatial_unit_not_silently_identified': True},
        'complete_small_register_controls': examples,
        'external_reference_inputs': {**references,
            'input_role': 'Inherited external constants: the Newton value explicitly calibrates N in this new member; hbar c converts the units.'},
        'scope': {'new_empirical_hardware_calibration': True, 'independent_Newton_prediction': False,
            'calibrated_weak_response_and_free_symbol_unit_equations_compatible': True,
            'full_dressed_photon_dispersion_or_finite_source_force_identification': False,
            'common_strong_geometry_Hamiltonian_and_spatial_refinement': False,
            'physical_graviton_or_nonlinear_Einstein_theory': False,
            'previous_fixed_source_failure_overturned': False, 'new_Nature_validation': False}}
