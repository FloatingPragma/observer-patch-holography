"""Exact complete-level Ramsey clocks and an all-transition lapse test."""
from fractions import Fraction as F
from math import factorial
import sympy as s
from model.openwave_native_force.force import source_matrices
from model.finite_area_response.response import Interval, spectral_enclosures
from .geometry import prepared


def sine(value):
    if not isinstance(value, Interval) or max(abs(value.lo),abs(value.hi)) > 8:
        raise ValueError('a certified phase interval in [-8,8] is required')
    square = value*value
    term = value
    total = term
    for k in range(1,41):
        term = -term*square/((2*k)*(2*k+1))
        total = total+term
    tail = max(abs(value.lo),abs(value.hi))**83/F(factorial(83))
    return total+Interval(-tail,tail)


def certificate():
    _, radii, _ = prepared(); locations = (15,23)
    symbol,matrix,_,_ = source_matrices()
    spectra = [spectral_enclosures(matrix.subs(symbol,s.Rational(radii[j]))) for j in locations]
    lapse = radii[15]/radii[23]
    rows = []
    for level in range(1,8):
        gaps = [spectrum[level]-spectrum[0] for spectrum in spectra]
        ratio = gaps[0]/gaps[1]
        mismatch = ratio-lapse
        assert mismatch.lo > 0 or mismatch.hi < 0
        probabilities = [(Interval.point(1)-sine(gap))/2 for gap in gaps]
        assert all(0<=p.lo<=p.hi<=1 and p.hi-p.lo<F(1,10**18) for p in probabilities)
        rows.append({'excited_level':level, 'gap_over_sqrt14_intervals':[gap.encoded() for gap in gaps],
                     'near_to_far_frequency_ratio_interval':ratio.encoded(),
                     'frequency_ratio_decimal_interval':[float(ratio.lo),float(ratio.hi)],
                     'ratio_minus_metric_lapse_interval':mismatch.encoded(),
                     'Ramsey_plus_Y_probability_intervals':[p.encoded() for p in probabilities],
                     'Ramsey_minus_Y_probability':'one minus plus-Y probability',
                     'all_six_other_level_outcomes':'zero, retained'})
    ratios = [Interval(*map(F,row['near_to_far_frequency_ratio_interval'])) for row in rows]
    assert all(a.hi<b.lo or b.hi<a.lo for i,a in enumerate(ratios) for b in ratios[i+1:])
    return {'native_local_Hamiltonian':'unchanged three-oscillator cutoff-two H(r)=sqrt(14)M(r)',
            'near_and_far_slice_indices':list(locations), 'complete_levels_per_clock':8,
            'joint_clock_carrier_dimension':64,
            'Ramsey_duration_in_unchanged_source_time':'1/sqrt(14)',
            'clock_preparation':'normalized ground/excited eigenstate superposition; spectral preparation is declared',
            'common_phase_apparatus':'Y-basis Ramsey readout at both locations',
            'rows':rows, 'single_metric_lapse_ratio_for_every_transition':False,
            'native_clocks_redefined_to_enforce_metric_ratio':False,
            'fixed_coordinate_devices_have_identical_proper_spatial_dimensions':False,
            'clock_comparison_refutes_equivalence_principle':False,
            'energy_transfer_between_clocks_or_universal_gravitational_redshift_proved':False,
            'all_complete_local_levels_and_reference_phase_retained':True}
