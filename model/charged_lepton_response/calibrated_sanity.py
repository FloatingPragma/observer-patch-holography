"""Explicit one-ratio phase calibration, with tau held out of the fit."""
from fractions import Fraction as F
from pathlib import Path
import json
import math

import numpy as np
import sympy as s

from model.native_flavor_orbit import intervals as iv
from model.standard_model_source_map import koide_source
from model.common_source import native_record_yukawa_field as field
from . import relativistic_mass

REFERENCE = Path(__file__).with_name('calibration_reference.json')


def _sqrt(box):
    return iv.I(iv.sqrt_rational(box.lo).lo, iv.sqrt_rational(box.hi).hi)


def choose_phase(muon_electron_ratio):
    """Use the muon/electron ratio alone; no tau, scale or detector data."""
    if type(muon_electron_ratio) not in (str, int, F):
        raise ValueError('an exact calibration ratio is required')
    r = F(muon_electron_ratio)
    if r <= 1:
        raise ValueError('the calibrated ordered ratio must exceed one')
    q = _sqrt(iv.I(r))
    u = (q-1)/(q+1)
    tangent = u*(_sqrt(iv.I(2))-1)/(_sqrt(iv.I(3))+_sqrt(3-u*u))
    # A readable rational unit-circle phase, then independently certified below.
    denominator = 10**18
    t = F(round((tangent.lo+tangent.hi)*denominator/2), denominator)
    x, y = (1-t*t)/(1+t*t), 2*t/(1+t*t)
    if x*x+y*y != 1:
        raise ArithmeticError('the selected phase is not on the exact unit circle')
    return {'half_angle_tangent': str(t), 'real': str(x), 'imaginary': str(y),
            'phase_radians_approx': 2*math.atan(float(t)),
            'calibration_input_muon_electron_ratio': str(r),
            'tau_used_to_choose_phase': False, 'phase_derived_from_repair': False}


def phase_value(phase):
    return s.Rational(phase['real'])+s.I*s.Rational(phase['imaginary'])


def exact_coefficient_boxes(phase):
    """Outward exact radical enclosures for the same circulant's eigenvalues."""
    x, y = F(phase['real']), F(phase['imaginary'])
    if x*x+y*y != 1 or y <= 0:
        raise ValueError('a positive exact unit-circle phase is required')
    a = _sqrt(iv.I(5))/6
    bx = _sqrt(iv.I(10))*x/12
    dy = _sqrt(iv.I(30))*y/12
    roots = [a-bx-dy, a-bx+dy, a+2*bx]
    if roots[0].lo <= 0 or any(left.hi >= right.lo for left, right in zip(roots, roots[1:])):
        raise ValueError('the ordered strictly positive chamber is required')
    masses = [root**2/3 for root in roots]
    ratios = [iv.I(1), masses[1]/masses[0], masses[2]/masses[0]]
    return masses, ratios


def conditional_tau_ratio(muon_interval):
    r = iv.I(*map(F, muon_interval))
    if r.lo <= 1:
        raise ValueError('ordered positive calibration data required')
    q = _sqrt(r)
    return (2*(1+q)+_sqrt(3*(1+4*q+r)))**2


def compare(reference, phase):
    """The held-out data can change the comparison, never the chosen source phase."""
    masses, ratios = exact_coefficient_boxes(phase)
    windows = {}
    for key in ('muon_electron_ratio', 'tau_electron_ratio'):
        datum = reference[key]
        center, uncertainty = F(datum['value']), F(datum['standard_uncertainty'])
        if uncertainty <= 0:
            raise ValueError('positive comparison uncertainty required')
        windows[key] = iv.I(center-uncertainty, center+uncertainty)
    mu, tau = ratios[1:]
    electron = F(reference['electron_MeV']['value'])
    tau_mev = tau*electron
    return {'source_action_coefficients': [x.strings() for x in masses],
            'source_action_ratios': [x.strings() for x in ratios],
            'calibrated_mu_ratio_inside_reported_one_sigma_window':
                windows['muon_electron_ratio'].lo <= mu.lo <= mu.hi <= windows['muon_electron_ratio'].hi,
            'held_out_tau_ratio_inside_reported_one_sigma_window':
                windows['tau_electron_ratio'].lo <= tau.lo <= tau.hi <= windows['tau_electron_ratio'].hi,
            'conditional_tau_ratio_from_mu_one_sigma_box':
                conditional_tau_ratio(windows['muon_electron_ratio'].strings()).strings(),
            'optional_calibrated_mass_MeV': [(x*electron).strings() for x in ratios],
            'common_action_energy_factor_MeV': (iv.I(electron)/masses[0]).strings(),
            'physical_loop_correction_bounds_supplied': False,
            'joint_statistical_significance_claimed': False,
            'comparison_proves_physical_lepton_identity': False,
            'tau_reference_consumed_by_source_phase_or_kernel': False,
            'tau_MeV_approx': float((tau_mev.lo+tau_mev.hi)/2)}


def certificate():
    reference = json.loads(REFERENCE.read_text())
    phase = choose_phase(reference['muon_electron_ratio']['value'])
    z = phase_value(phase)
    source = koide_source.exact_source(z)
    if not source['physical_action_Koide_ratio_certified']:
        raise ArithmeticError('calibration escaped the positive Koide chamber')
    k, y = relativistic_mass.source_kernel(z)
    controls = field.dispersion_certificate(np.array(k.evalf(17), complex), speed=1)
    if controls['dispersion_spectral_error'] > 1e-12:
        raise ArithmeticError('the calibrated phase did not survive the native field kernel')
    result = compare(reference, phase)
    if not result['calibrated_mu_ratio_inside_reported_one_sigma_window']:
        raise ArithmeticError('phase calibration failed its declared precision')
    return {'classification': 'explicitly calibrated free-action sanity test',
            'phase': phase, 'exact_positive_root_Koide_ratio': '2/3',
            'same_native_180_dimensional_free_field_control': controls,
            'comparison': result, 'reference': reference,
            'inputs': ['muon/electron ratio for phase', 'electron MeV for optional common scale'],
            'held_out_quantity': 'tau/electron ratio',
            'phase_calibration_is_source_derivation': False,
            'free_action_is_interacting_measured_particle_theory': False,
            'physical_particle_assignments': []}
