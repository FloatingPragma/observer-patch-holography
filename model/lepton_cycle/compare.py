"""Frozen architecture comparison; empirical targets cannot enter the model."""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import mpmath as mp

from . import composition, sensitivity
from tools import strict_json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ACCURACY = F(1, 10**6)  # Declared channel accuracy, independent of mass residuals.


def fmt(x):
    return mp.nstr(x, 25)


def architecture():
    with mp.workdps(60):
        theta, balance = mp.mpf(2)/9, mp.sqrt(2)
        mu, tau = sensitivity.ratios(theta, balance)
        e = [x*x for x in sensitivity.spectrum(theta, balance, mp.sqrt(mp.mpf(5)/108))]
        j = sensitivity.logarithmic_jacobian(theta, balance)
        return {'theta': '2/9', 'balance': 'sqrt(2)',
                'mode_assignment': {'electron': 1, 'muon': 2, 'tau': 0},
                'muon_electron_ratio': fmt(mu), 'tau_electron_ratio': fmt(tau),
                'effect_eigenvalues_k0_k1_k2': list(map(fmt, e)),
                'log_ratio_jacobian_theta_balance': [[fmt(j[k,l]) for l in (0,1)] for k in (0,1)],
                'physical_pole_mass_map_proved': False}


def ratio_box(numerator, denominator, numerator_error, denominator_error):
    """Algebraic enclosure conditional on an input box; no coverage claim."""
    a, b, u, v = map(mp.mpf, (numerator, denominator, numerator_error, denominator_error))
    if not all(mp.isfinite(x) for x in (a,b,u,v)) or u < 0 or v < 0 or a < u or b <= v:
        raise ValueError('positive denominator and nonnegative bounded numerator required')
    return [(a-u)/(b+v), (a+u)/(b-v)]


def preparation_sensitivity(eta):
    """Weyl/variational effect bound for a fixed analyzer and independent source."""
    with mp.workdps(60):
        eta = mp.mpf(eta)
        e = [mp.mpf(x) for x in architecture()['effect_eigenvalues_k0_k1_k2']]
        if eta < 0 or eta >= min(e) or 2*eta >= min(abs(x-y) for i,x in enumerate(e) for y in e[i+1:]):
            raise ValueError('preparation error must preserve positivity and separated mode labels')
        return {'source_pair_trace_distance_bound': fmt(eta),
                'effect_operator_norm_bound': fmt(eta),
                'muon_electron_rate_ratio_interval': list(map(fmt, ratio_box(e[2],e[1],eta,eta))),
                'tau_electron_rate_ratio_interval': list(map(fmt, ratio_box(e[0],e[1],eta,eta))),
                'scope': 'Fixed analyzer; source independent of register/reference. These are limiting effect-rate enclosures, not finite-channel mass-error bars or physical mass uncertainties.'}


def finite_design():
    designs = {kind: composition.implementation_design(ACCURACY, cycle=kind) for kind in ('weak','mass')}
    return {'target_half_diamond_error': str(ACCURACY), 'declared_time': '1',
            'round_selection': 'Ceiling of twice the proved coefficient divided by the declared channel accuracy; no target masses consumed.',
            'designs': designs,
            'preparation_error_examples': [composition.error_budget(designs[k]['rounds'], cycle=k,
                preparation_error=eta) for k in ('weak','mass') for eta in (F(0), F(1,10**12), F(1,10**10))],
            'scope': 'Arbitrary untouched reference, complete channels and independent fresh preparation. Conditional service design, not measured hardware accuracy. Channel distance does not identify an additive physical mass or ratio correction.'}


def comparison(inputs=None):
    inputs = strict_json.load(HERE/'experimental_inputs.json') if inputs is None else inputs
    with mp.workdps(60):
        model = architecture()
        values = {}
        for key in ('muon_electron_ratio', 'electron_MeV', 'tau_MeV'):
            row = inputs[key]
            x,u = mp.mpf(row['value']),mp.mpf(row['standard_uncertainty'])
            if not mp.isfinite(x) or not mp.isfinite(u) or x <= 0 or not 0 < u < x:
                raise ValueError('positive finite reference and standard uncertainty required')
            values[key] = (x,u)
        mu,um = values['muon_electron_ratio']; e,ue = values['electron_MeV']; tau,ut = values['tau_MeV']
        rm,rt = mp.mpf(model['muon_electron_ratio']),mp.mpf(model['tau_electron_ratio'])
        tau_ref = tau/e
        # For d(tau/e) = d(tau)/e - tau*d(e)/e², Cauchy-Schwarz
        # gives the attainable range over all PSD two-variable covariances.
        a,b = ut/e,tau*ue/e**2
        sigma_range = [abs(a-b),a+b]
        delta_tau = rt-tau_ref
        return {'architecture': model, 'experimental_inputs': inputs,
            'muon': {'reference_ratio': fmt(mu), 'reference_standard_uncertainty': fmt(um),
                     'residual_ratio': fmt(rm-mu), 'relative_residual_ppm': fmt((rm/mu-1)*10**6),
                     'reference_only_standardized_discrepancy': fmt((rm-mu)/um)},
            'tau': {'reference_ratio': fmt(tau_ref), 'reference_ratio_linearized_sd_range': list(map(fmt,sigma_range)),
                    'residual_ratio': fmt(delta_tau), 'relative_residual_ppm': fmt((rt/tau_ref-1)*10**6),
                    'reference_only_standardized_discrepancy_range': list(map(fmt,sorted([delta_tau/sigma_range[1],delta_tau/sigma_range[0]]))),
                    'conditional_one_standard_uncertainty_input_box_ratio': list(map(fmt,ratio_box(tau,e,ut,ue))),
                    'electron_calibrated_candidate_MeV': fmt(rt*e), 'candidate_minus_reference_MeV': fmt(rt*e-tau)},
            'uncertainty_contract': {'physical_theory_standard_uncertainty': None, 'joint_chi_squared': None,
                'joint_covariance_available': False, 'sigma_range_method': 'Linearized covariance-compatible PSD/Cauchy bounds using the supplied marginal uncertainties.',
                'interpretation': 'Standardized discrepancies divide by reference uncertainty only. They are not a physical-theory significance or a joint fit. Input boxes have no asserted frequentist coverage.'},
            'implementation': finite_design(),
            'preparation_sensitivity': [preparation_sensitivity(x) for x in ('0.0000000001','0.00000001')],
            'exposure': {'comparison_type': 'conditional retrospective postdiction',
                'parameters_fitted_in_this_comparison': [], 'targets_read_by_architecture': [],
                'target_selected_round_count': False, 'phase_displacement_applied': False,
                'architectural_choices': ['independent fresh source pairs','equal analyzer weights','theta=2/9','balance=sqrt(2)','squared-eigenvalue rate interpretation','electron/muon/tau mode assignment'],
                'empirical_calibrations': ['electron energy only for the explicitly calibrated MeV display'],
                'prior_empirical_exposure': 'Particle masses and Koide relations were known during architecture development; no blind parameter-selection claim.',
                'unproved_physical_maps': ['register effect rates to physical interacting lepton mass poles','declared time to a physical energy clock']}}


def report():
    result = comparison()
    paths = ['compare.py','test_compare.py','COMPARISON.md','experimental_inputs.json',
             'composition.py','sensitivity.py','PROOF.md','results.json']
    result['schema'] = 'oph.lepton_physical_comparison.v1'
    result['source_bindings'] = {str((HERE/p).relative_to(ROOT)):hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths}
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true'); group.add_argument('--check',action='store_true')
    args=parser.parse_args(); result=report(); path=HERE/'physical_comparison.json'
    if args.write: path.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    elif strict_json.load(path)!=result: raise ValueError('lepton comparison differs from replay')
    print(json.dumps({'muon':result['muon'],'tau':result['tau']},indent=2))


if __name__=='__main__': main()
