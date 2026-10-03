"""Conditional top coordinate from independently fixed P and declared criticality.

No particle-mass target, flavor-dependent P closure, or quark template is
consumed. The external energy unit, D10 gauge recipe, criticality boundary,
transport truncations and incomplete physical matching are explicit in the
adjacent specification. Numerical convergence is not a theory error bound.
"""
from __future__ import annotations

import argparse
from contextlib import ExitStack
from dataclasses import asdict
from decimal import Decimal
import hashlib
import json
import math
from pathlib import Path
import sys
from unittest.mock import patch

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, root

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SPEC_PATH = HERE / 'conditional_top_fixed_p_spec.json'
OUT_PATH = HERE.parent / 'runs/calibration/conditional_top_fixed_p.json'
sys.path.insert(0, str(ROOT / 'code/P_derivation'))
sys.path.insert(0, str(HERE))
from paper_math import PaperMathContext
import derive_d11_criticality_boundary_scan as legacy
import sm_two_loop_rge_engine as rge

PIN_PATHS = (
    'code/P_derivation/codata_2022_alpha_fixture.json',
    'code/P_derivation/paper_math.py',
    'code/particles/calibration/derive_d11_criticality_boundary_scan.py',
    'code/particles/calibration/sm_two_loop_rge_engine.py',
    'code/particles/calibration/conditional_top_fixed_p_spec.json',
    'code/particles/calibration/conditional_top_fixed_p.py',
)
FORBIDDEN_METHODS = (
    'alpha_external_from_d10', 'structured_thomson_running',
    'structured_thomson_running_asymptotic', 'diagonal_quark_masses',
    'charged_lepton_masses',
)


def forbidden(*args, **kwargs):
    raise AssertionError('Flavor-dependent P closure or mass template was called')


def source_pins():
    return {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in PIN_PATHS}


def source_packet(spec):
    fixture = json.loads((ROOT / spec['alpha_fixture']).read_text())
    alpha = fixture['inverse_fine_structure_constant']
    if (alpha['value'] != spec['inverse_alpha'] or
            alpha['standard_uncertainty'] != spec['inverse_alpha_standard_uncertainty']):
        raise ValueError('Independent alpha fixture differs from specification')
    d = spec['D10']
    with ExitStack() as stack:
        for name in FORBIDDEN_METHODS:
            if hasattr(PaperMathContext, name):
                stack.enter_context(patch.object(PaperMathContext, name, forbidden))
        # The constructor ordinarily initializes an unrelated flavor helper.
        # Suppress even that unused initialization on this independent-P route.
        stack.enter_context(patch.object(PaperMathContext, '_derive_stage5_integer_vectors', lambda self: {}))
        ctx = PaperMathContext(precision=d['precision'], su2_cutoff=d['su2_cutoff'], su3_cutoff=d['su3_cutoff'])
        p = ctx.p_from_inverse_alpha(Decimal(alpha['value']))
        packet = ctx.build_d10_from_p(p)
    return p, packet


def critical_roots(gauges):
    """All real roots in h=yt^2, at lambda=0, of the declared two-loop beta."""
    g1, g2, g3 = gauges
    a, b, c = g1*g1, g2*g2, g3*g3
    k = 1/(16*math.pi**2)
    coefficients = [30*k*k,
        -6*k + k*k*(-32*c-8*a/5),
        k*k*(-9*b*b/4-171*a*a/100+63*a*b/10),
        k*3*(2*b*b+(b+3*a/5)**2)/8 + k*k*(305*b**3/16-3411*a**3/2000-289*b*b*a/80-1677*b*a*a/400)]
    roots = np.roots(coefficients)
    if any(abs(z.imag) > 1e-9 for z in roots):
        raise ValueError('Criticality root inventory changed')
    return sorted(float(z.real) for z in roots)


def critical_yt(gauges, bracket):
    f = lambda y: rge.beta_2loop(*gauges, y, 0.)[4]
    lo, hi = bracket
    if f(lo)*f(hi) >= 0:
        raise ValueError('Perturbative critical root is not bracketed')
    y = brentq(f, lo, hi, xtol=1e-14)
    positive = [h for h in critical_roots(gauges) if h > 0]
    if not positive or abs(y*y-positive[0]) > 1e-10:
        raise ValueError('Criticality selected a nonperturbative root')
    return y


def build_payload():
    spec = json.loads(SPEC_PATH.read_text())
    p, d10 = source_packet(spec)
    E = float(spec['E_star_display_GeV'])
    packet = {key: str(value) for key, value in asdict(d10).items()}
    with ExitStack() as stack:
        for name, value in {
            'E_STAR_DISPLAY_GEV': E, 'ALPHA_Y_MZ': float(d10.alpha_y_mz),
            'ALPHA_2_MZ': float(d10.alpha2_mz), 'ALPHA_3_MZ': float(d10.alpha3_mz),
            'MU_Z_OVER_E': float(d10.mz_run),
        }.items():
            stack.enter_context(patch.object(legacy, name, value))
        scales = legacy.source_scales(float(p), float(d10.alpha_u))
        boundary = scales['log_midpoint_half_turn']
        v = scales['v_transmutation_gev']
        mz = scales['mz_run_gev']
        hybrid = {}
        for loops in (1, 2):
            rows = []
            for steps in spec['transport']['hybrid_RK4_steps']:
                row = legacy.run_boundary(boundary, v, mz, n_steps=steps, loops=loops)
                row['top_QCD_pole_proxy_GeV'] = row.pop('mt_pole_gev')
                row['top_running_mass_coordinate_GeV'] = row.pop('mt_msbar_gev')
                row['Higgs_tree_proxy_GeV'] = row.pop('mh_tree_gev')
                row['top_running_coordinate_over_Estar'] = row['top_running_mass_coordinate_GeV']/E
                row['top_QCD_pole_proxy_over_Estar'] = row['top_QCD_pole_proxy_GeV']/E
                row['RK4_steps'] = steps
                rows.append(row)
            hybrid[str(loops)] = rows
        low = np.sqrt(4*math.pi*np.array([5*float(d10.alpha_y_mz)/3, float(d10.alpha2_mz), float(d10.alpha3_mz)]))
        guess = np.array(legacy.gauge_couplings(boundary, mz))*[math.sqrt(5/3), 1, 1]

        def integrate(g, tolerance):
            y0 = np.r_[g, critical_yt(g, spec['boundary']['yukawa_bracket']), 0.]
            solution = solve_ivp(lambda t, y: rge.beta_2loop(*y),
                (math.log(boundary), math.log(50.)), y0,
                rtol=tolerance, atol=tolerance*spec['transport']['coupled_atol_over_rtol'],
                dense_output=True, method=spec['transport']['coupled_method'])
            if not solution.success:
                raise ValueError('Coupled integration failed')
            return solution

        coupled = []
        for tolerance in spec['transport']['coupled_rtol']:
            def residual(g):
                return integrate(g, tolerance).sol(math.log(mz))[:3]-low
            shot = root(residual, guess, tol=1e-10)
            if not shot.success or np.max(abs(residual(shot.x))) > 1e-9:
                raise ValueError('Low-anchor gauge shooting failed')
            guess = shot.x
            sol = integrate(shot.x, tolerance)
            fn = lambda mu: sol.sol(math.log(mu))[3]*v/math.sqrt(2)-mu
            mass = brentq(fn, 50., v, xtol=1e-11)
            _, _, g3, yt, lam = sol.sol(math.log(mass))
            if lam <= 0:
                raise ValueError('Higgs tree proxy is not real and positive')
            alpha_s = g3*g3/(4*math.pi)
            proxy = legacy.top_pole_from_msbar(mass, alpha_s, n_l=spec['readout']['n_l'])
            coupled.append({
                'relative_tolerance': tolerance, 'boundary_gauges': shot.x.tolist(),
                'low_anchor_gauge_residual': (sol.sol(math.log(mz))[:3]-low).tolist(),
                'boundary_top_yukawa': critical_yt(shot.x, spec['boundary']['yukawa_bracket']),
                'critical_polynomial_roots_yt_squared': critical_roots(shot.x),
                'top_running_mass_coordinate_GeV': mass, 'top_QCD_pole_proxy_GeV': proxy,
                'top_running_coordinate_over_Estar': mass/E,
                'top_QCD_pole_proxy_over_Estar': proxy/E,
                'Higgs_tree_proxy_GeV': math.sqrt(2*lam)*v,
                'alpha_s_self_scale': alpha_s, 'lambda_self_scale': lam,
                'yukawa_self_scale': yt, 'top_fixed_point_residual_GeV': fn(mass),
                'beta_lambda_boundary': rge.beta_2loop(*shot.x, critical_yt(shot.x, spec['boundary']['yukawa_bracket']), 0.)[4],
            })
    return {
        'schema': 'oph.conditional_top_fixed_p.v1', 'specification': spec,
        'scope': spec['scope'], 'source_pins': source_pins(),
        'P_from_independent_alpha': str(p), 'D10': packet,
        'source_scales_GeV': scales, 'hybrid_rows': hybrid,
        'fully_coupled_two_loop_same_low_gauge_anchor': coupled,
        'numerical_diagnostics': {
            'hybrid_top_step_differences_GeV': {key: abs(rows[1]['top_running_mass_coordinate_GeV']-rows[0]['top_running_mass_coordinate_GeV']) for key, rows in hybrid.items()},
            'coupled_top_tolerance_difference_GeV': abs(coupled[1]['top_running_mass_coordinate_GeV']-coupled[0]['top_running_mass_coordinate_GeV']),
            'hybrid_to_coupled_top_shift_GeV': hybrid['2'][-1]['top_running_mass_coordinate_GeV']-coupled[-1]['top_running_mass_coordinate_GeV'],
            'interpretation': spec['uncertainty'],
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUT_PATH)
    args = parser.parse_args()
    payload = build_payload()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+'\n')
    print(json.dumps({'receipt': str(args.output), 'coupled_top_coordinate_GeV': payload['fully_coupled_two_loop_same_low_gauge_anchor'][-1]['top_running_mass_coordinate_GeV']}))


if __name__ == '__main__':
    main()
