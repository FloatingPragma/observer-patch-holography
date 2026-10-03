"""Independent electroweak checks for the declared source-scale VEV attachment.

The two-loop running uses squared couplings and Radau, independently of
SMDR's C Runge--Kutta implementation.  The charged-current coefficient is
evaluated from Martin--Robertson, arXiv:1907.02500, in the zero-light-Yukawa
limit.  It is a local Wilson coefficient, not a simulated muon lifetime.
This module never imports the matching producer or an experimental target.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from verify_top_positive_capacity_band import beta as squared_beta

K = 1.0 / (16.0 * math.pi**2)
ROOT = HERE.parents[2]
SPEC = HERE / 'source_ew_vev_matching_spec.json'
OUT = ROOT / 'code/particles/runs/calibration/source_ew_vev_matching.json'
REVIEWED_SPEC_SHA256 = 'edc9be3eba98cb4b06da54854c5ed8c268fd53c5f6dc0dafb1124b4bd551f466'
REVIEWED_C_SOURCE_SHA256 = 'b36f333bf227e069d8ffea2410b793ab5847a49d5e30df2222f07eb074e85805'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def finite_number(value):
    require(isinstance(value, (int, float)) and not isinstance(value, bool),
            'invalid numeric type')
    require(math.isfinite(value), 'nonfinite number')
    return float(value)


def strict_load(path):
    def pairs(items):
        result = {}
        for key, value in items:
            require(key not in result, 'duplicate JSON key')
            result[key] = value
        return result

    def invalid(value):
        raise ValueError('nonfinite JSON constant')

    return json.loads(Path(path).read_text(), object_pairs_hook=pairs,
                      parse_constant=invalid)


def same_structure(actual, expected):
    """Exact typed metadata comparison, including boolean/number separation."""
    require(json.dumps(actual, sort_keys=True, allow_nan=False)
            == json.dumps(expected, sort_keys=True, allow_nan=False),
            'reviewed matching metadata changed')


def near(actual, expected, *, atol=1e-9, rtol=2e-10):
    require(abs(finite_number(actual) - expected) <= atol + rtol * abs(expected),
            'independent EW numerical mismatch')


def higgs_anomalous_dimension(squared):
    """Landau-gauge gamma=-d log(v)/d log(Q), through two loops.

    Top-only specialization of the Higgs anomalous dimension used by SMDR,
    with gY rather than GUT-normalized hypercharge.  Its two-loop gauge and
    Yukawa terms are retained; treating the source VEV as constant fails
    this check even when the other five running parameters are correct.
    """
    a, b, c, h, lam = squared
    gamma1 = 3*h - 9*b/4 - 3*a/4
    gamma2 = (-271*b*b/32 + 9*a*b/16 + 431*a*a/96 + 6*lam*lam
              + 45*b*h/8 + 20*c*h + 85*a*h/24 - 27*h*h/4)
    return K*gamma1 + K*K*gamma2


def running_parameters(inputs, scale):
    """Propagate one fixed boundary; do not reidentify v at the new scale."""
    q0 = finite_number(inputs['Q0_GeV'])
    v0 = finite_number(inputs['v_Q0_GeV'])
    scale = finite_number(scale)
    require(q0 > 0 and v0 > 0 and scale > 0, 'positive scale and VEV required')
    require(inputs['non_top_yukawas'] == 0.0, 'top-only specialization required')
    require(inputs['RG_loop_order'] == 2, 'two-loop running required')
    state = [finite_number(inputs[key])**2 for key in
             ('gY_Q0', 'g2_Q0', 'g3_Q0', 'yt_Q0')]
    state += [finite_number(inputs['lambda_Q0']), math.log(v0)]

    def rhs(t, values):
        return np.r_[squared_beta(t, values[:5]),
                     -higgs_anomalous_dimension(values[:5])]

    if scale == q0:
        result = np.asarray(state)
    else:
        solution = solve_ivp(rhs, (math.log(q0), math.log(scale)), state,
                             method='Radau', rtol=2e-12, atol=2e-14)
        require(solution.success, 'independent RG integration failed')
        result = solution.y[:, -1]
    require(all(result[:4] > 0), 'nonpositive squared running coupling')
    return dict(zip(('gY', 'g2', 'g3', 'yt'), map(math.sqrt, result[:4])),
                **{'lambda': float(result[4]), 'v': math.exp(result[5]), 'Q': scale})


def charged_current_coefficient(parameters, loop_order):
    """GF= (1+Delta-r-tilde)/(sqrt(2)*v²), at zero or one loop.

    The finite empirical muon-mass correction in the stock SMDR routine is
    absent in this zero-muon theory.  No hadronic alpha input is required.
    A(x)=x(log(x/Q²)-1) is the finite tadpole integral convention.
    """
    require(type(loop_order) is int and loop_order in (0, 1),
            'independent GF check supports zero/one loop')
    gp, g, yt, lam, v, q = [finite_number(parameters[key]) for key in
                            ('gY', 'g2', 'yt', 'lambda', 'v', 'Q')]
    require(min(gp, g, yt, lam, v, q) > 0, 'positive matching parameters required')
    tree = 1 / (math.sqrt(2)*v*v)
    if loop_order == 0:
        return tree
    w = g*g*v*v/4
    z = (g*g + gp*gp)*v*v/4
    h = 2*lam*v*v
    t = yt*yt*v*v/2
    require(abs(w-h) > 1e-10*max(w, h), 'degenerate W/H limit not implemented')
    require(abs(z-w) > 1e-10*max(z, w), 'degenerate Z/W limit not implemented')

    def tadpole(x):
        return x*(math.log(x/(q*q))-1)

    delta = (3*yt*yt*tadpole(t)/t
             + 3*(g*g-gp*gp)*(tadpole(z)-tadpole(w))/(4*(z-w))
             + (3*g*g*tadpole(h)/4 + 3*(6*lam-g*g)*tadpole(w))/(w-h)
             - 3*g*g/8 - gp*gp/8 - lam + 3*yt*yt/2)
    return tree*(1 + K*delta)


def verify_gf_and_running(inputs, forward):
    """Check every scheduled scale/order, including finite-muon removal.

    This independently verifies the zero/one-loop charged-current values and
    two-loop transport.  Two-loop matching requires the separate pinned C
    replay; this function does not reclassify it as independently calculated.
    """
    require(inputs['scale_multipliers'] == [1, 2, 4], 'scale controls changed')
    rows = forward['rows']
    inventory = [(row['scale_multiplier'], row['matching_order']) for row in rows]
    require(inventory == [(scale, order) for scale in (1, 2, 4)
                          for order in (0, 1, 2)], 'incomplete matching controls')
    expected = {scale: running_parameters(inputs, scale*inputs['Q0_GeV'])
                for scale in (1, 2, 4)}
    for row in rows:
        working = row['working']
        ref = expected[row['scale_multiplier']]
        for key, value in ref.items():
            near(working[key], value)
        near(working['mt_MSbar'], ref['yt']*ref['v']/math.sqrt(2))
        removed = 0.00000051862/(math.sqrt(2)*ref['v']**2)
        near(row['GF_removed_finite_muon_term_GeVm2'], removed, atol=1e-21)
        near(row['GF_local_zero_muon_GeVm2'],
             row['GF_raw_GeVm2'] - removed, atol=1e-18, rtol=1e-12)
        require(finite_number(row['GF_local_zero_muon_GeVm2']) > 0,
                'positive local charged-current coefficient required')
        order = row['matching_order']
        if order <= 1:
            near(row['GF_local_zero_muon_GeVm2'],
                 charged_current_coefficient(ref, order), atol=1e-16, rtol=2e-10)
        if order == 0:
            tree = {'top_method0': ref['yt']*ref['v']/math.sqrt(2),
                    'top_method1': ref['yt']*ref['v']/math.sqrt(2),
                    'higgs': math.sqrt(2*ref['lambda'])*ref['v'],
                    'W': ref['g2']*ref['v']/2,
                    'Z': math.hypot(ref['g2'], ref['gY'])*ref['v']/2}
            for name, mass in tree.items():
                near(row[name]['mass_GeV'], mass)
                near(row[name]['width_GeV'], 0, atol=1e-12, rtol=0)
        for name in ('W', 'Z'):
            particle = row[name]
            mass = finite_number(particle['mass_GeV'])
            width = finite_number(particle['width_GeV'])
            # SMDR returns sqrt(s)=M-i Gamma/2.  PDG running-width BW uses
            # Re(s)=M_BW²/(1+Gamma_BW²/M_BW²), not this same M.
            re_s = mass*mass-width*width/4
            require(re_s > 0, 'positive real squared pole required')
            old_width = mass*width/math.sqrt(re_s)
            bw_mass = math.sqrt(re_s + old_width*old_width)
            near(particle['BW_mass_GeV'], bw_mass)
            near(particle['BW_width_GeV'], old_width*bw_mass/math.sqrt(re_s))
    return {'transport_scales_checked': 3, 'matching_rows_checked': 9,
            'independent_GF_orders': [0, 1],
            'independent_two_loop_matching_claimed': False,
            'physical_source_attachment_derived': False}


def verify(receipt):
    """Certify source ancestry and independent lower-order consequences.

    The retained multiloop payload is byte-bound, and must additionally be
    reproduced by the pinned upstream C implementation.  Its digest is not
    treated as a proof of the two-loop integrals.
    """
    require(hashlib.sha256(SPEC.read_bytes()).hexdigest() == REVIEWED_SPEC_SHA256,
            'unreviewed matching specification')
    spec = strict_load(SPEC)
    require(set(receipt) == {'schema', 'specification', 'scope', 'inputs',
                             'source_pins', 'forward'}, 'matching receipt keys')
    same_structure(receipt['schema'], 'oph.source_ew_vev_matching.v1')
    same_structure(receipt['specification'], spec)
    same_structure(receipt['scope'], spec['scope'])
    source = ROOT / spec['source_receipt']
    paths = [HERE/'source_ew_vev_matching.py', SPEC,
             HERE/'source_ew_vev_matching_smdr.c',
             HERE/'sm_two_loop_rge_engine.py', source]
    pins = {path.relative_to(ROOT).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in paths}
    same_structure(receipt['source_pins'], pins)
    require(pins['code/particles/calibration/source_ew_vev_matching_smdr.c']
            == REVIEWED_C_SOURCE_SHA256, 'unreviewed C matching driver')
    prior = strict_load(source)
    scales = prior['source_scales_GeV']
    branch = prior['fully_coupled_two_loop_same_low_gauge_anchor'][1]
    uv = finite_number(scales['log_midpoint_half_turn'])
    q0 = finite_number(scales['mz_run_gev'])
    v0 = finite_number(scales['v_transmutation_gev'])
    gauges = [finite_number(x) for x in branch['boundary_gauges']]
    top = finite_number(branch['boundary_top_yukawa'])
    require(len(gauges) == 3 and min(gauges) > 0 and top > 0 and uv > q0 > 0,
            'invalid source boundary')
    initial = [3*gauges[0]**2/5, gauges[1]**2, gauges[2]**2, top*top, 0.0]
    flow = solve_ivp(squared_beta, (math.log(uv), math.log(q0)), initial,
                     method='Radau', rtol=2e-12, atol=2e-14)
    require(flow.success, 'independent UV transport failed')
    a, b, c, h, lam = flow.y[:, -1]
    numeric = {'Q0_GeV': q0, 'v_Q0_GeV': v0, 'gY_Q0': math.sqrt(a),
               'g2_Q0': math.sqrt(b), 'g3_Q0': math.sqrt(c),
               'yt_Q0': math.sqrt(h), 'lambda_Q0': float(lam)}
    metadata = {'schema': 'oph.smdr_conditional_matching_input.v1',
                'non_top_yukawas': 0.0, 'm2_minimum_order': 2,
                'RG_loop_order': 2, 'scale_multipliers': [1, 2, 4],
                'SMDR_commit': spec['upstream']['commit'],
                'attachment': spec['hypotheses']['vev_attachment'],
                'criticality_input': {'scale_GeV': uv,
                                     'initial': [*gauges, top, 0.0]},
                'source_file': spec['source_receipt'],
                'consumed_fields': spec['source_boundary_fields']}
    inputs = receipt['inputs']
    require(set(inputs) == set(numeric) | set(metadata), 'input inventory changed')
    for key, value in numeric.items():
        near(inputs[key], value, atol=2e-11, rtol=2e-10)
    for key, value in metadata.items():
        same_structure(inputs[key], value)
    forward = receipt['forward']
    require(set(forward) == {'schema', 'm2_minimum_controls', 'rows'},
            'forward payload inventory changed')
    same_structure(forward['schema'], 'oph.smdr_forward_matching.v1')
    digest = hashlib.sha256(json.dumps(forward, sort_keys=True,
                            separators=(',', ':'), allow_nan=False).encode()).hexdigest()
    require(digest == spec['retained_forward_payload']['sha256'],
            'retained multiloop payload changed')
    result = verify_gf_and_running(inputs, forward)
    from verify_source_ew_top_pole import verify as verify_top
    require(verify_top(receipt), 'independent top/tadpole check failed')
    return dict(result, source_boundary_reconstructed=True,
                full_multiloop_replay_required=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('receipt', type=Path, nargs='?', default=OUT)
    args = parser.parse_args()
    print(json.dumps(verify(strict_load(args.receipt)), sort_keys=True))
