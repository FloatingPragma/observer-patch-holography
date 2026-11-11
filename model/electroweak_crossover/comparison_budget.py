"""Covariance-compatible delta-method errors for the shared EW reference.

The extrema are exact for the serialized rational coefficients and supplied
marginal variance estimates. They are not nonlinear confidence intervals.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

from . import reduction as rd

HERE = Path(__file__).resolve().parent


def rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, float, str, F)):
        raise ValueError('finite rational-compatible scalar required')
    try:
        return F(str(value)) if isinstance(value, float) else F(value)
    except (ValueError, OverflowError, ZeroDivisionError) as exc:
        raise ValueError('finite rational-compatible scalar required') from exc


def solve(matrix, rhs):
    """Exact elimination for the small fixed-design normal equations."""
    n = len(rhs)
    if len(matrix) != n or any(len(row) != n for row in matrix):
        raise ValueError('square linear system required')
    a = [[rational(x) for x in row] + [rational(y)] for row, y in zip(matrix, rhs)]
    for j in range(n):
        pivot = next((i for i in range(j, n) if a[i][j]), None)
        if pivot is None:
            raise ValueError('rank-deficient fit design')
        a[j], a[pivot] = a[pivot], a[j]
        d = a[j][j]
        a[j] = [x / d for x in a[j]]
        for i in range(n):
            if i != j:
                d = a[i][j]
                a[i] = [x - d*y for x, y in zip(a[i], a[j])]
    return [row[-1] for row in a]


def intercept_weights(points):
    """Influence of each peak on c0 in c0+C/beta+size offsets."""
    if not points:
        raise ValueError('nonempty peak table required')
    sizes = sorted({rational(p['size']) for p in points})
    slope = len({rational(p['beta_G']) for p in points}) > 1
    design, errors = [], []
    for p in points:
        beta, error = rational(p['beta_G']), rational(p['err'])
        if beta <= 0 or error <= 0:
            raise ValueError('positive beta and marginal error required')
        design.append([F(1)] + ([1/beta] if slope else []) +
                      [F(rational(p['size']) == size) for size in sizes[:-1]])
        errors.append(error)
    n = len(design[0])
    normal = [[sum((row[i]*row[j]/error**2 for row, error in zip(design, errors)), F(0))
               for j in range(n)] for i in range(n)]
    z = solve(normal, [1] + [0]*(n-1))
    weights = [sum((x*y for x, y in zip(row, z)), F(0))/error**2
               for row, error in zip(design, errors)]
    assert all(sum((w*row[j] for w, row in zip(weights, design)), F(0)) == (j == 0)
               for j in range(n))
    return weights


def variance_range(coefficients, errors, groups=None):
    """Sharp PSD variance extrema; distinct groups have zero cross covariance.

    A covariance matrix is a Gram matrix. Within each group, lengths
    |a_i| sigma_i form vectors with freely varying angles. Their sum has
    norm between max(0,2*max(length)-sum(length)) and sum(length).
    """
    if len(coefficients) != len(errors) or not coefficients:
        raise ValueError('equal nonempty coefficient and error vectors required')
    a, sig = list(map(rational, coefficients)), list(map(rational, errors))
    if any(x < 0 for x in sig):
        raise ValueError('nonnegative marginal errors required')
    if groups is None:
        groups = [0]*len(a)
    if len(groups) != len(a):
        raise ValueError('one covariance group per variable required')
    blocks = {}
    for x, s, group in zip(a, sig, groups):
        blocks.setdefault(group, []).append(abs(x)*s)
    lo = hi = F(0)
    for lengths in blocks.values():
        total = sum(lengths, F(0))
        lo += max(F(0), 2*max(lengths)-total)**2
        hi += total**2
    return lo, hi


def encode_range(bounds):
    lo, hi = bounds
    def sqrt_box(x):
        scale = 2**80
        k = math.isqrt(x.numerator*scale**2 // x.denominator)
        lower = F(k, scale)
        upper = lower if lower*lower == x else F(k+1, scale)
        return lower, upper
    low, high = sqrt_box(lo)[0], sqrt_box(hi)[1]
    return {'variance_bounds_exact': [str(lo), str(hi)],
            'sd_bounds_exact_GeV': [str(low), str(high)],
            'sd_range_GeV': [float(low), float(high)]}


def slope_at(function, t):
    """Five-point numerical derivative with a separately reported step check."""
    def estimate(h):
        return (function(t-2*h)-8*function(t-h)+8*function(t+h)-function(t+2*h))/(12*h)
    coarse, fine = estimate(.02), estimate(.01)
    if not math.isfinite(fine) or fine <= 0 or abs(coarse-fine) > 1e-8*abs(fine):
        raise ValueError('positive stable local temperature derivative required')
    return fine, abs(coarse-fine)


def input_custody():
    manifest_path = HERE/'inputs/input_manifest.json'
    manifest = json.loads(manifest_path.read_text())
    hashes = {'inputs/input_manifest.json': hashlib.sha256(manifest_path.read_bytes()).hexdigest()}
    numerical = []
    for row in manifest['files']:
        relative = Path(row['path'])
        if relative.is_absolute() or '..' in relative.parts:
            raise ValueError('escaping comparison input')
        path = HERE/'inputs'/relative
        raw = path.read_bytes()
        if len(raw) != row['bytes'] or hashlib.sha256(raw).hexdigest() != row['sha256']:
            raise ValueError('comparison input differs: '+str(relative))
        if not str(relative).endswith('.raw'):
            numerical.append(str(relative))
        hashes['inputs/'+str(relative)] = row['sha256']
    # These producers return only component means and errors; the paired
    # deletion vectors are local temporaries, not an exported covariance.
    oph = (HERE/'inputs/EW/code/ew_oph.py.raw').read_text()
    wham = (HERE/'inputs/EW/code/ew_wham.py.raw').read_text()
    if 'full, err = rs.jackknife(fall)' not in oph or 'return full, err' not in wham:
        raise ValueError('retained jackknife export contract changed')
    if any(Path(name).suffix != '.json' for name in numerical):
        raise ValueError('inspect additional numerical containers before using the missing-replica premise')
    forbidden = {'replicas', 'jackknife_replicas', 'paired_peak_replicas', 'covariance', 'samples', 'S2', 'S4'}
    def inspect(value):
        if isinstance(value, dict):
            if forbidden.intersection(value):
                raise ValueError('inspect retained joint data before using covariance bounds')
            for item in value.values():
                inspect(item)
        elif isinstance(value, list):
            for item in value:
                inspect(item)
    for name in numerical:
        inspect(json.loads((HERE/'inputs'/name).read_text()))
    return {'retained_numerical_files': numerical, 'paired_peak_replicas_retained': False,
            'raw_simulation_series_retained': False,
            'source_pairing': 'y peak and every T peak evaluated on the same deleted blocks within each run set',
            'export': 'full vector and marginal errors only; paired block vectors discarded',
            'separate_set_independence': 'declared statistical sampling assumption, not inferred from marginal errors'}, hashes


def analyse(production, rows, fit_y, finite_size):
    points = fit_y['points']
    sets = production['chain_point']
    keys = []
    for point in points:
        matches = [key for key, value in sets.items()
                   if value['bG'] == point['beta_G'] and value['L'] == point['n_side']
                   and abs(value['yscan']['x']-.296) < 1e-9]
        if len(matches) != 1:
            raise ValueError('each fitted peak needs one actual run set')
        keys.append(matches[0])
    reference = next((key for key in keys if sets[key]['bG'] == 8 and sets[key]['L'] == 24), None)
    if reference is None:
        raise ValueError('shared reference ensemble absent')
    weights = intercept_weights(points)
    intercept = sum((w*rational(p['value']) for w, p in zip(weights, points)), F(0))
    if abs(float(intercept)-fit_y['continuum']) > 1e-11:
        raise ValueError('influence weights differ from the owning fit')
    sigmas = [rational(p['err']) for p in points]
    custody, hashes = input_custody()
    out = {}
    for name, row in rows.items():
        if 'T_c' not in row:
            continue
        function = lambda t: rd.reduce(t, **row['inputs'], linear=True)['y']
        at_ref, check_ref = slope_at(function, row['T_max_lattice'])
        at_c, check_c = slope_at(function, row['T_c'])
        numerator = [w-F(key == reference) for w, key in zip(weights, keys)]
        coeff = [x/rational(at_c) for x in numerator] + [rational(at_ref)/rational(at_c)]
        errors = sigmas + [rational(row['T_max_lattice_err'])]
        groups = keys + [row['set_used']]
        block = variance_range(coeff, errors, groups)
        all_psd = variance_range(coeff, errors)
        independent = sum((a*a*s*s for a, s in zip(coeff, errors)), F(0))
        out[name] = {'T_c_GeV': row['T_c'], 'reference_set': reference,
                     'temperature_set': row['set_used'],
                     'variables': ['y_max:'+key for key in keys]+['T_max:'+row['set_used']],
                     'influence_coefficients_exact': list(map(str, coeff)),
                     'marginal_errors_exact': list(map(str, errors)),
                     'independent_run_sets': encode_range(block),
                     'unrestricted_cross_set_covariance': encode_range(all_psd),
                     'all_variables_independent_control': encode_range((independent, independent)),
                     'uncorrected_independence_control_GeV': row['T_c_err_lattice'],
                     'derivative_dy_dT': {'at_T_lat': at_ref, 'at_T_c': at_c,
                         'step_halving_difference': [check_ref, check_c],
                         'rigorous_derivative_remainder_certified': False},
                     'reduction_allowance_GeV': row['truncation_error'],
                     'nonlinear_sampling_distribution_identified': False,
                     'confidence_coverage_certified': False}
    # If a*T_lat fluctuation equals the reference-y fluctuation, they
    # cancel exactly before propagation: only the fitted intercept varies.
    cancellation = [w-F(key == reference)+F(key == reference) for w, key in zip(weights, keys)]
    checks = {'intercept_weights_preserve_constant': sum(weights) == 1,
              'shared_reference_cancels_exactly_under_exact_map_control': cancellation == weights,
              'reference_correlation_not_imputed': custody['paired_peak_replicas_retained'] is False,
              'all_ranges_ordered_and_contained': all(
                  F(r['unrestricted_cross_set_covariance']['variance_bounds_exact'][0]) <=
                  F(r['independent_run_sets']['variance_bounds_exact'][0]) <=
                  F(r['independent_run_sets']['variance_bounds_exact'][1]) <=
                  F(r['unrestricted_cross_set_covariance']['variance_bounds_exact'][1]) for r in out.values())}
    for name in ('production.json', 'comparison_budget.py', 'certificate.py', 'reduction.py',
                 'finite_size.py', 'test_comparison_budget.py', 'PROOF.md', 'README.md'):
        hashes[name] = hashlib.sha256((HERE/name).read_bytes()).hexdigest()
    return {'schema': 'oph.electroweak_comparison_budget.v1',
            'scope': 'Sharp PSD delta-method variance ranges conditional on supplied marginal estimates and serialized numerical influence coefficients; no nonlinear or coverage bound',
            'data_custody': custody, 'input_sha256': hashes,
            'fit_intercept': {'value_exact': str(intercept), 'weights_exact': list(map(str, weights)),
                'set_order': keys, 'fit_error_residual_scale': fit_y['error_scale'],
                'residual_scale_is_a_covariance_measurement': False},
            'rows': out, 'finite_design_sensitivity_GeV':
                finite_size['matched_size_cutoff']['T_max']['sensitivity_envelope'],
            'unbounded_effects': ['omitted cutoff powers', 'infinite-volume limit',
                'perturbative reduction remainder', 'omitted simulated hypercharge field',
                'source-to-physical coupling and clock attachment'],
            'checks': checks, 'all_checks_pass': all(checks.values())}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    from .certificate import certificate
    from .__main__ import same
    value = certificate()['comparison_budget']
    path = HERE/'comparison_budget.json'
    if not value['all_checks_pass']:
        raise SystemExit('comparison budget checks failed')
    if args.write:
        path.write_text(json.dumps(value, indent=2, sort_keys=True)+'\n')
    elif not same(json.loads(path.read_text()), value):
        raise SystemExit('comparison budget differs from recomputation')
    print('PASS: shared-reference covariance budget and source bindings')


if __name__ == '__main__':
    main()
