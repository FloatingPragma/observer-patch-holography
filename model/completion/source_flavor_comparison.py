"""Exact downstream lepton-ratio comparison, without fitting or producer imports.

The two inputs are retained source intervals and an empirical common-scale
parameter packet. Disjoint ratios reject their direct same-scale identification;
they do not construct a physical scale map or exclude other source laws.
"""
import argparse
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
PATH = Path('model/completion/source_flavor_comparison.json')
SOURCE = 'model/native_flavor_orbit/results.json'
BENCHMARK = 'model/common_source/calibrated_matter_parameters.json'
SCHEMA = 'oph.source_flavor_comparison.v1'
POLICIES = ('always', 'C', 'E', 'CE')
RATIOS = (('lightest_over_heaviest', 0, 2),
          ('middle_over_heaviest', 1, 2),
          ('middle_over_lightest', 1, 0))
DIRECT = (SOURCE, BENCHMARK, 'tools/strict_json.py',
          'model/completion/source_flavor_comparison.py',
          'model/completion/test_source_flavor_comparison.py')
SCOPE = {
    'retained_source_intervals_used_without_producer_import': True,
    'empirical_input_used_only_in_downstream_comparison': True,
    'all_four_declared_duties_retained': True,
    'common_positive_normalization_cancels_from_ratios': True,
    'same_current_canonical_kinetic_basis_changes_preserve_ratios': True,
    'source_parameters_fitted_or_modified': False,
    'same_scale_physical_identification_established': False,
    'renormalization_group_or_threshold_map_computed': False,
    'running_to_pole_mass_map_computed': False,
    'statistical_confidence_or_joint_likelihood_computed': False,
    'physical_mass_prediction_or_unconditional_nature_exclusion': False,
    'all_source_angles_or_feedback_laws_excluded': False,
    'interacting_continuum_or_Einstein_join_established': False,
}
CORRECTION_MODEL = {
    'formula': 'y_i = alpha * sigma_i * (1 + e_i)',
    'domain': 'alpha > 0; |e_i| <= delta < 1; the same alpha for all three families',
    'meaning': 'y_i is the fixed empirical running Yukawa; sigma_i is the source singular value',
    'ratio_bound': '(1-delta)/(1+delta) <= (y_i/y_j)/(sigma_i/sigma_j) <= (1+delta)/(1-delta)',
    'necessary_bound': 'delta >= (R-1)/(R+1), R=max(1,target/upper,lower/target), for each positive ratio interval',
    'combination': 'take the maximum necessary lower bound over the three retained ratios',
    'zero_case': 'a zero source singular value cannot become positive under any finite multiplicative correction',
    'not_asserted': 'the lower bound is not a computed RG correction, fitted correction, or sufficient matching condition',
}
PRIMARY_SOURCES = [
    {'role': 'published common-scale running-parameter and matching conventions',
     'url': 'https://arxiv.org/html/2510.01312v2#S2'},
    {'role': 'SMDR renormalization-group and on-shell matching formulation',
     'url': 'https://arxiv.org/abs/1907.02500'},
    {'role': 'electron experimental mass convention; not an additional numerical input',
     'url': 'https://pdg.lbl.gov/2025/listings/rpp2025-list-electron.pdf'},
    {'role': 'muon experimental mass convention; not an additional numerical input',
     'url': 'https://pdg.lbl.gov/2025/listings/rpp2025-list-muon.pdf'},
    {'role': 'tau experimental mass convention; not an additional numerical input',
     'url': 'https://pdg.lbl.gov/2025/listings/rpp2025-list-tau.pdf'},
]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def rational(value):
    require(type(value) is str, 'exact rational strings are required')
    try:
        return Fraction(value)
    except (ValueError, ZeroDivisionError) as error:
        raise ValueError('finite rational string required') from error


def ratio_interval(numerator, denominator):
    """Outward exact division of nonnegative intervals, or an explicit zero pole."""
    require(len(numerator) == len(denominator) == 2, 'two interval endpoints required')
    a, b = map(rational, numerator)
    c, d = map(rational, denominator)
    require(0 <= a <= b and 0 <= c <= d, 'ordered nonnegative intervals required')
    if c == d == 0:
        return None
    require(c > 0, 'denominator interval must be positive or exactly zero')
    return [str(a/d), str(b/c)]


def necessary_correction(interval, target):
    """Necessary uniform relative-error bound, with an explicit exact-zero case."""
    target = rational(target)
    require(target > 0, 'positive empirical ratio required')
    require(type(interval) is list and len(interval) == 2, 'ratio interval required')
    lo, hi = map(rational, interval)
    require(0 <= lo <= hi, 'ordered nonnegative ratio required')
    if lo == hi == 0:
        return None
    require(lo > 0, 'positive or exactly zero source ratio required')
    stretch = max(Fraction(1), target/hi, lo/target)
    return str((stretch-1)/(stretch+1))


def source_bindings(root=ROOT):
    root = Path(root).resolve()
    return {name: hashlib.sha256((root/name).read_bytes()).hexdigest()
            for name in sorted(DIRECT)}


def _inputs(root):
    source, benchmark = (strict_json.load(root/name) for name in (SOURCE, BENCHMARK))
    require(type(source) is dict and source.get('schema') == 'oph.native_flavor_orbit.v1'
            and source.get('all_pass') is True
            and source.get('canonical_physical_premises_discharged') == [],
            'retained finite source receipt required')
    spectrum = source.get('current_gated_joint_spectrum')
    require(type(spectrum) is dict and spectrum.get('declared_angle_exact') == 'pi/5'
            and type(spectrum.get('declared_reset_port')) is int
            and spectrum['declared_reset_port'] == 2
            and spectrum.get('joint_action_uses_declared_static_mean_feedback_service') is True
            and spectrum.get('joint_action_equals_finite_correlated_memory_channel') is False,
            'declared source protocol or static-mean domain changed')
    rows = spectrum.get('policies')
    require(type(rows) is list and all(type(row) is dict for row in rows)
            and [row.get('policy') for row in rows] == list(POLICIES),
            'all four ordered source policies required')
    require(type(benchmark) is dict
            and benchmark.get('schema') == 'oph.calibrated_matter.parameters.v1'
            and benchmark.get('role') == 'empirical_calibration_sanity_check_not_parameter_derivation',
            'empirical benchmark packet required')
    scale = benchmark.get('common_scale')
    require(type(scale) is dict and scale.get('scheme') == 'MSbar Standard Model'
            and scale.get('mu') == 'M_Z', 'common MSbar M_Z benchmark required')
    try:
        values = benchmark['matched_SM']['yukawa']['charged_lepton']
    except (KeyError, TypeError) as error:
        raise ValueError('three benchmark charged-lepton values required') from error
    require(type(values) is list and len(values) == 3, 'three benchmark charged-lepton values required')
    y = list(map(rational, values))
    require(0 < y[0] < y[1] < y[2], 'ordered positive benchmark values required')
    return spectrum, benchmark, y


def build(root=ROOT):
    """Recompute exact comparison arithmetic from two retained files only."""
    root = Path(root).resolve()
    pins = source_bindings(root)
    spectrum, benchmark, y = _inputs(root)
    comparisons = []
    for policy in spectrum['policies']:
        channels = policy.get('channels')
        require(type(channels) is list and all(type(row) is dict for row in channels)
                and [row.get('channel') for row in channels] == ['up', 'down', 'charged_lepton'],
                'complete source channel order required')
        channel = channels[2]
        boxes = channel.get('dimensionless_singular_values')
        require(type(boxes) is list and len(boxes) == 3
                and all(type(box) is list and len(box) == 2 for box in boxes),
                'all three source singular-value intervals required')
        bounds = [tuple(map(rational, box)) for box in boxes]
        require(all(0 <= a <= b for a, b in bounds)
                and all(bounds[i][1] < bounds[i+1][0] for i in (0, 1)),
                'ordered disjoint source singular-value intervals required')
        rank = 3 if policy['policy'] in ('always', 'C') else 2
        require(type(channel.get('rank')) is int and channel['rank'] == rank
                and (bounds[0][0] > 0 if rank == 3 else bounds[0] == (0, 0)),
                'source rank or exact zero direction changed')
        ratios = {}
        corrections = []
        for name, i, j in RATIOS:
            target = str(y[i]/y[j])
            interval = ratio_interval(boxes[i], boxes[j])
            correction = None if interval is None else necessary_correction(interval, target)
            if correction is not None:
                corrections.append(Fraction(correction))
            ratios[name] = {
                'source_interval_exact': interval,
                'undefined_reason': 'exactly zero source denominator' if interval is None else None,
                'benchmark_exact': target,
                'benchmark_outside_source_interval': None if interval is None else not (
                    rational(interval[0]) <= rational(target) <= rational(interval[1])),
                'necessary_delta_lower_bound_exact': correction,
            }
        zero = rank == 2
        comparisons.append({
            'policy': policy['policy'], 'source_rank': rank,
            'source_singular_value_intervals_exact': deepcopy(boxes),
            'ratios': ratios,
            'direct_common_scale_identification_excluded': any(
                row['benchmark_outside_source_interval'] is True for row in ratios.values()),
            'exact_zero_source_mode_against_positive_benchmark': zero,
            'finite_multiplicative_correction_possible': not zero,
            'necessary_uniform_delta_lower_bound_exact': None if zero else str(max(corrections)),
            'correction_bound_is_sufficient_for_joint_matching': False,
        })
    result = {
        'schema': SCHEMA, 'source_bindings': pins,
        'calculation_passed': True,
        'canonical_physical_premises_discharged': [],
        'physical_predictions_added': 0,
        'source_protocol': {
            'angle_exact': spectrum['declared_angle_exact'], 'reset_port': 2,
            'policies': list(POLICIES),
            'coefficient': 'stationary first-order real-projector feedback, canonical kinetic form and unit Higgs norm',
        },
        'benchmark': {
            'common_scale': deepcopy(benchmark['common_scale']),
            'charged_lepton_yukawas_exact': list(map(str, y)),
            'values_are': 'fixed empirical calibration constants; not exact measurements or a joint likelihood',
        },
        'conditional_identification': {
            'family_order': 'ascending source singular values identified with electron, muon, tau',
            'normalization': 'one common positive scale, including a common Higgs norm or clock conversion',
            'required_transport': 'same physical fields, canonical kinetic normalization, renormalization scheme and scale',
            'not_supplied': 'source-to-MSbar scale, RG and threshold transport, physical species and running-to-pole matching',
        },
        'correction_model': deepcopy(CORRECTION_MODEL),
        'comparisons': comparisons,
        'all_policies_conditionally_disjoint': all(
            row['direct_common_scale_identification_excluded'] for row in comparisons),
        'scope': deepcopy(SCOPE),
        'primary_source_metadata': deepcopy(PRIMARY_SOURCES),
        'external_data_fetched_during_calculation': False,
    }
    require(source_bindings(root) == pins, 'input changed during exact comparison')
    return result


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def validate(packet, root=ROOT):
    require(type(packet) is dict and packet.get('schema') == SCHEMA, 'wrong comparison schema')
    require(canonical(packet) == canonical(build(root)), 'comparison differs from exact retained-input replay')
    return {'policies': list(POLICIES), 'canonical_physical_premises_discharged': [],
            'all_policies_conditionally_disjoint': packet['all_policies_conditionally_disjoint'],
            'physical_predictions_added': 0}


def write(path=None, root=ROOT):
    target = Path(root)/PATH if path is None else Path(path)
    result = build(root)
    target.write_text(canonical(result))
    return result


def check(path=None, root=ROOT):
    target = Path(root)/PATH if path is None else Path(path)
    return validate(strict_json.load(target), root)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--write', nargs='?', const=str(ROOT/PATH))
    mode.add_argument('--check', nargs='?', const=str(ROOT/PATH))
    args = parser.parse_args()
    if args.write:
        result = write(args.write)
        print(json.dumps({'written': args.write, 'policies': list(POLICIES),
                          'all_policies_conditionally_disjoint': result['all_policies_conditionally_disjoint']}))
    else:
        print(json.dumps(check(args.check), indent=2))


if __name__ == '__main__':
    main()
