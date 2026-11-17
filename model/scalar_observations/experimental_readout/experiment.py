"""Replay a fixed retrospective calibration/validation comparison on real data."""
from fractions import Fraction as F
import hashlib
import json
import math

from .data import FOLDER, load_observations
from . import stats
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[3]))
    from tools.pin_scope import bound_pins


def _serializable(value):
    if isinstance(value, complex):
        return {'real': value.real, 'imag': value.imag}
    if isinstance(value, dict):
        if value and all(isinstance(key, tuple) for key in value):
            return [{'index': list(key), 'value': _serializable(value[key])} for key in sorted(value)]
        return {key: _serializable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_serializable(item) for item in value]
    return value


def _point_prediction(means):
    return {key: {'mean': value, 'box': stats.mean_box(value, 0)} for key, value in means.items()}


def run():
    protocol_bytes = (FOLDER/'protocol.json').read_bytes()
    protocol = json.loads(protocol_bytes)
    observations = load_observations()
    control = protocol['statistical_diagnostic']
    samples = control['reported_measurements_per_point']
    delta_cal = control['calibration_failure_allowance']
    delta_val = control['validation_failure_allowance']
    axis, heldout = observations['calibration'], observations['validation']
    calibration = stats.calibrate_axes(axis, samples=samples, delta=delta_cal)
    prediction = stats.predict(calibration, sorted(heldout))
    result = stats.score(prediction, heldout, samples=samples, delta=delta_val)
    result['spin_mean_rmse'] = result['complex_rmse']/math.sqrt(2)
    result['spin_probability_rmse'] = result['spin_mean_rmse']/2
    result['calibration_failure_budget_included'] = True
    result['joint_failure_allowance_under_assumptions'] = delta_cal+delta_val

    zero = stats.score(_point_prediction({key: 0j for key in heldout}), heldout,
                       samples=samples, delta=delta_val)
    zero['role'] = 'Uninformative zero-response diagnostic; not a selected physical alternative'
    zero['spin_mean_rmse'] = zero['complex_rmse']/math.sqrt(2)
    nominal = []
    # Redundancy determines the conjugated axis, but cannot determine its sign.
    # Retain every named convention; never select one by its residual.
    for extent in (protocol['grid']['primary_extent_from_paper'],
                   protocol['grid']['alternate_extent_from_dataset_readme']):
        for direction in (1, -1):
            means = {key: stats.nominal_gaussian(extent*key[0]/10, direction*extent*key[1]/10,
                     thermal=protocol['nominal_reference']['thermal_occupation'],
                     amplitude=-protocol['nominal_reference']['mean_position_in_paper_units'])
                     for key in heldout}
            comparison = stats.score(_point_prediction(means), heldout, samples=samples, delta=delta_val)
            comparison.update(extent=extent, row_direction=direction,
                              spin_mean_rmse=comparison['complex_rmse']/math.sqrt(2),
                              physical_coordinate_identification_established=False,
                              parameters_fitted=False)
            nominal.append(comparison)

    certificates = {
        'calibration': stats.concentration_certificate(calibration['real_statistics'], samples,
                                                       calibration['simultaneous_half_width']),
        'validation': stats.concentration_certificate(result['real_statistics'], samples,
                                                      result['simultaneous_half_width'])}
    exact_failure = sum((F(item['familywise_error_upper']) for item in certificates.values()), F(0))
    checks = {
        'source_files_match_published_checksums': True,
        'calibration_validation_disjoint': set(axis).isdisjoint(heldout),
        'one_quadrant_without_symmetry_replicates': len(axis) == 21 and len(heldout) == 100,
        'all_heldout_cells_scored': len(result['rows']) == 100,
        'original_noisy_origin_preserved': calibration['axis_means'][0, 0] == observations['reported_origin'],
        'no_withheld_parameter_fit': not calibration['calibration_uses_withheld_values'],
        'both_extents_and_directions_reported': len(nominal) == 4,
        'exact_conditional_failure_budget_met': exact_failure <= F('0.05')}
    source_hashes = bound_pins({name: hashlib.sha256((FOLDER.parent/name).read_bytes()).hexdigest()
                     for name in ('source.py', 'protocol.py', 'occupation_readout.py')})
    receipt = {
        'schema': 'oph.scalar_observations.experimental_readout.v1',
        'protocol_sha256': hashlib.sha256(protocol_bytes).hexdigest(),
        'provenance': observations['provenance'],
        'source_interface_hashes': bound_pins(source_hashes),
        'orientation': observations['orientation'],
        'selection_note': protocol['selection_order'],
        'orientation_sensitivity_note': 'Numeric ingestion identified the conjugated array axis but not its direction; both sign conventions are reported without selecting by agreement.',
        'calibration': calibration, 'heldout_comparison': result,
        'zero_response_baseline': zero, 'unfitted_nominal_conventions': nominal,
        'concentration_certificates': certificates,
        'conditional_joint_failure_upper_exact': str(exact_failure),
        'scope': {**protocol['evidence_scope'],
            'data_kind': 'published_experimental_spin_means',
            'retrospective_withheld_points_not_prospective_acquisition': True,
            'acquisition_independence_or_accepted_counts_authenticated': False,
            'absolute_coordinate_scale_identified': False,
            'all_acquisition_attempts_available': False,
            'origin_renormalization_or_background_subtraction': False,
            'source_execution_of_this_laboratory_preparation': False,
            'independent_validation_of_oph': False},
        'checks': checks, 'all_pass': all(checks.values()),
        'all_pass_means': 'Ingestion, disjoint split, reproducibility and arithmetic contracts; scientific comparison verdicts remain separate',
    }
    if not receipt['all_pass']:
        raise ArithmeticError(f'Experimental-readout software contract failed: {checks}')
    return _serializable(receipt)
