"""Predict two-order constraint disks from axes, then score frozen interior means."""
import csv
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import numpy as np

from ..measurement_bounds import experiment as detector_source
from ..observer.scalar_response.experiment import canonical, compare_receipt
from ..scalar_observations.experimental_readout import data
from . import bounds, source, randomized, randomized_source

ROOT = Path(__file__).resolve().parents[2]
FOLDER = Path(__file__).resolve().parent
SOURCE_PATHS = tuple(sorted(set((*detector_source.source_bindings(),
    'model/primitive_observer/runtime.py',
    'model/common_source/open_protocol.py',
    'model/calibrated_measurement/bounds.py',
    'model/scalar_observations/experimental_readout/data.py',
    'model/scalar_observations/experimental_readout/provenance.json',
    *(f'model/scalar_observations/experimental_readout/data/{name}' for name in
      ('cf_0ms_real.csv', 'cf_0ms_imag.csv', 'SOURCE_README.md')),
    *(f'model/characteristic_constraints/{name}' for name in
      ('__init__.py', '__main__.py', 'bounds.py', 'source.py', 'randomized.py',
       'randomized_source.py', 'experiment.py', 'protocol.json'))))))


def exact_data():
    loaded = data.load_observations()  # Authenticate bytes before reading exact decimal strings.
    grids = []
    for name in ('cf_0ms_real.csv', 'cf_0ms_imag.csv'):
        with (data.FOLDER/'data'/name).open(newline='') as stream:
            grids.append([[F(x) for x in row] for row in csv.reader(stream)])
    transpose = loaded['orientation']['imaginary_coordinate_array_axis'] == 1
    values = {}
    for j in range(11):
        for i in range(11):
            row, col = (10+i, 10+j) if transpose else (10+j, 10+i)
            values[i, j] = tuple(grid[row][col] for grid in grids)
    axes = {key: value for key, value in values.items() if 0 in key}
    interior = {key: value for key, value in values.items() if 0 not in key}
    return axes, interior, loaded['orientation']


def predictions(axes, extent):
    """This function has no access to the off-axis outcomes."""
    if type(extent) is not int or extent not in (4, 5):
        raise ValueError('retain both declared extent scenarios: 4 and 5')
    if set(axes) != {(i, 0) for i in range(11)} | {(0, j) for j in range(11)}:
        raise ValueError('exactly the 21 axis records required')
    rows = []
    for i in range(1, 11):
        for j in range(1, 11):
            a, b = axes[i, 0], axes[0, j]
            if bounds.norm2(a) >= 1 or bounds.norm2(b) >= 1:
                raise ValueError('axis plug-in means outside this strict-interior protocol; do not clip')
            product = bounds.mul(a, b)
            rows.append({'i': i, 'j': j, 'phase': str(F(extent**2*i*j, 100)),
                'axis_product': list(map(str, product)),
                'radius_squared': str((1-bounds.norm2(a))*(1-bounds.norm2(b)))})
    return rows


def score(axes, interior, prediction):
    rows = []
    for row in prediction:
        i, j = row['i'], row['j']
        a, b, z = axes[i, 0], axes[0, j], interior[i, j]
        for sign in (-1, 1):
            phase = sign*F(row['phase'])
            g = bounds.gram(a, b, z, phase)
            distance = abs(g[1, 2]-g[1, 0]*g[0, 2])
            radius = np.sqrt(float(F(row['radius_squared'])))
            certificate = bounds.witness_floor(a, b, z, phase)
            rows.append({'i': i, 'j': j, 'product_order_sign': sign,
                'minimum_Gram_eigenvalue_diagnostic': float(np.linalg.eigvalsh(g)[0]),
                'disk_excess_diagnostic': float(distance-radius),
                'witness': certificate})
    certified = [row for row in rows if F(row['witness']['complex_entry_correction_lower_bound']) > 0]
    worst = max(rows, key=lambda row: F(row['witness']['complex_entry_correction_lower_bound']))
    return {'complex_interior_points': len(interior), 'order_checks': len(rows),
            'certified_incompatible_plugin_Grams': len(certified),
            'worst_witness': worst, 'checks': rows}


def run():
    protocol = json.loads((FOLDER/'protocol.json').read_text())
    execution = source.execute()
    axes, interior, orientation = exact_data()
    # Freeze both prospective functions before inspecting interior responses here.
    predictions_by_extent = {str(e): predictions(axes, e) for e in protocol['coordinate_extents']}
    sealed = canonical(predictions_by_extent)
    scored = {str(e): score(axes, interior, predictions_by_extent[str(e)]) for e in protocol['coordinate_extents']}
    return {'schema': 'oph.characteristic_constraints.v1',
        'source_sha256': {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in SOURCE_PATHS},
        'source_execution': execution, 'orientation': orientation,
        'randomized_continuing_execution': randomized_source.execute(),
        'future_randomized_designs': randomized.empirical_designs(scored, axes, interior),
        'prediction_sha256': hashlib.sha256(sealed.encode()).hexdigest(),
        'axis_only_predictions': predictions_by_extent, 'empirical_comparison': scored,
        'reported_origin_unchanged': list(map(str, axes[0, 0])),
        'scope': {'new_fitted_parameters': 0, 'Gaussian_state_assumed': False,
            'oscillator_dynamics_assumed': False, 'laboratory_Weyl_instrument_assumed': True,
            'common_normalized_preparation_assumed': True,
            'statistical_rejection_claimed': False, 'OPH_specific_confirmation': False,
            'prospective_preregistration': False, 'raw_trials_available': False,
            'accepted_counts_available': False, 'old_scores_changed': False,
            'randomized_future_common_state_across_trials_assumed': False,
            'randomized_future_independent_trials_assumed': False,
            'randomized_future_independent_setting_choice_assumed': True,
            'randomized_future_physical_services_certified': False}}


def compare(expected, actual):
    compare_receipt(expected, actual)
