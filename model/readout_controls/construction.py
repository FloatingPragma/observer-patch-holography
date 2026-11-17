"""Join source-executed readout controls to an unchanged public dataset."""
import hashlib
from pathlib import Path

from tools import strict_json
from ..measurement_bounds import experiment as detector_source
from ..compiled_native_interface import construction as compiled
from ..observer.scalar_response.experiment import compare_receipt
from . import readout, prospective


ROOT = Path(__file__).resolve().parents[2]
FOLDER = Path(__file__).resolve().parent
SOURCE_PATHS = tuple(sorted(set((*detector_source.source_bindings(), *compiled.SOURCE_PATHS,
    'model/calibrated_measurement/bounds.py', 'model/anharmonic_observations/data.py',
    'model/anharmonic_observations/provenance.json',
    'model/readout_controls/__init__.py', 'model/readout_controls/__main__.py',
    'model/readout_controls/construction.py', 'model/readout_controls/readout.py',
    'model/readout_controls/prospective.py', 'model/readout_controls/protocol.json',
    *(f'model/anharmonic_observations/data/{name}' for name in readout.data.PUBLISHED_MD5)))))


def run():
    protocol = strict_json.load(FOLDER/'protocol.json')
    source = readout.execute(protocol)
    discrimination = prospective.execute(protocol['prospective']['native_attempts'])
    # No observed response enters either executable prediction.
    empirical = readout.compare_data()
    return {'schema': 'oph.readout_controls.v1',
            'source_sha256': {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in SOURCE_PATHS},
            'null_source_execution': source, 'empirical_null_comparison': empirical,
            'prospective_native_discrimination': discrimination,
            'scope': {'new_OPH_specific_empirical_confirmation': False,
                      'source_selected_laboratory_action': False,
                      'existing_anharmonic_scores_changed': False,
                      'physical_service_error_budget_available': False,
                      'new_fitted_parameters': 0}}


def compare(expected, actual):
    # Exact fractions, source hashes, records and statistical guarantees must
    # agree; only executable numerical diagnostics use the shared tolerance.
    compare_receipt(expected, actual)
