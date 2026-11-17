"""One reproducible receipt for the exact and bounded local-clock instances."""
from __future__ import annotations

from fractions import Fraction
import hashlib
from pathlib import Path

from . import core, federation


ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(sorted(set((*core.SOURCE_PATHS, *federation.SOURCE_PATHS,
    'model/observer_time_fields/__main__.py', 'model/observer_time_fields/receipt.py',
    'model/observer_time_fields/federation.py'))))


def json_value(value):
    """Make the producer's exact runtime types agree with strict JSON replay."""
    if type(value) is Fraction:
        return str(value)
    if type(value) in (tuple, list):
        return [json_value(item) for item in value]
    if type(value) is dict:
        if any(type(key) is not str for key in value):
            raise TypeError('receipt keys must be strings')
        return {key: json_value(item) for key, item in value.items()}
    if value is None or type(value) in (str, bool, int):
        return value
    raise TypeError(f'unsupported nonexact receipt value: {type(value).__name__}')


def run():
    elementary = core.summary()
    joined = federation.summary()
    result = {
        'schema': 'oph.observer_time_fields.v1',
        'source_sha256': {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
                          for path in SOURCE_PATHS},
        'elementary_exact_instance': elementary,
        'federation_complete_instance': joined,
        'effective_logarithm': {
            'time_unit': 'incident repair attempt, including zero-current outcomes',
            'input_operator': 'T=(I-A)^(-1/2) B (I-A)^(-1/2)',
            'approximation': 'H_m=sum_(j=1)^m (I-T)^j/j',
            'degree': 512,
            'elementary_operator_remainder_upper': core.logarithm_remainder(512),
            'federation_transfer_lower': Fraction(1, 21),
            'federation_operator_remainder_upper': Fraction(20, 21)**513/Fraction(513, 21),
            'bound_is_full_operator_not_only_gap': True,
            'source_instance_clock_and_dictionary_required': True,
            'Hamiltonian_is_a_native_unitary_time_generator': False,
        },
        'scope': {
            'same_source_local_clock_and_field_words_connected': True,
            'complete_initial_labels_success_and_timeout_records': True,
            'adaptive_internal_native_intervention_included': True,
            'uniform_sector_preparation_supplied': True,
            'complete_authentic_incident_record_access_supplied': True,
            'physical_initialization_or_record_service_error_certified': False,
            'coherent_C_E_completion_selected': False,
            'joint_spatial_refinement_or_relativistic_limit_established': False,
            'new_empirical_observations': False,
        },
    }
    return json_value(result)
