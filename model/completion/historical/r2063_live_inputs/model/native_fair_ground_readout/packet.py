"""Replay the fair-source ground, local phase and complete scored readout."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from tools import strict_json
from model.native_energy_measurement import packet as measurement_owner
from model.native_energy_measurement import measurement, spectral
from . import preparation, reference

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
# Full 66-word replay established at 97633a11f; audited metadata-only renewals
# preserve its scientific payload. A live input check remains mandatory.
MEASUREMENT_RECEIPT_SHA256 = '6c4f66f29e34eac29bac5d7ebe1c857606888eede03113761556c41eb1b0ffeb'
OWNED_FILES = tuple('model/native_fair_ground_readout/' + name for name in (
    '__init__.py', 'preparation.py', 'reference.py', 'packet.py',
    'test_preparation.py', 'test_reference.py', 'test_packet.py', 'PROOF.md', 'README.md'))
ADDITIONAL_FILES = (
    'model/native_exchange_clock/PROOF.md', 'model/federated_observer/SCALING.md',
    'model/native_quantum_reconstruction/PROOF.md', 'model/native_rooted_refinement/PROOF.md',
    'model/native_energy_measurement/results.json')


def encode(value):
    if type(value) is F:
        return str(value)
    if type(value) in (tuple, list):
        return [encode(item) for item in value]
    if type(value) is dict and all(type(key) is str for key in value):
        return {key: encode(item) for key, item in value.items()}
    if value is None or type(value) in (str, int, bool):
        return value
    raise TypeError('exact typed scientific values required')


def canonical(value):
    return json.dumps(encode(value), sort_keys=True, indent=2, allow_nan=False) + '\n'


def source_bindings():
    paths = set(measurement_owner.source_bindings()) | set(OWNED_FILES) | set(ADDITIONAL_FILES)
    return {path: hashlib.sha256((ROOT / path).read_bytes()).hexdigest() for path in sorted(paths)}


def owning_measurement():
    digest = hashlib.sha256(measurement_owner.OUTPUT.read_bytes()).hexdigest()
    if digest != MEASUREMENT_RECEIPT_SHA256:
        raise ArithmeticError('the audited complete measurement receipt SHA differs')
    retained = strict_json.load(measurement_owner.OUTPUT)
    expected = {
        'source_bindings': measurement_owner.source_bindings(),
        'spectral_certificate': spectral.certify(),
        'complete_instrument_budget': measurement.budget(),
        'counted_energy_label_intervals': spectral.energy_labels()}
    for name, value in expected.items():
        if canonical(retained.get(name)) != canonical(value):
            raise ArithmeticError('the complete measurement live sources or consumed bounds differ')
    labels = spectral.energy_labels()
    if any(not 0 <= low <= high <= 30 for low, high in labels):
        raise ArithmeticError('the complete energy score labels escaped their certified range')
    return {
        'audited_full_66_word_receipt_sha256': digest,
        'live_source_map_and_consumed_bounds_recomputed': True,
        'every_retained_word_and_diagnostic_authenticated_by_exact_receipt_sha': True,
        'full_66_word_execution_repeated_here': False}


def produce():
    before = source_bindings()
    owner = owning_measurement()
    value = {
        'schema': 'oph.native_fair_ground_readout.v1',
        'source_bindings': before,
        'inherited_complete_energy_measurement': owner,
        'thermal_vacuum_test': reference.thermal_witness(),
        'executed_fair_ground_preparation': preparation.run(),
        'example_ground_history_bound': reference.history_budget(
            (1, 1), (F(1, 100), F(1, 100)), 630),
        'regional_ground_weights_V240_q80_r2': reference.dicke_weights(240, 80, 2),
        'regional_ground_limit_trace_distance_upper': reference.regional_error(240, 80, 2, F(1, 3)),
        'local_nonseparating_ground_witness_V240_q80': reference.local_ground_witness(240, 80),
        'scope': {
            'common_ground_local_history_bound': True,
            'finite_native_ground_preparation': True,
            'complete_finite_excitation_readout_contrast': True,
            'ground_equals_faithful_transfer_reference': False,
            'growing_coherent_preparation_executed': False,
            'all_level_port_routing_limit_proved': False,
            'passive_repair_equals_controlled_unitary': False,
            'relativistic_3_plus_1_continuum_or_particles_identified': False}}
    if source_bindings() != before:
        raise ArithmeticError('source bytes changed during the ground/readout replay')
    return encode(value)


def check_payload(retained, fresh):
    if canonical(retained) != canonical(fresh):
        raise ArithmeticError('the ground/readout receipt differs from a fresh owning replay')
    return True


def verify(path=OUTPUT):
    fresh = produce()
    check_payload(strict_json.load(Path(path)), fresh)
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.write:
        value = produce()
        OUTPUT.write_text(canonical(value), encoding='utf-8')
        check_payload(strict_json.load(OUTPUT), value)
        print('Fair ground and complete excitation readout receipt written.')
    else:
        verify()
        print('Fair ground and complete excitation readout receipt verified.')


if __name__ == '__main__':
    main()
