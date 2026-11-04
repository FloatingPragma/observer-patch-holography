"""Source-bound replay of the finite moving scattering packet observation."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from tools import strict_json
from tools.build_completion_program import integration_source_paths
from ..native_spin_observer.precision import generator_norms
from . import bounds, native


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results.json'


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def serialize(value):
    if type(value) is F:
        return str(value)
    if type(value) is dict:
        return {key: serialize(item) for key, item in value.items()}
    if type(value) in (list, tuple):
        return [serialize(item) for item in value]
    return value


def bindings():
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': 'model/scattering_packet_readout/PROOF.md',
        'source_receipts': ['model/moving_contact_scattering/results.json',
                            'model/native_spin_observer/results.json'],
        'independent_checks': ['model/scattering_packet_readout/test_contour.py',
                               'model/scattering_packet_readout/test_cook.py']})
    paths.add('model/scattering_packet_readout/README.md')
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


_IMPORTED_BINDINGS = bindings()


def control_certificate():
    rows = []
    for interacting in (True, False):
        result = native.execute(interacting)
        report = result['report']
        exact = native.literal_certificate(result['word'], interacting)
        probability = report['pointer_probabilities'][0]
        window = (F(739, 1000), F(741, 1000)) if interacting else (F(499999, 10**6), F(500001, 10**6))
        if (not float(window[0]) < probability < float(window[1]) or
                report['complete_vector_error'] >= 1e-3 or
                report['reference_gram_error'] >= 1e-8 or
                report['failure_probability'] >= 1e-20):
            raise ArithmeticError('executed source control exceeded its numerical admission')
        rows.append({**exact, 'interacting': interacting,
                     'joint_quantum_rows': 32, 'complete_pointer_slots': 6,
                     'first_export_records': len(result['first_failure_events']),
                     'first_export_outcome_slots': sum(len(x) for x in result['first_failure_events']),
                     'all_failure_arrays_regenerated': True,
                     'diagnostic_Y_plus_window': list(window),
                     'diagnostic_complete_vector_error_ceiling': F(1, 1000),
                     'diagnostic_reference_gram_error_ceiling': F(1, 10**8),
                     'diagnostic_export_failure_mass_ceiling': F(1, 10**20),
                     'diagnostic_ceilings_are_hardware_guarantees': False})
    return rows


def build():
    before = bindings()
    if before != _IMPORTED_BINDINGS:
        raise ValueError('source changed after import; use a fresh process')
    if sum(generator_norms()) != 29:
        raise ArithmeticError('the exact refining native port coefficient bound changed')
    result = serialize({'schema': 'oph.scattering_packet_readout.v1',
                        'source_bindings': before,
                        'packet': bounds.profile_certificate(),
                        'analytic': bounds.analytic_certificate(),
                        'finite_observation': bounds.plan(),
                        'sampling': bounds.sampling_certificate(),
                        'executed_small_native_controls': control_certificate(),
                        'scope': {'particle_number': 2, 'spatial_dimensions': 3,
                                  'finite_complete_scattering_observation_theorem': True,
                                  'all_eight_valleys_and_all_source_outputs_retained': True,
                                  'same_preparation_and_readout_in_free_control': True,
                                  'analytic_scattering_phase_used_as_native_gate': False,
                                  'fine_spacing_native_program_executed': False,
                                  'finite_spacing_error_sent_to_zero_in_fixed_binary64_ports': False,
                                  'passive_background_drift_or_repair_added_to_control_words': False,
                                  'physical_particle_identification': False,
                                  'many_body_ground_vacuum_constructed': False,
                                  'relativistic_interacting_QFT_constructed': False,
                                  'Einstein_dynamics_constructed': False,
                                  'completion_gate_accepted': False}})
    if before != bindings():
        raise ValueError('source changed during replay')
    return result


def validate(value, replay=None):
    actual = build() if replay is None else replay
    if canonical(value) != canonical(actual):
        raise ValueError('packet proof, complete outcomes, source or scope differs from replay')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--write', type=Path, nargs='?', const=OUTPUT)
    modes.add_argument('--check', type=Path, nargs='?', const=OUTPUT)
    args = parser.parse_args()
    result = build()
    if args.check:
        validate(strict_json.load(args.check), result)
    else:
        args.write.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+'\n')
    print('PASS finite moving scattering packet observation and complete native controls')
