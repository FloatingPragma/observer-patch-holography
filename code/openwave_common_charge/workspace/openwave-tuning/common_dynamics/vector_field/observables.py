"""Actual charge-current and field-energy records of the coupled member."""
from pathlib import Path
import hashlib
import json

import numpy as np

import autonomous as source

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def run():
    paths = [Path(__file__), HERE / 'AUTONOMOUS.json', HERE / 'JOINT_RECORDS.json',
             HERE / 'COMPOSITION.json', Path(source.mb.common.__file__)]
    pins = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    records = json.loads((HERE / 'JOINT_RECORDS.json').read_text())
    output = []
    for name, record in records.items():
        charge = int(name.split(':')[0])
        current = np.zeros(288)
        signs = np.zeros(3)
        for row in record['joint_position_field_event_records']:
            x = row['position']
            for event, probability in enumerate(row['event_probabilities']):
                raw = source.mb.common.raw_record(source.mb.common.encode_position(x), event, 'C')
                j = charge * raw['current']
                current[raw['edge_index']] += probability * j
                signs[j + 1] += probability
        assert abs(signs.sum() - 1) < 2e-10
        output.append({'member': name, 'native_current_expectations288': current.tolist(),
                       'native_signed_current_probabilities_minus_zero_plus': signs.tolist(),
                       'mean_electric_energy': record['electric_energy_expectation'],
                       'full_outcomes_preserved': True})
    baseline, coupled = output[0], output[3]
    force_response = max(abs(a - b) for a, b in zip(baseline['native_current_expectations288'], coupled['native_current_expectations288']))
    assert force_response > 1e-6
    autonomous = json.loads((HERE / 'AUTONOMOUS.json').read_text())
    complete = json.loads((HERE / 'COMPOSITION.json').read_text())
    error = complete['three_step_full_record_TV_bound_upper']
    position_lower = autonomous['interaction_position_TV_at_step3'] - 2 * error
    energy_lower = autonomous['interaction_electric_expectation_change'] - 2 * error * 9216
    assert position_lower > 0 and energy_lower > 0
    result = {'schema': 'oph.local12.coupled_record_observables.v1',
              'checks': {'full_joint_source_records_consumed': True,
                         'actual_native_current_decoder': True,
                         'field_changes_source_current_measurement': True,
                         'position_response_exceeds_complete_native_error': True,
                         'electric_response_exceeds_complete_native_error': True,
                         'source_bytes_unchanged': True},
              'records': output, 'maximum_edge_current_response': force_response,
              'certified_position_TV_response_lower': position_lower,
              'certified_electric_expectation_change_lower': energy_lower,
              'current_measurement_scope': 'In the one-particle/one-hole prepared domain, complete native mode occupations decode the source/event address and the frozen inverse-recoding current. Every zero-current and failure flag remains. A force or acceleration dictionary needs the shared physical length/time and inertial mass attachment; this receipt measures current response.',
              'energy_dictionary': 'E_field=sum centered(F_e)^2/2 is the named electric-kick observable. It is not asserted to equal a conserved total quasienergy.',
              'error_dictionary': 'Worst-case total-variation errors include source initialization, full3step native word and all3200 native reads. Both compared programs are bounded by the reported coupled-member bound. Field energy lies between0and9216, including every modular label; assign every abort flag an explicit value in that interval.',
              'source_sha256': pins, 'official_rows_changed': []}
    assert pins == {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (HERE / 'OBSERVABLES.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ('checks', 'maximum_edge_current_response',
                                           'certified_position_TV_response_lower', 'certified_electric_expectation_change_lower')}))


if __name__ == '__main__':
    run()
