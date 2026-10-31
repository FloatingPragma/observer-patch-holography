"""Append an operational co-driven clock to the retained complete field state."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import numpy as np

from . import clock

ROOT = Path(__file__).resolve().parents[2]
FIELD_REPORT = ROOT/'model/defect_field_bridge/dynamics_results.json'
FIELD_ARRAYS = FIELD_REPORT.with_suffix('.npz')
OUTPUT = Path(__file__).with_name('results.json')
ARRAYS = OUTPUT.with_suffix('.npz')
PROBABILITY_GRID = 10**9


def append_to_field(report, field_arrays):
    """Retain old quantum factors and tensor every branch with its clock.

    Clock preparation is at the start. Since it has no feedback until its
    final reads and commutes with nonclock operations, its source execution
    can be contracted separately without changing the joint instrument.
    """
    from ..defect_field_bridge import dynamics_packet
    dynamics_packet.validate_arrays(report, field_arrays)
    duration = F(str(report['duration']))
    first, second = report['ramsey_indices']
    rows, arrays = [], {}
    inherited = F(str(report['ideal_contracted_source_factors']['normalized_input_probability_TV_error_below']))
    for law, branches in report['laws'].items():
        for branch in branches:
            currents = () if branch['abort'] else (branch['current'],)
            run = clock.history(currents, duration)
            tag = branch['array_prefix']
            values = np.stack([run['outputs'][x, y] for x in (0, 1) for y in (0, 1)]).reshape(2, 2, 4, 4)
            arrays[tag+'_clock_work_and_purifiers'] = values
            old = field_arrays[tag+'_all_pointer_probabilities']
            old_mass = float(old.sum())
            marginal = {key: old_mass*value for key, value in run['probabilities'].items()}
            if abs(sum(marginal.values())-old_mass) > float(clock.SERVICE_ERROR)*max(1.,old_mass):
                raise ArithmeticError('adding the clock changed an earlier record marginal')
            result = {'law': law, 'native_record': {'coin': branch['coin'], 'raw': branch['raw'],
                                                  'current': branch['current'], 'abort': branch['abort']},
                      'array_prefix': tag, 'old_pointer_shape': list(old.shape),
                      'new_clock_pointer_shape': [2,2],
                      'prior_probability_mass': old_mass,
                      'joint_clock_marginals': {str(x)+str(y): marginal[x,y] for x,y in marginal},
                      'source_clock_resources': run['source_cost'],
                      'clock_output_work_dimension': 4, 'untouched_clock_purifier_columns': 4,
                      'clock_preparations': 1, 'clock_resets': 0,
                      'fresh_clock_hands': 2, 'fresh_clock_purifiers': 2,
                      'fresh_clock_pointer_bits': 2,
                      'full_clock_instrument_input_columns': 4,
                      'clock_source_instrument_error': run['complete_instrument_frobenius_error'],
                      'expected_effective_time_for_comparison_only': run['effective_time_for_comparison_only']}
            weight = F(str(branch['ideal_native_weight']))
            if weight > 0:
                # This is conditional-frequency inference, not a normalized
                # replacement of any stored quantum output or probability.
                field_error = F(str(branch['joint_probability_error_below']))+inherited
                conditional_error = field_error/weight
                clock_error = conditional_error+clock.SERVICE_ERROR*(1+conditional_error)
                # A declared rational measurement grid makes the exact inverse
                # portable; its worst binary-marginal rounding error is added.
                observed_clock = {key: F(round(value/float(weight)*PROBABILITY_GRID),PROBABILITY_GRID)
                                  for key,value in marginal.items()}
                slot = branch['internal_Ramsey_slot']
                selected = float(old[0,:,first,:].sum() if slot == 0 else old[0,:,:,first].sum())
                observed_field = F(round(selected/float(weight)*PROBABILITY_GRID),PROBABILITY_GRID)
                decoded = clock.infer_energy_gap(observed_field, observed_clock,
                                                field_error=conditional_error+F(1,2*PROBABILITY_GRID),
                                                clock_error=clock_error+F(1,PROBABILITY_GRID))
                decoded['record_probability_grid'] = str(F(1,PROBABILITY_GRID))
                decoded['rounded_clock_probabilities'] = {str(x)+str(y):str(value) for (x,y),value in observed_clock.items()}
                decoded['rounded_field_Ramsey_probability'] = str(observed_field)
                lo, hi = map(F,decoded['energy_gap_GeV_interval'])
                true_gap = F(str(report['specified_energies'][second]))-F(str(report['specified_energies'][first]))
                if not lo <= true_gap <= hi:
                    raise ArithmeticError('record-only clock/field inverse misses its comparison gap')
                result['record_only_inverse'] = decoded
            else:
                result['record_only_inverse'] = {'status':'no conditioning on a zero ideal native branch'}
            rows.append(result)
    return rows, arrays


def build():
    from ..defect_field_bridge import dynamics_packet
    report = json.loads(FIELD_REPORT.read_text())
    field_arrays = dynamics_packet.load_arrays(report, FIELD_ARRAYS)
    rows, arrays = append_to_field(report, field_arrays)
    files = [Path(__file__), Path(clock.__file__), Path(__file__).with_name('README.md'),
             Path(__file__).with_name('test_clock.py'), Path(__file__).with_name('test_packet.py'),
             FIELD_REPORT, FIELD_ARRAYS,
             ROOT/'model/defect_field_bridge/dynamics.py',
             ROOT/'model/defect_field_bridge/dynamics_packet.py',
             ROOT/'model/common_source/mass_response/history.py',
             ROOT/'model/common_source/mass_response/model.py',
             ROOT/'model/native_modular_reference/local_preparation.py']
    packet = {'schema':'oph.defect_field_clock.v1',
              'source_sha256': {p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
              'calibration': clock.calibration(), 'SI_unit_map': clock.si_calibration(),
              'ensemble_sampling_example': clock.sampling_envelope(10**7),
              'joint_state': 'For each retained native record: sum(source_factor tensor internal_factor) tensor clock_work_and_purifiers; all factors are kept without normalization.',
              'parent_quantum_arrays': FIELD_ARRAYS.relative_to(ROOT).as_posix(),
              'native_records':len(rows),'terminal_clock_records':4*len(rows),
              'rows':rows,
              'same_action_convention':'Every field term and H_clock share lambda_j=1+j/4. The clock is prepared once before the field history; only the final two hand readings are made.',
              'noncommuting_piecewise_field_actions_allowed':True,
              'time_from_native_pulse_count_or_host_runtime':False,
              'species_dependent_energy_scale':False,
              'absolute_energy_zero_or_vacuum_offset_observed':False,
              'bounds':'The retained result inverts reconstructed ensemble probabilities, not a single-shot time reading. Exact rational inversion is conditional on the stated numerical source-service probability admission and phase windows. The separate finite-shot envelope can be added for independent repeated preparations; no shots or physical hardware observations were collected.',
              'physical_clock_fit_established':False}
    return packet, arrays


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true');parser.add_argument('--check',action='store_true')
    args=parser.parse_args();packet,arrays=build()
    if args.write:
        OUTPUT.write_text(json.dumps(packet,indent=2)+'\n');np.savez_compressed(ARRAYS,**arrays)
    elif args.check:
        from ..defect_field_benchmark.benchmark import compare_receipt
        compare_receipt(packet,json.loads(OUTPUT.read_text()))
        with np.load(ARRAYS,allow_pickle=False) as data:
            if set(data.files)!=set(arrays):raise ValueError('clock output array keys differ')
            for name,value in arrays.items():
                if not np.allclose(value,data[name],atol=2e-11,rtol=2e-11):
                    raise ValueError('retained complete clock state differs: '+name)
    else:
        print(json.dumps(packet,indent=2))
    print('Complete co-driven field clock and record-only energy inverse verified.')


if __name__=='__main__':main()
