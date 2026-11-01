"""Invert the retained native/colored source pointer into a gap interval.

The inversion uses a stated first Ramsey lobe, common clock and numerical
execution error. It does not use a target mass to choose a phase alias.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from .benchmark import compare_receipt
from .quasiparticles import quasiparticles

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/'model/defect_field_bridge/dynamics_results.json'
ARRAYS = REPORT.with_suffix('.npz')
RECEIPT = Path(__file__).with_name('readout_results.json')


def infer_gap(probability, duration, error):
    """Invert cos²(t gap/2), given 0 <= t gap <= pi and |p-p_true|<=error."""
    if (any(isinstance(x, (bool, str)) or not math.isfinite(float(x))
            for x in (probability, duration, error))
            or not 0 <= probability <= 1 or duration <= 0 or not 0 <= error < 1):
        raise ValueError('finite probability, positive clock and error in [0,1) required')
    lower_p = max(0., probability-error)
    upper_p = min(1., probability+error)
    inverse = lambda value: 2*math.acos(math.sqrt(value))/duration
    return {'central_gap_GeV': inverse(probability),
            'interval_GeV': [inverse(upper_p), inverse(lower_p)],
            'alias_assumption': '0 <= duration*gap <= pi',
            'error_type': 'bounded numerical source-execution error, not sampling uncertainty'}


def infer_charge(probability, error):
    """Invert the separate neutral pointer p0=(1-sin(q))/2 on |q|<=pi/2."""
    if (any(isinstance(x, (bool, str)) or not math.isfinite(float(x)) for x in (probability, error))
            or not 0 <= probability <= 1 or not 0 <= error < 1):
        raise ValueError('finite probability and error in [0,1) required')
    inverse = lambda value: math.asin(1-2*value)
    return {'central_charge': inverse(probability),
            'charge_interval': [inverse(min(1., probability+error)),
                                inverse(max(0., probability-error))],
            'charge_window': '|q|<=pi/2; all seven model charges lie in [-1,1]',
            'probe_law': 'p0=(1-sin(q))/2; specified unit phase and pi/2 quadrature'}


def certificate():
    from ..defect_field_bridge import dynamics_packet

    report = json.loads(REPORT.read_text())
    arrays = dynamics_packet.load_arrays(report, ARRAYS)
    dynamics_packet.validate_arrays(report, arrays)
    inherited_error = float(report['ideal_contracted_source_factors']['normalized_input_probability_TV_error_below'])
    first, second = report['ramsey_indices']
    chart = quasiparticles()
    if ([first, second] != list(chart['ramsey_indices'])
            or report['specified_charges'][first] != report['specified_charges'][second]):
        raise ArithmeticError('source readout changed its equal-charge energy chart')
    programmed_gap = report['specified_energies'][second]-report['specified_energies'][first]
    action_gap = float(chart['energies'][second]-chart['energies'][first])
    rows, zero_rows, total_rows = [], 0, 0
    for law, branches in report['laws'].items():
        for branch in branches:
            total_rows += 1
            probability = arrays[branch['array_prefix']+'_all_pointer_probabilities']
            ideal_weight = branch['ideal_native_weight']
            if ideal_weight == 0:
                zero_rows += 1
                continue
            slot = branch['internal_Ramsey_slot']
            selected = float(probability[0, :, first, :].sum() if slot == 0
                             else probability[0, :, :, first].sum())
            error = (float(branch['joint_probability_error_below'])+inherited_error)/ideal_weight
            observed = selected/ideal_weight
            inferred = infer_gap(observed, branch['effective_duration'], error)
            lower, upper = inferred['interval_GeV']
            if not lower <= programmed_gap <= upper:
                raise ArithmeticError('source-inferred interval misses the programmed gap')
            if not lower-2e-8 <= action_gap <= upper+2e-8:
                raise ArithmeticError('source-inferred interval misses the admitted action gap')
            charge_probability = float(probability[:, 0].sum())/ideal_weight
            charge_inferred = infer_charge(charge_probability, error)
            if not charge_inferred['charge_interval'][0] <= chart['charge6'][first]/6 <= charge_inferred['charge_interval'][1]:
                raise ArithmeticError('source-inferred charge interval misses the actual charge')
            rows.append({'law': law, 'coin': branch['coin'], 'raw': branch['raw'],
                         'native_current': branch['current'], 'native_joint_weight': ideal_weight,
                         'conditioning': 'divide the retained joint pointer probability by exact native weight',
                         'observed_conditional_probability': observed,
                         'conditional_probability_error_bound': error,
                         'effective_duration_GeV_inverse': branch['effective_duration'],
                         'separate_neutral_charge_pointer_zero_probability': charge_probability,
                         **inferred, **charge_inferred})
    if len(rows) != 4 or total_rows != 24:
        raise ArithmeticError('the full C/E source record domain changed')
    paths = [Path(__file__), REPORT, ARRAYS, ROOT/'model/defect_field_bridge/dynamics_controls.json',
             ROOT/'model/defect_field_bridge/dynamics.py', ROOT/'model/defect_field_bridge/dynamics_packet.py',
             ROOT/'model/defect_field_benchmark/quasiparticles.py']
    return {'schema': 'oph.defect_field_benchmark.source_gap_readout.v1',
            'source_sha256': {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
            'all_native_raw_slots_retained_in_source_arrays': total_rows,
            'zero_ideal_native_slots_validated_without_conditioning': zero_rows,
            'inherited_probability_TV_error_added_before_conditioning': inherited_error,
            'nonzero_rows': rows,
            'programmed_gap_GeV_for_comparison_only': programmed_gap,
            'unquantized_action_gap_GeV_for_comparison_only': action_gap,
            'gap_used_to_choose_inverse_alias': False,
            'common_energy_offset_measured': False,
            'common_electric_phase_cancels_on_neutral_pair': True,
            'charge_recovered_from_separate_QND_pointer': True,
            'physical_detector_sampling_performed': False,
            'meaning': 'complete floating source amplitudes recover the calibrated internal gap in the declared clock and first-lobe window'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = certificate()
    if args.write:
        RECEIPT.write_text(json.dumps(result, indent=2)+'\n')
    elif args.check:
        compare_receipt(result, json.loads(RECEIPT.read_text()))
    else:
        print(json.dumps(result, indent=2))
    print('complete source gap readout verified')


if __name__ == '__main__':
    main()
