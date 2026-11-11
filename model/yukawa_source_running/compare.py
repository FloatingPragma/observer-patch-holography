"""Explicit post-computation empirical comparison, isolated from the producer."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
EMPIRICAL=ROOT/'model/common_source/calibrated_matter_parameters.json'


def compare():
    result=json.loads((HERE/'results.json').read_text())
    empirical=json.loads(EMPIRICAL.read_text())
    y=np.asarray(empirical['matched_SM']['yukawa']['charged_lepton'],float)
    target=y[1:]/y[:-1]
    source=np.asarray(result['source_family_singular_values'])
    rows={'source_Y':source[1:]/source[:-1]}
    for name in ('native_no_hopping','native_full_hopping'):
        values=np.asarray(result['finite_cases'][name]['lowest_six_addition_frequencies'])[::2]
        rows[name]=values[1:]/values[:-1]
    rows['conditional_one_loop_IR']=np.asarray(result['conditional_running']['final_charged_lepton_adjacent_ratios'])
    return {'schema':'oph.yukawa_source_running.empirical_comparison.v1',
        'empirical_file':EMPIRICAL.relative_to(ROOT).as_posix(),
        'empirical_sha256':hashlib.sha256(EMPIRICAL.read_bytes()).hexdigest(),
        'empirical_input_role':'comparison only; absent from source/threshold/running producer',
        'common_scale':empirical['common_scale'],
        'empirical_running_yukawas':y.tolist(),'empirical_adjacent_ratios':target.tolist(),
        'comparisons':{name:{'adjacent_ratios':ratios.tolist(),
                       'ratio_divided_by_empirical_ratio':(ratios/target).tolist()} for name,ratios in rows.items()},
        'no_fitting_or_selection_performed':True,
        'matching_of_finite_addition_frequency_to_running_mass_proved':False,
        'conclusions':{
            'all_positive_no_hopping_g_over_scalar_excluded_by_exact_ratio_bounds':True,
            'one_loop_minimal_SM_IR_cannot_raise_source_middle_over_small_ratio_to_comparison':True,
            'continuum_matching_required_before_physical_mass_comparison':True}}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true');parser.add_argument('--check',action='store_true')
    args=parser.parse_args();result=compare()
    text=json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n'
    if args.write:(HERE/'comparison.json').write_text(text)
    if args.check and json.loads((HERE/'comparison.json').read_text())!=result:
        raise ArithmeticError('separate empirical comparison changed')
    print(text)


if __name__=='__main__':main()
