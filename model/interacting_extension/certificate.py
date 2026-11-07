"""Replay complete interacting detector records and their finite-shot bound."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import numpy as np

from ._source import ROOT, verify_sources
from .higgs_join.native_inputs import SOURCE_NAMES
from .interacting_records import run
from .recorded_response import spectral_shift_tv_squared

HERE=Path(__file__).resolve().parent
OUTPUT=HERE/'results.json'


def source_bindings():
    pins=verify_sources(ROOT,SOURCE_NAMES)
    paths=list(HERE.glob('*.py'))+[HERE/'README.md',HERE/'PROOF.md',HERE/'source_manifest.json']
    paths+=list((HERE/'higgs_join').glob('*.py'))
    paths += [HERE/'higgs_join'/name for name in
              ('joined_basis.json','joined_exact_operators.json','joined_operator_receipt.json',
               'joined_operators.npz','results.json')]
    paths += [HERE/'dependencies/exact_family_operators.json']
    pins.update({p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
                 for p in paths})
    return dict(sorted(pins.items()))


def report(export=None):
    before=source_bindings()
    result,arrays=run(export)
    p=arrays['probabilities'];w=arrays['matrix_weights']
    assert p.shape==(129,2,6,6,2,2) and w.shape==(1596,6,6)
    assert np.isfinite(p).all() and np.min(p)>=0 and np.max(p)<=1
    assert np.max(abs(p.sum(axis=-1)-1))<1e-12
    assert np.max(abs(w-w.swapaxes(1,2).conj()))<1e-12
    assert np.linalg.eigvalsh(w).min()>-1e-12
    bits=np.unpackbits(arrays['packed_binary_records'],axis=-1,count=4096)
    assert bits.size==12681216
    assert np.array_equal(bits.sum(axis=-1),arrays['plus_counts'])
    windows=((2.57e-8,2.60e-8),(.00100665,.00100666),(.4520146,.4520148))
    first=result['first_six_numerically_grouped_lines'][:3]
    assert all(lo<row['frequency_mean']<hi for row,(lo,hi) in zip(first,windows,strict=True))
    assert all(row['multiplicity']==2 for row in first)
    assert result['full_backaction_normalization_error']<2e-10
    assert result['independent_full_joined_control']['maximum_complex_signal_or_norm_error']<2e-10
    assert result['spectral_sum_rule_error']<2e-10
    assert max(result['charge_offset_audit']['consistent_common_gauge_shift_residuals'])<1e-11
    assert result['records']['observed_maximum_quadrature_error']<.1
    times=tuple(F(i,10) for i in range(129))
    tv2=spectral_shift_tv_squared(times,4096,24,F(3,10**8))
    assert tv2<F(1,2500)**2
    if source_bindings()!=before:
        raise RuntimeError('interacting-record source changed during replay')
    return {
        'schema':'oph.interacting_extension.records.v1',
        'source_sha256':before,
        'scope':{
            'declared_lepton_only_gauge_Higgs_model':True,
            'physical_mass_identification':False,'native_ground_preparation_derived':False,
            'physical_clock_calibrated':False,'chiral_continuum_derived':False,
            'Einstein_limit_derived':False,'empirical_mass_inputs':[]},
        'complete_carrier':{'neutral_dimension':438,'charged_dimension':1596,
                            'creation_probes':6,'signed_response_rows':1596},
        'records':{'probability_axes':['time','sign','left_probe','right_probe','quadrature','outcome'],
                   'probability_shape':list(p.shape),'shots_per_diagonal_setting':4096,
                   'simulated_binary_outcomes':int(bits.size),'all_packed_bits_match_counts':True,
                   'all_spectral_residues_retained':True,'blocked_creation_retained':True,
                   'complete_Kraus_actions_checked':True,'forward_evolution_only':True,
                   'complete_gauge_offset_cancellation_checked':True,
                   'joint_Hoeffding_failure_probability':'1/100'},
        'frequency_diagnostics':{
            'acceptance_windows':[['0.0000000257','0.0000000260'],
                                  ['0.00100665','0.00100666'],['0.4520146','0.4520148']],
            'windows_are_rigorous_eigenvalue_enclosures':False,
            'first_three_spin_multiplicities':[2,2,2]},
        'finite_record_information':{
            'record_only_bound_not_a_limit_on_specification_readout':True,
            'time_schedule':'i/10, i=0,...,128','independent_prepared_trials':True,
            'hypotheses':'H1 changes by -delta*P, [P,H1]=0, P orthogonal; all probes/state/other terms fixed',
            'absolute_delta_upper':'3/100000000','record_TV_squared_upper':str(tv2),
            'record_TV_strict_upper':'1/2500',
            'bound_requires_no_frequency_fit':True,
            'quoted_tiny_frequency_is_interval_certified':False}}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write',nargs='?',const=OUTPUT,type=Path)
    mode.add_argument('--check',nargs='?',const=OUTPUT,type=Path)
    parser.add_argument('--export',type=Path)
    args=parser.parse_args();current=report(args.export)
    if args.write is not None:args.write.write_text(json.dumps(current,indent=2,sort_keys=True)+'\n')
    elif json.loads(args.check.read_text())!=current:raise ValueError('interacting record receipt differs')
    print('Complete interacting records, source bindings and finite-shot bound verified.')


if __name__=='__main__':main()
