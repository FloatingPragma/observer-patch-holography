"""Fresh exact-carrier and complete-response replay for the finite Higgs join."""
from __future__ import annotations

import argparse
from contextlib import redirect_stdout
from fractions import Fraction
import hashlib
import io
import json
from pathlib import Path
import tempfile

import numpy as np
import sympy as sp

from ..dependencies.source_reference import exact_source_packet
from . import audit_exact_join
from .interacting_carrier import HERE,ROOT,build_join,primary_family_matrices,produce_operators
from .run_joined_response import addition_data,block_eigh,operator_checks,schur_checks

OUTPUT=HERE/'results.json'


def source_bindings():
    paths={p for p in HERE.iterdir() if p.suffix in {'.py','.json','.npz','.md'} and p!=OUTPUT}
    paths.update(HERE.parent/n for n in ('__init__.py','_source.py','source_manifest.json',
        'dependencies/__init__.py','dependencies/source_reference.py','dependencies/exact_family_operators.json'))
    paths.add(ROOT/'model/requirements.txt')
    return {p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}


def _window(value,lower,upper):
    if not float(Fraction(lower))<value<float(Fraction(upper)):
        raise ArithmeticError('configured numerical response left its replay window')
    return {'lower':lower,'upper':upper}


def report():
    before=source_bindings()
    source=exact_source_packet(ROOT)
    expressions=json.loads((HERE.parent/'dependencies/exact_family_operators.json').read_text())
    for name in ('Y','D','J'):
        if any(sp.simplify(v)!=0 for v in sp.Matrix(expressions[name])-source[name]):
            raise ArithmeticError('family expression is not its source record value')
    with tempfile.TemporaryDirectory(prefix='oph-higgs-replay-') as directory:
        output=Path(directory)
        with redirect_stdout(io.StringIO()):
            exact=produce_operators(ROOT,output,3)
            audit=audit_exact_join.run(ROOT,HERE,output)
        for name in ('joined_basis.json','joined_exact_operators.json'):
            if (output/name).read_bytes()!=(HERE/name).read_bytes():
                raise ArithmeticError('exact source reconstruction differs: '+name)
        with np.load(output/'joined_operators.npz',allow_pickle=False) as fresh, \
             np.load(HERE/'joined_operators.npz',allow_pickle=False) as saved:
            if set(fresh.files)!=set(saved.files):raise ArithmeticError('operator array inventory')
            for name in fresh.files:
                if fresh[name].shape!=saved[name].shape or not np.allclose(fresh[name],saved[name],rtol=0,atol=1e-14):
                    raise ArithmeticError('normalized exact operator differs: '+name)
    y,sigma=primary_family_matrices()
    joined=build_join(HERE,family0=y,family1=sigma,yukawa=y)
    eigen={q:block_eigh(joined['sectors'][q]) for q in (0,1)}
    checks={q:operator_checks(joined['sectors'][q]) for q in (0,1)}
    arrays={}
    addition=addition_data(joined,eigen,'control',arrays)
    schur=schur_checks(joined['sectors'][1],arrays,'control')
    if eigen[0][3]['ground_numerical_multiplicity']!=1 or eigen[1][3]['ground_numerical_multiplicity']!=2:
        raise ArithmeticError('configured numerical ground multiplicity changed')
    frequencies=arrays['control_addition_frequencies']
    windows=[_window(frequencies[i],lo,hi) for i,lo,hi in (
        (0,'2579/100000000000','2581/100000000000'),
        (2,'1006653/1000000000','1006654/1000000000'),
        (4,'4520147/10000000','4520148/10000000'))]
    self_energy=_window(schur['resolvent_checks'][0]['self_energy_operator_norm'],
                        '1793879/1000000','1793880/1000000')
    if schur['coupling_gram_min_eigenvalue']<=1.8:
        raise ArithmeticError('configured Higgs coupling is not full column rank')
    if checks[1]['Higgs_hopping_commutator_frobenius']<=12.5:
        raise ArithmeticError('configured Higgs and hopping do not exhibit the expected interaction')
    expected={'gauss':48816,'spin':36612,'car':15768,'reference':12204}
    actual={'gauss':audit['all_lifted_24_Gauss_constraints_checked'],
            'spin':audit['all_raw_literal_spin_CAR_Higgs_Ward_checks'],
            'car':audit['all_raw_Q0_addition_CAR_anticommutators_checked'],
            'reference':audit['all_explicit_graded_reference_probe_CAR_sign_checks']}
    if actual!=expected:raise ArithmeticError('complete exact audit coverage differs')
    result={
        'schema':'oph.interacting.higgs_join.replay.v1','source_sha256':before,
        'exact':{'complete_Gauss_dimensions':exact['complete_Gauss_dimensions'],
                 'operator_complement_columns_checked':exact['complete_operator_leakage_columns_checked'],
                 'full_basis_checks':actual,'source_family_expressions_reconstructed':True},
        'configured_response':{'all_final_states_retained':addition['all_Q1_eigenstate_rows_retained'],
            'all_initial_final_pairs_retained':addition['all_initial_final_pairs_retained'],
            'probe_count':len(joined['addition']),
            'first_three_distinct_frequency_replay_windows':windows,
            'self_energy_norm_at_1_plus_half_i_replay_window':self_energy,
            'complete_eigensystem_residual_tolerance':2e-11,
            'complete_spectral_sum_tolerance':1e-11,
            'Schur_resolvent_tolerance':1e-10,
            'replay_windows_are_rigorous_interval_enclosures':False},
        'scope':{'finite_declared_control':True,'source_selected_couplings':False,
                 'source_selected_common_vacuum':False,'scalar_occupancy_cap':1,
                 'quark_pair_vertices_enabled':False,'physical_clock_or_GeV_unit_derived':False,
                 'measured_particle_masses_derived':False,'interacting_chiral_continuum_derived':False}}
    if source_bindings()!=before:raise ArithmeticError('source changed during replay')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write',type=Path,nargs='?',const=OUTPUT)
    mode.add_argument('--check',type=Path,nargs='?',const=OUTPUT)
    args=parser.parse_args();result=report()
    if args.write is not None:
        args.write.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
    elif json.loads(args.check.read_text())!=result:
        raise ArithmeticError('finite Higgs replay receipt differs')
    print('Verified exact finite Higgs carrier and complete configured response.')


if __name__=='__main__':
    main()
