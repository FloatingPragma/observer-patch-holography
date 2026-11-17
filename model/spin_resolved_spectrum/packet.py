"""Cold replay of the complete ground-family finite spectrum."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from tools import strict_json
from tools.build_completion_program import integration_source_paths
from . import spin,spectrum


ROOT=Path(__file__).resolve().parents[2]
OUTPUT=Path(__file__).resolve().parent/'results.json'


def source_bindings():
    paths=integration_source_paths(ROOT,{
        'domain_and_proof':'model/spin_resolved_spectrum/README.md',
        'source_receipts':[spectrum.PRIOR,spectrum.SCALAR],
        'independent_checks':['model/spin_resolved_spectrum/test_spin.py',
                              'model/spin_resolved_spectrum/test_spectrum.py',
                              'model/spin_resolved_spectrum/test_packet.py'],
    })
    paths.add('paper/tex_fragments/SPIN_RESOLVED_SPECTRUM.tex')
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def validate(value):
    if value.get('schema')!='oph.spin_resolved_spectrum.v1' or value.get('all_pass') is not True:
        raise ValueError('invalid spin-resolved spectrum receipt')
    s=value['source_spin'];r=value['resolved_spectrum']
    if (s['ground_family']['creator_columns_checked']!=80 or
            s['ground_family']['annihilator_columns_checked']!=80 or
            s['direct_twirled_moments']['all_16_by_16_spatial_and_spin_entries_checked'] is not True or
            r['reachable_N5_dimension']!=1256 or r['dark_N5_dimension']!=3112 or
            r['mixed_retarded_distinct_signed_poles']!=79 or
            r['mixed_retarded_minimal_cyclic_dimension']!=632 or
            len(r['all_78_positive_addition_lines'])!=78 or
            r['addition_line_counts_by_total_hidden_spin']!={'5/2':1,'3/2':77}):
        raise ValueError('complete source ground family or spectral scope changed')
    for key in ('dark_energy_eigenvalues_or_full_energy_degeneracies_computed',
                'arbitrary_ground_selection_adds_frequencies',
                'native_full_quintet_preparation_or_79_line_resolution_executed',
                'total_hidden_spin_identified_with_particle_Lorentz_spin',
                'physical_particle_species_or_relativistic_rest_masses_identified'):
        if r[key] is not False:raise ValueError('finite internal-spin result was promoted')
    c=r['finite_gap_dictionary']
    if (c['excluded_for_every_ground_choice_and_linear_single_CAR_probe'] is not True or
            c['common_positive_unit_factor_only'] is not True or
            c['eligible_physical_comparisons_added']!=0 or
            c['unseen_composites_other_particle_numbers_or_continuum_excluded'] is not False):
        raise ValueError('ground-family dictionary scope changed')


def canonical(value):return json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n'


def _produce():
    pins=source_bindings()
    value={'schema':'oph.spin_resolved_spectrum.v1','all_pass':True,
           'source_spin':spin.certificate(),'resolved_spectrum':spectrum.certificate(),
           'source_bindings':pins,
           'validation_boundary':'Exact full-source ladders and coherent moments extend the authenticated finite spectrum to every ground vector and linear CAR probe. Internal total-spin labels do not identify particle Lorentz spin, relativistic masses or species; dark-state energies and continuum limits remain uncomputed.'}
    if source_bindings()!=pins:raise ValueError('sources changed during spin-spectrum replay')
    validate(value)
    return value


def report():
    pins=source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-spin-spectrum-') as temp:
        path=Path(temp)/'worker.json';environment=dict(os.environ)
        environment['OPENBLAS_NUM_THREADS']=environment['OMP_NUM_THREADS']='1'
        worker=subprocess.run([sys.executable,'-W','error','-X',
                               'pycache_prefix='+str(Path(temp)/'pycache'),'-m',
                               'model.spin_resolved_spectrum.packet','--worker',str(path)],
                              cwd=ROOT,env=environment,capture_output=True,text=True)
        if worker.returncode:raise RuntimeError('fresh spin-spectrum worker failed: '+worker.stderr)
        result=strict_json.loads(path.read_text())
    if source_bindings()!=pins or result['source_bindings']!=pins:
        raise ValueError('sources changed across spin-spectrum replay')
    validate(result)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    modes=parser.add_mutually_exclusive_group()
    modes.add_argument('--write',type=Path,nargs='?',const=OUTPUT)
    modes.add_argument('--check',type=Path,nargs='?',const=OUTPUT)
    modes.add_argument('--worker',type=Path,help=argparse.SUPPRESS)
    args=parser.parse_args();value=_produce() if args.worker else report()
    if args.check:
        expected=strict_json.loads(args.check.read_text());validate(expected)
        if canonical(expected)!=canonical(value):raise ValueError('spin-spectrum receipt changed')
    elif args.write or args.worker:(args.write or args.worker).write_text(canonical(value))
    if not args.worker:
        print(json.dumps({'all_pass':True,'reachable_N5_dimension':1256,
                          'mixed_reference_signed_poles':79,'source_bindings':len(value['source_bindings'])}))


if __name__=='__main__':main()
