"""Cold replay of the complete finite matrix spectral response."""
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
from . import matrix, identification


ROOT=Path(__file__).resolve().parents[2]
OUTPUT=Path(__file__).resolve().parent/'results.json'


def source_bindings():
    paths=integration_source_paths(ROOT,{
        'domain_and_proof':'model/matrix_spectral_response/README.md',
        'source_receipts':['model/retarded_matter_response/results.json'],
        'independent_checks':['model/matrix_spectral_response/test_matrix.py'],
    })
    paths.add('paper/tex_fragments/MATRIX_SPECTRAL_RESPONSE.tex')
    paths.add(identification.TARGET_PATH)
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def prior_certificate():
    prior=strict_json.loads((ROOT/'model/retarded_matter_response/results.json').read_text())
    if prior.get('all_pass') is not True:raise ValueError('prior response did not pass')
    for path,digest in prior['source_bindings'].items():
        if hashlib.sha256((ROOT/path).read_bytes()).hexdigest()!=digest:
            raise ValueError('prior response source changed: '+path)
    old=prior['resolved_local_spectrum']
    if old['resolved_distinct_positive_lines']!=78 or old['unresolved_local_probe_weight']!='0':
        raise ValueError('prior local spectrum is incomplete')
    return {'all_prior_response_source_hashes_verified':True,
            'prior_local_spectrum_lines':78,
            'prior_full_source_dimension':560,
            'new_matrix_recurrence_and_residues_replayed':True}


def validate(value):
    if value.get('schema')!='oph.matrix_spectral_response.v1' or value.get('all_pass') is not True:
        raise ValueError('invalid matrix response receipt')
    result=value['matrix_response']
    if (result['all_eight_columns_annihilated_exactly'] is not True or
            result['moment_matrices_in_span_I_A']!=78 or
            result['observed_cyclic_dimension']!=312 or
            result['orthogonal_unobserved_dimension']!=248 or
            result['band_line_counts']!={'plus':40,'minus':38} or
            len(result['matrix_residue_rows'])!=78):
        raise ValueError('complete matrix source or cyclic dimension changed')
    for key in ('all_560_eigenvalues_or_degeneracies_resolved',
                'gap_bound_applied_to_unobserved_sector',
                'infinite_volume_locality_or_gap_inferred_from_cube',
                'common_Lorentz_metric_or_physical_masses_identified',
                'native_measurement_executed'):
        if result[key] is not False:raise ValueError('matrix response exceeds its proved scope')
    if result['static_exact']['floating_linear_solver_used'] is not False:
        raise ValueError('static exact arithmetic replaced by a numerical solver')
    comparison=value.get('finite_gap_dictionary')
    if not isinstance(comparison,dict):raise ValueError('finite-gap dictionary test missing')
    if (comparison['finite_gap_as_rest_energy_dictionary_rejected'] is not True or
            comparison['eligible_physical_comparisons_added']!=0 or
            comparison['whole_OPH_model_falsified'] is not False):
        raise ValueError('finite-gap dictionary boundary changed')
    for key in ('unit_factor_or_source_parameters_fitted',
                'finite_gaps_proved_to_be_relativistic_rest_masses',
                'unobserved_sectors_composites_or_continuum_excluded'):
        if comparison[key] is not False:raise ValueError('finite-gap test exceeds its scope')


def canonical(value):return json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n'


def _produce():
    pins=source_bindings()
    value={'schema':'oph.matrix_spectral_response.v1','all_pass':True,
           'prior_proof':prior_certificate(),'matrix_response':matrix.certificate(),
           'finite_gap_dictionary':identification.certificate(),
           'source_bindings':pins,
           'validation_boundary':'Complete source recurrence and rational matrix resolvent, all 78 rank-four residues, exact static coefficients and a reducing-space frequency bound. The 248 orthogonal source dimensions, thermodynamic spectra, Lorentz geometry and native measurement remain outside this proof.'}
    if source_bindings()!=pins:raise ValueError('source changed during matrix response replay')
    validate(value)
    return value


def report():
    pins=source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-matrix-spectral-') as temp:
        path=Path(temp)/'worker.json';environment=dict(os.environ)
        environment['OPENBLAS_NUM_THREADS']=environment['OMP_NUM_THREADS']='1'
        worker=subprocess.run([sys.executable,'-W','error','-X',
                               'pycache_prefix='+str(Path(temp)/'pycache'),'-m',
                               'model.matrix_spectral_response.packet','--worker',str(path)],
                              cwd=ROOT,env=environment,capture_output=True,text=True)
        if worker.returncode:raise RuntimeError('fresh matrix worker failed: '+worker.stderr)
        result=strict_json.loads(path.read_text())
    if source_bindings()!=pins or result['source_bindings']!=pins:
        raise ValueError('source changed across fresh matrix response replay')
    validate(result)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    modes=parser.add_mutually_exclusive_group()
    modes.add_argument('--write',type=Path,nargs='?',const=OUTPUT)
    modes.add_argument('--check',type=Path,nargs='?',const=OUTPUT)
    modes.add_argument('--worker',type=Path,help=argparse.SUPPRESS)
    args=parser.parse_args()
    value=_produce() if args.worker else report()
    if args.check:
        expected=strict_json.loads(args.check.read_text());validate(expected)
        if canonical(expected)!=canonical(value):raise ValueError('matrix response receipt changed')
    elif args.write or args.worker:(args.write or args.worker).write_text(canonical(value))
    if not args.worker:
        print(json.dumps({'all_pass':True,'observed_cyclic_dimension':312,
                          'rank_four_spectral_lines':78,'source_bindings':len(value['source_bindings'])}))


if __name__=='__main__':main()
