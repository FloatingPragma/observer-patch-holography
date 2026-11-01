"""Cold source-bound replay for infinite-ground inverse-response locality."""
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
from . import algebra,bounds,connected,euclidean


ROOT=Path(__file__).resolve().parents[2]
OUTPUT=Path(__file__).resolve().parent/'results.json'


def source_bindings():
    paths=integration_source_paths(ROOT,{
        'domain_and_proof':'model/retarded_spatial_locality/README.md',
        'source_receipts':['model/thermodynamic_ground/results.json',
                           'model/field_symmetry/results.json'],
        'independent_checks':[],
    })
    paths.update(['model/full_fock_field/README.md',
                  'model/full_fock_field/LOCAL_SIMULATION.md',
                  'model/full_fock_field/bounds.py',
                  'model/observer_particles/geometric_matter_join.py',
                  'model/relative_observer_geometry/intrinsic_geometry.py',
                  'paper/tex_fragments/RETARDED_SPATIAL_LOCALITY.tex'])
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def prior_proof():
    for package in ('thermodynamic_ground','field_symmetry'):
        prior=strict_json.loads((ROOT/'model'/package/'results.json').read_text())
        bound=prior.get('source_bindings',prior.get('source_sha256'))
        if not bound:raise ValueError('prior proof lacks source custody')
        for path,digest in bound.items():
            if hashlib.sha256((ROOT/path).read_bytes()).hexdigest()!=digest:
                raise ValueError('prior proof source changed: '+path)
        if package=='thermodynamic_ground':
            if (prior.get('all_pass') is not True or
                    prior['existence_theorem']['stationary_under_original_thermodynamic_dynamics'] is not True or
                    prior['existence_theorem']['translation_invariant_ground_state_exists_by_Folner_average'] is not True or
                    prior['density']['number_density_lower']!='1/4' or
                    prior['density']['number_density_upper']!='1'):
                raise ValueError('prior ground or density domain changed')
        elif (prior.get('schema')!='oph.field_symmetry.v1' or
                prior['symmetry']['frame']['all_integer_open_boxes'] is not True or
                prior['symmetry']['frame']['coordinate_period']!=4 or
                not all(row['spin_eta_algebra_exact'] is True
                        for row in prior['symmetry']['full_carrier_checks'])):
            raise ValueError('prior all-box spin symmetry contract changed')
    return {'all_prior_source_hashes_verified':True,
            'ground':'original-mu0 translation-averaged full-Fock ground subsequence at fixed spacing',
            'density_interval':['1/4','1'],
            'prior_expensive_proofs_reexecuted':False,
            'limiting_state_values_computed':False}


def validate(value):
    if value.get('schema')!='oph.retarded_spatial_locality.v1' or value.get('all_pass') is not True:
        raise ValueError('invalid inverse-locality receipt')
    a=value['source_algebra'];b=value['spatial_bounds'];c=value['connected_field']
    e=value['Euclidean_hierarchy']
    if (a['strictly_nonzero_infinite_ground_response_memory'] is not True or
            a['selected_ground_symmetry']['memory_Gram_lower_factor']!='7/64' or
            a['selected_ground_symmetry']['memory_Gram_upper_factor']!='1/4'):
        raise ValueError('selected-ground memory theorem changed')
    for key in ('eta_zero_uniformity_claimed','limiting_state_or_symbol_values_computed',
                'common_Lorentz_metric_or_Einstein_premise_completed','nontrivial_continuum_scattering_claimed'):
        if b[key] is not False:raise ValueError('complex-frequency domain was promoted')
    if (a['selected_ground_symmetry']['ground_state_functional_faithfulness_inferred'] is not False or
            a['operator_norm_large_frequency_limit_claimed'] is not False):
        raise ValueError('unproved state or limit upgrade')
    for key in ('non_Gaussian_multitime_ground_field_proved',
                'linear_CAR_quadratic_evolution_excluded',
                'quadratic_exclusion_includes_pairing_and_auxiliaries'):
        if c[key] is not True:raise ValueError('connected-field theorem changed')
    for key in ('quadratic_exclusion_requires_Gaussian_reference',
                'nonlinear_observable_encodings_excluded',
                'local_state_faithfulness_assumed',
                'continuum_limit_or_scattering_proved',
                'physical_particle_species_identified'):
        if c[key] is not False:raise ValueError('connected-field scope changed')
    if (c['uniform_bounds']['second_derivative_norm']!='169/64' or
            c['uniform_bounds']['time_radius']!='4/169' or
            c['uniform_bounds']['variance_t_squared_lower']!='1/64' or
            c['local_identity']['formal_Wick_residual_terms']!=0 or
            c['local_identity']['anomalous_pairings_retained'] is not True):
        raise ValueError('connected-field constants or Wick identity changed')
    for key in ('all_ordered_local_Euclidean_npoint_limits_exist',
                'same_selected_ground_subsequence',
                'positive_time_reflection_is_word_Gram'):
        if e[key] is not True:raise ValueError('Euclidean hierarchy theorem changed')
    for key in ('complex_time_boundary_rate_computed',
                'limiting_hierarchy_values_computed',
                'bosonic_permutation_invariance_inferred',
                'Euclidean_O4_covariance_inferred',
                'continuum_OS_axioms_completed'):
        if e[key] is not False:raise ValueError('ordered time hierarchy was promoted')
    if (e['all_principal_minors_checked']!=31 or
            e['odd_word_cross_entry']!='-1/64' or
            e['signed_polynomial_norm_squared']!='209/256'):
        raise ValueError('ordered CAR reflection witness changed')


def canonical(value):return json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n'


def _produce():
    pins=source_bindings()
    value={'schema':'oph.retarded_spatial_locality.v1','all_pass':True,
           'prior_proof':prior_proof(),'source_algebra':algebra.certificate(),
           'connected_field':connected.certificate(),
           'Euclidean_hierarchy':euclidean.certificate(),
           'spatial_bounds':bounds.certificate(),
           'proof_scope':'Exact local CAR/frame identities, connected four-point and ordered-reflection witnesses, and rational analytic constants are replayed. The README and paper prove the infinite-operator, time-semigroup hierarchy, bounded-weight and Fourier steps. Finite tests do not compute an infinite ground or certify a real-axis or spatial-continuum limit.',
           'source_bindings':pins}
    if source_bindings()!=pins:raise ValueError('sources changed during replay')
    validate(value)
    return value


def report():
    pins=source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-inverse-locality-') as temp:
        path=Path(temp)/'worker.json';environment=dict(os.environ)
        environment['OPENBLAS_NUM_THREADS']=environment['OMP_NUM_THREADS']='1'
        worker=subprocess.run([sys.executable,'-W','error','-X',
                               'pycache_prefix='+str(Path(temp)/'pycache'),'-m',
                               'model.retarded_spatial_locality.packet','--worker',str(path)],
                              cwd=ROOT,env=environment,capture_output=True,text=True)
        if worker.returncode:raise RuntimeError('fresh locality worker failed: '+worker.stderr)
        value=strict_json.loads(path.read_text())
    if source_bindings()!=pins or value['source_bindings']!=pins:
        raise ValueError('sources changed across cold replay')
    validate(value)
    return value


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--write',type=Path,nargs='?',const=OUTPUT)
    group.add_argument('--check',type=Path,nargs='?',const=OUTPUT)
    group.add_argument('--worker',type=Path,help=argparse.SUPPRESS)
    args=parser.parse_args();value=_produce() if args.worker else report()
    if args.check:
        expected=strict_json.loads(args.check.read_text());validate(expected)
        if canonical(expected)!=canonical(value):raise ValueError('locality receipt or sources changed')
    elif args.write or args.worker:(args.write or args.worker).write_text(canonical(value))
    if not args.worker:
        print(json.dumps({'all_pass':True,'ground_memory_factor':['7/64','1/4'],
                          'non_Gaussian_ground_field':True,
                          'ordered_Euclidean_time_hierarchy':True,
                          'domain':'fixed-spacing Imz>=eta>0',
                          'source_bindings':len(value['source_bindings'])}))


if __name__=='__main__':main()
