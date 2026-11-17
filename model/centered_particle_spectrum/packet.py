"""Cold deterministic replay of the fixed integer spectral witness."""
import argparse,csv,hashlib,io,json,os,subprocess,sys,tempfile
from pathlib import Path
import numpy as np
from tools import strict_json
from tools.build_completion_program import integration_source_paths
from ..centered_source_ground import cube as inherited_ground
from . import certificate,comparison,source

ROOT=Path(__file__).resolve().parents[2]
DIRECTORY=Path(__file__).resolve().parent
OUTPUT=DIRECTORY/'results.json'
WITNESS=DIRECTORY/'witness.npz'
TABLE=DIRECTORY/'lines.csv'
PRIOR=ROOT/'model/centered_source_ground/results.json'


def canonical(value):return json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n'


def load_witness(path=WITNESS):
    with np.load(path,allow_pickle=False) as archive:
        return {name:archive[name] for name in archive.files}


def source_bindings():
    paths=integration_source_paths(ROOT,{
        'domain_and_proof':'model/centered_particle_spectrum/README.md',
        'source_receipts':['model/centered_source_ground/results.json'],
        'independent_checks':[],
    })
    paths.update({'paper/tex_fragments/CENTERED_PARTICLE_SPECTRUM.tex',
                  'model/centered_particle_spectrum/witness.npz',comparison.TARGET})
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


def inherited_inputs():
    prior=strict_json.loads(PRIOR.read_text())
    manifest=prior.get('source_bindings')
    if not isinstance(manifest,dict):raise ValueError('historical source manifest missing')
    mismatches=[name for name,digest in manifest.items()
                if not (ROOT/name).is_file() or hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest]
    return {'mathematical_ground_input':'existing all-number centered ground selection plus explicitly imported Lieb uniqueness theorem',
            'operator_bound_proof':'model/interacting_vacuum/global_ground.py',
            'ground_selection_proof':'model/centered_source_ground/cube.py',
            'current_Lieb_source_hypotheses':inherited_ground.lieb_premises(),
            'current_exact_sector_bound_arithmetic':inherited_ground.sector_bounds(),
            'old_all_number_Gram_execution_repeated_here':False,
            'historical_receipt_path':PRIOR.relative_to(ROOT).as_posix(),
            'historical_receipt_all_pass_flag':prior.get('all_pass'),
            'historical_receipt_source_custody_valid':not mismatches,
            'historical_receipt_source_mismatches':mismatches,
            'historical_custody_is_new_certificate_admission':False}


def validate(value):
    if value.get('schema')!='oph.centered_particle_spectrum.v1' or value.get('all_pass') is not True:
        raise ValueError('invalid centered particle spectrum receipt')
    if value.get('all_pass_scope')!='new exact finite spectral certificate; inherited receipt custody is separately reported':
        raise ValueError('new and inherited custody scope must remain separate')
    s=value['spectrum'];rows=s['rows']
    if (s['positive_visible_lines_per_local_spin_charge_branch']!=228 or len(rows)!=228
            or s['signed_local_retarded_atoms']!=456 or s['unresolved_local_weight']!='0'
            or s['reachable_N7_dimension_all_linear_CAR_probes']!=1824
            or s['reachable_N9_dimension_all_linear_CAR_probes']!=1824
            or s['unobserved_N7_dimension']!=9616 or s['unobserved_N9_dimension']!=9616):
        raise ValueError('complete centered observable spectrum count changed')
    for i,row in enumerate(rows):
        if (row['index']!=i or row['strictly_positive_weight_proved'] is not True
                or row['visible_spatial_multiplicity']!=4
                or len(row['removal_spatial_residue_intervals'])!=8
                or any(len(r)!=8 for r in row['removal_spatial_residue_intervals'])):
            raise ValueError('incomplete positive spectral residue')
    for key in ('global_neutral_gap_certified','all_N7_N9_eigenvalues_resolved',
                'finite_energies_identified_as_relativistic_rest_masses',
                'physical_particle_species_identified','native_observer_measurement_executed',
                'continuum_interaction_or_gravity_proved_by_this_packet'):
        if s[key] is not False:raise ValueError('finite spectral result promoted beyond its domain')
    if value['comparison']['all_three_distinct_frequency_assignments_excluded'] is not True:
        raise ValueError('finite dictionary exclusion missing')
    if (type(value['comparison']['new_eligible_empirical_comparisons']) is not int
            or value['comparison']['new_eligible_empirical_comparisons']!=0):
        raise ValueError('retrospective comparison promoted')


def _produce():
    pins=source_bindings();spectrum=certificate.verify(load_witness())
    result={'schema':'oph.centered_particle_spectrum.v1','all_pass':True,
            'all_pass_scope':'new exact finite spectral certificate; inherited receipt custody is separately reported',
            'source_symmetry':source.symmetry_certificate(),'inherited_inputs':inherited_inputs(),
            'spectrum':spectrum,'comparison':comparison.certificate(spectrum),
            'source_bindings':pins}
    validate(result)
    if pins!=source_bindings():raise ValueError('source changed during exact finite replay')
    return result


def report():
    pins=source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-centered-spectrum-') as temp:
        destination=Path(temp)/'worker.json';environment=dict(os.environ)
        environment['OPENBLAS_NUM_THREADS']=environment['OMP_NUM_THREADS']='1'
        process=subprocess.run([sys.executable,'-W','error','-X',
            'pycache_prefix='+str(Path(temp)/'pycache'),'-m',
            'model.centered_particle_spectrum.packet','--worker',str(destination)],
            cwd=ROOT,env=environment,capture_output=True,text=True)
        if process.returncode:raise RuntimeError('cold exact spectral replay failed:\n'+process.stderr)
        result=strict_json.loads(destination.read_text())
    if pins!=source_bindings() or result['source_bindings']!=pins:
        raise ValueError('source changed across exact cold replay')
    validate(result);return result


def table(value):
    buffer=io.StringIO();writer=csv.writer(buffer,lineterminator='\n')
    writer.writerow(('index','energy_lower','energy_upper','weight_lower','weight_upper','visible_spatial_rank'))
    for row in value['spectrum']['rows']:
        writer.writerow((row['index'],*row['excitation_energy_interval'],
                         *row['local_branch_weight_interval'],row['visible_spatial_multiplicity']))
    return buffer.getvalue()


def main():
    parser=argparse.ArgumentParser(description=__doc__);group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
    group.add_argument('--worker',type=Path);group.add_argument('--propose',action='store_true')
    args=parser.parse_args()
    if args.propose:
        from .proposal import propose
        witness=propose();certificate.verify(witness)
        np.savez_compressed(WITNESS,**witness)
        print('accepted a freshly proposed integer witness; regenerate receipt with --write');return
    if args.worker:args.worker.write_text(canonical(_produce()));return
    value=report()
    if args.write:OUTPUT.write_text(canonical(value));TABLE.write_text(table(value))
    else:
        retained=strict_json.loads(OUTPUT.read_text());validate(retained)
        if canonical(value)!=canonical(retained) or TABLE.read_text()!=table(value):
            raise ValueError('centered spectrum receipt or complete table changed')
    print(json.dumps({'all_pass':True,'visible_positive_lines':228,'signed_atoms':456,
                      'source_bindings':len(value['source_bindings']),
                      'historical_custody_valid':value['inherited_inputs']['historical_receipt_source_custody_valid']}))


if __name__=='__main__':main()
