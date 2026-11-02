"""Fresh-process reproducible neutral-spectrum certificate."""
import argparse,csv,hashlib,io,json,os,subprocess,sys,tempfile
from pathlib import Path
import numpy as np
from tools import strict_json
from tools.build_completion_program import integration_source_paths
from ..centered_source_ground import cube as ground_input
from . import certificate,source
ROOT=Path(__file__).resolve().parents[2]
DIRECTORY=Path(__file__).resolve().parent
OUTPUT=DIRECTORY/'results.json'
WITNESS=DIRECTORY/'witness.npz'
TABLE=DIRECTORY/'lines.csv'
PRIORS=('model/centered_particle_spectrum/results.json','model/centered_source_ground/results.json')


def canonical(value):return json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n'


def bindings():
    paths=integration_source_paths(ROOT,{'domain_and_proof':'model/centered_neutral_spectrum/README.md','source_receipts':list(PRIORS),'independent_checks':[]})
    paths.update({'paper/tex_fragments/CENTERED_NEUTRAL_SPECTRUM.tex','model/centered_neutral_spectrum/witness.npz'})
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def inherited():
    rows=[]
    for path in PRIORS:
        old=strict_json.loads((ROOT/path).read_text());pins=old['source_bindings']
        bad=[p for p,h in pins.items() if not (ROOT/p).is_file() or hashlib.sha256((ROOT/p).read_bytes()).hexdigest()!=h]
        rows.append({'path':path,'historical_source_custody_valid':not bad,'mismatches':bad})
    return {'receipts':rows,'current_Lieb_hypotheses':ground_input.lieb_premises(),
            'current_charged_floor_arithmetic':ground_input.sector_bounds(),
            'old_all_number_experiments_repeated':False,
            'old_integer_ground_vector_directly_verified':True,
            'historical_custody_is_new_admission':False}


def load_witness(path=WITNESS):
    with np.load(path,allow_pickle=False) as z:return {k:z[k] for k in z.files}


def validate(result):
    if result.get('schema')!='oph.centered_neutral_spectrum.v1' or result.get('all_pass') is not True:
        raise ValueError('invalid neutral-spectrum receipt')
    if result.get('all_pass_scope')!='new complete finite neutral-spectrum certificate; inherited custody reported separately':
        raise ValueError('new/inherited admission scope changed')
    r=result['spectrum']
    for key in ('finite_gap_is_relativistic_mass','physical_species_identified','native_detector_executed','continuum_mass_gap_or_gravity_proved'):
        if r[key] is not False:raise ValueError('finite result promoted beyond its domain')
    if r['global']['gap_certifies_full_65536_state_source'] is not True:raise ValueError('full finite gap certificate absent')
    if r['global']['continuum_gap_claimed'] is not False:raise ValueError('finite gap promoted to continuum')
    for row in r['responses']:
        if (type(row['positive_lines']) is not int or row['positive_lines']!=135
                or type(row['exact_dark_lines']) is not int or row['exact_dark_lines']!=1 or row['unresolved_weight']!='0'):
            raise ValueError('incomplete neutral response')
        if len(row['rows'])!=136:raise ValueError('missing local candidates')
        for line in row['rows']:
            if type(line['strictly_positive']) is not bool or type(line['exactly_dark']) is not bool or line['strictly_positive']==line['exactly_dark']:
                raise ValueError('line weight classification invalid')


def produce():
    pins=bindings();spectrum=certificate.verify(load_witness())
    result={'schema':'oph.centered_neutral_spectrum.v1','all_pass':True,
            'all_pass_scope':'new complete finite neutral-spectrum certificate; inherited custody reported separately',
            'spectrum':spectrum,'source_checks':source.source_certificate(),'inherited_inputs':inherited(),'source_bindings':pins}
    validate(result)
    if pins!=bindings():raise ValueError('source changed during neutral replay')
    return result


def report():
    pins=bindings()
    with tempfile.TemporaryDirectory(prefix='oph-neutral-') as temp:
        dest=Path(temp)/'worker.json';env=dict(os.environ)
        env['OPENBLAS_NUM_THREADS']=env['OMP_NUM_THREADS']='1'
        proc=subprocess.run([sys.executable,'-W','error','-X','pycache_prefix='+str(Path(temp)/'pycache'),'-m','model.centered_neutral_spectrum.packet','--worker',str(dest)],cwd=ROOT,env=env,capture_output=True,text=True)
        if proc.returncode:raise RuntimeError('cold neutral replay failed:\n'+proc.stderr)
        result=strict_json.loads(dest.read_text())
    if pins!=bindings() or result['source_bindings']!=pins:raise ValueError('source changed across cold replay')
    validate(result);return result


def table(result):
    out=io.StringIO();writer=csv.writer(out,lineterminator='\n')
    writer.writerow(('probe','index','energy_lower','energy_upper','weight_lower','weight_upper','exactly_dark','visible_spatial_rank','divergence_weight_lower','divergence_weight_upper'))
    for probe in result['spectrum']['responses']:
        for row in probe['rows']:
            writer.writerow((probe['probe'],row['index'],*row['energy_interval'],*row['weight_interval'],row['exactly_dark'],row['visible_spatial_rank'],*row['current_divergence_weight_interval']))
    return out.getvalue()


def main():
    parser=argparse.ArgumentParser(description=__doc__);g=parser.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');g.add_argument('--worker',type=Path);g.add_argument('--propose',action='store_true')
    args=parser.parse_args()
    if args.propose:
        from .proposal import propose
        witness=propose();certificate.verify(witness);np.savez_compressed(WITNESS,**witness)
        print('new integer witness admitted; regenerate receipt');return
    if args.worker:args.worker.write_text(canonical(produce()));return
    result=report()
    if args.write:OUTPUT.write_text(canonical(result));TABLE.write_text(table(result))
    else:
        old=strict_json.loads(OUTPUT.read_text());validate(old)
        if canonical(old)!=canonical(result) or TABLE.read_text()!=table(result):raise ValueError('neutral receipt/table changed')
    print(json.dumps({'all_pass':True,'neutral_gap':result['spectrum']['global']['neutral_gap_interval'],'source_bindings':len(result['source_bindings'])}))


if __name__=='__main__':main()
