"""Seal the executed row30 candidate and unchanged imported source bytes."""
from pathlib import Path
from fractions import Fraction
import hashlib
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(name):return json.loads((HERE/name).read_text())
def save(name,data):(HERE/name).write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')

def main():
    result=load('RESULT.json');extra=load('ADDITIONAL_CONTROLS.json');rest=load('REST_CONTROL.json')
    assert result['full_criterion_candidate'] and not result['full_criterion_met']
    assert len(result['finite_native_instruments'])==18 and len(result['continuum_spectroscopy'])==72
    assert len(extra['finite_native_budgets'])==6 and len(rest['rows'])==3
    assert extra['missing_export_accounting_rejected'] and extra['all_provider_pins_unchanged']
    assert rest['same_native_preparation_clock_readout'] and rest['provider_pins_unchanged']
    sources=result['source_sha256']
    for rel,h in sources.items():assert sha(ROOT/rel)==h,rel
    imports=load('IMPORTS.json')['files_sha256']
    for rel,h in imports.items():assert sources[rel]==h and sha(ROOT/rel)==h,rel
    summary={
      'candidate_id':result['candidate_id'],'full_criterion_candidate':True,
      'full_criterion_met':False,'awaiting':'Two producer-independent full-criterion reviews',
      'native_controlled_primitive_complete_columns':sum(x['complete_columns'] for x in result['native_primitive_checks']),
      'finite_source_preparation_instruments':18,'complete_native_fine_record_settings':288,
      'native_rest_frequency_record_settings':12,'continuum_cases':72,
      'largest_controlled_native_matrix_residual':max(x['operator_residual'] for x in result['native_primitive_checks']),
      'largest_native_preparation_residual':max(x['native_preparation_error'] for x in result['finite_native_instruments']),
      'formal_family_last_size':8192,
      'formal_family_last_complete_TV_upper':float(Fraction(result['formal_complete_native_family'][-1]['native_vs_continuum_uncapped_exact'])),
      'formal_family_last_pulses_upper':result['formal_complete_native_family'][-1]['whole_experiment_pulses_upper'],
      'formal_family_last_memory_bits':result['formal_complete_native_family'][-1]['total_memory_bits'],
      'whole_large_native_words_executed':False,'source_bytes_unchanged':True,
      'runtime_import_files':len(imports),'producer_files':len(sources),'official_rows_changed':[]}
    save('SUMMARY.json',summary)
    names=['native_spectroscopy.py','run.py','additional_controls.py','rest_control.py','seal.py',
           'README.md','PROOF.md','CONTRACT.json','CRITERION.json','RESULT.json','COMPLETE_RECORDS.json',
           'ADDITIONAL_CONTROLS.json','REST_CONTROL.json','IMPORTS.json','SUMMARY.json',
           'RUN.log','ADDITIONAL_RUN.log','REST_RUN.log']
    names += [str(p.relative_to(HERE)) for p in sorted((HERE/'upstream').iterdir()) if p.is_file()]
    files={name:sha(HERE/name) for name in names}
    save('MANIFEST.json',{'schema':'oph.row30.native_kg_completion.manifest.v1',
         'candidate_id':result['candidate_id'],'criterion_sha256':load('CRITERION.json')['criterion_sha256'],
         'files_sha256':files,'source_sha256':sources,'provider_bytes_unchanged':True,
         'full_criterion_candidate':True,'full_criterion_met':False,
         'two_independent_full_reviews':False,'official_rows_changed':[]})
    save('VALIDATION.json',{'row':30,'candidate_id':result['candidate_id'],
         'criterion_sha256':load('CRITERION.json')['criterion_sha256'],
         'manifest_sha256':sha(HERE/'MANIFEST.json'),'full_criterion_met':False,'reviews':[],
         'pending':'Two independent full-criterion reviews of exact source/preparation/clock/record contract',
         'evidence_files':{str((HERE/p).relative_to(ROOT)):h for p,h in files.items()}})
    print(json.dumps({'manifest_sha256':sha(HERE/'MANIFEST.json'),'files':len(files),
         'producer_pins':len(sources),'summary':summary},indent=2))

if __name__=='__main__':main()
