"""Bind the completed row candidate to the finished observer-only audit."""
from pathlib import Path
import hashlib
import json
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
AUDIT=ROOT/'openwave-tuning/next/CHAIN3/row1_integrity_audit'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(n,x):(HERE/n).write_text(json.dumps(x,indent=2)+'\n')

def main():
    result=json.loads((HERE/'RESULT.json').read_text());vacuum=json.loads((HERE/'VACUUM.json').read_text())
    assert all(result['checks'].values()) and all(result['controls'].values())
    source={**result['source_sha256'],**vacuum['source_sha256']}
    for name,h in source.items():assert sha(ROOT/name)==h,(name,'source changed')
    audit=json.loads((AUDIT/'MANIFEST.json').read_text())
    assert audit['verdict']=='VALIDATED_WITH_RECORD_ONLY_FACADE'
    source[str((AUDIT/'MANIFEST.json').relative_to(ROOT))]=sha(AUDIT/'MANIFEST.json')
    for name,h in audit['files_sha256'].items():
        p=AUDIT/name;assert sha(p)==h
        if str(p.relative_to(ROOT)) in source:assert source[str(p.relative_to(ROOT))]==h
        source[str(p.relative_to(ROOT))]=h
    data=json.loads((AUDIT/'RESULT.json').read_text())
    assert data['complete_joint_laws']['independently_regenerated_joint_maps']==768
    assert data['independent_pointer_laws']['every_effect_sum_equals_identity']
    assert data['observer_only_facade']['pointer_snapshot_defect']['facade_verdict']['compatible_pointer_initials']==[0,9]
    for name,h in json.loads((AUDIT/'PINNED_INPUTS.json').read_text())['files_sha256'].items():
        assert sha(ROOT/name)==h
        if name in source:assert source[name]==h
        source[name]=h
    criterion=result['criterion'];ch=hashlib.sha256(criterion.encode()).hexdigest()
    save('CRITERION.json',{'criterion':criterion,'criterion_sha256':ch,
        'hash_rule':'SHA256 of the exact UTF-8 criterion string, without appended newline',
        'regime':'static','upstream_file_sha256':sha(HERE/'upstream/MODELS.md'),
        'precedent_sha256':sha(HERE/'upstream/m5_1_winding.py'),
        'no_added_mass_spin_or_exchange_gate':True})
    save('CUSTODY.json',{'observer_only_audit_manifest_sha256':sha(AUDIT/'MANIFEST.json'),
         'source_sha256':source,'all_native_and_facade_bytes_match_executed_row_result':True,
         'independent_complete_laws':768,'effect_sum_is_identity':True,
         'pending_full_row_reviews':True})
    excluded={'MANIFEST.json','VALIDATION.json','STATUS.json'}
    files={str(p.relative_to(ROOT)):sha(p) for p in sorted(HERE.rglob('*'))
           if p.is_file() and p.name not in excluded and '__pycache__' not in p.parts}
    save('MANIFEST.json',{'schema':'oph.row01.full_criterion_candidate.manifest.v1',
         'candidate_id':result['candidate_id'],'files':files,'source_sha256':source,
         'criterion_sha256':ch,'full_criterion_met':False,
         'scope':'Concrete full-criterion candidate pending two producer-independent full-row decisions; no official status change.',
         'observer_only_audit_manifest_sha256':sha(AUDIT/'MANIFEST.json')})
    save('VALIDATION.json',{'full_criterion_met':False,'candidate_id':result['candidate_id'],
         'evidence_files':{**files,str((HERE/'MANIFEST.json').relative_to(ROOT)):sha(HERE/'MANIFEST.json')},
         'criterion_sha256':ch,'reviews':[],
         'pending':'Two producer-independent full-criterion reviews on this frozen source/readout contract.'})
    save('STATUS.json',{'row':1,'state':'awaiting_independent_full_row_reviews',
         'candidate_id':result['candidate_id'],'manifest_sha256':sha(HERE/'MANIFEST.json'),
         'full_criterion_met':False,'scientific_execution_complete':True,
         'external_dependency_wait':None,'next_action':'Root and a second producer-independent reviewer assess the full criterion; fix concrete findings additively.'})
    print(json.dumps({'manifest_sha256':sha(HERE/'MANIFEST.json'),'core_files':len(files),
         'source_pins':len(source),'observer_audit':sha(AUDIT/'MANIFEST.json'),'criterion_sha256':ch}))

if __name__=='__main__':main()
