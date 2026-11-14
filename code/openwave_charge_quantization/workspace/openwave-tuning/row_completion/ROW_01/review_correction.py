"""Test and seal the additive first-review correction, preserving the base."""
from copy import deepcopy
from pathlib import Path
import hashlib
import json
import topology
import reviewed
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(n,x):(HERE/n).write_text(json.dumps(x,indent=2)+'\n')

def main():
    base=HERE/'MANIFEST.json';expected='c4f6db87b1d2c7ddd4caa4c8ac5a95f09c33679f060d50845aa0bf42f5e54d12'
    assert sha(base)==expected
    manifest=json.loads(base.read_text())
    for key in ('files','source_sha256'):
        for relative,h in manifest[key].items():assert sha(ROOT/relative)==h
    shear=[[int(i==j) for j in range(11)] for i in range(11)];shear[0][1]=1
    assert topology.determinant(shear)==1 and shear[0][1]!=0
    selected={}
    for row in json.loads((HERE/'EXAMPLES.json').read_text()):
        key=(row['topological_degree'],row['claim']['meter']['reason'])
        selected.setdefault(key,row)
    out=[]
    for key,row in selected.items():
        verdict=reviewed.admit_and_decode(row['claim'],row['G'],row['region'],2)
        assert verdict['admitted'] and verdict['raw_spatial_pointer_degree']==row['raw_pointer_record_degree']
        out.append({'degree':key[0],'outcome':key[1],'admitted':True,'raw_degree':verdict['raw_spatial_pointer_degree']})
    vacuum=json.loads((HERE/'VACUUM.json').read_text())
    good=reviewed.admit_and_decode(vacuum['claim'],vacuum['G'],(0,),2)
    assert good['admitted'] and good['readout']==0
    bad=deepcopy(vacuum['claim']);bad['meter']['records']=[]
    rejected=reviewed.admit_and_decode(bad,vacuum['G'],(0,),2)
    assert not rejected['admitted']
    receipt={'base_manifest_sha256':expected,'abstract_GL11_shear':shear,
         'shear_determinant':1,'abstract_dictionary_uniqueness_claim':False,
         'native_source_and_local_B_dictionary_changed':False,'guarded_decoder_cases':out,
         'genuine_zero_readout':good,'empty_success_rejected':rejected,
         'base_and_provider_bytes_unchanged':True}
    save('REVIEW_CORRECTION_RECEIPT.json',receipt)
    names=['REVIEW_CORRECTION_01.md','reviewed.py','review_correction.py','REVIEW_CORRECTION_RECEIPT.json']
    files={str((HERE/n).relative_to(ROOT)):sha(HERE/n) for n in names}
    save('REVIEW_FIX_MANIFEST.json',{'base_manifest_sha256':expected,'files':files,
         'criterion_sha256':manifest['criterion_sha256'],'supersedes':'PROOF.md §2 abstract dictionary uniqueness sentence; raw decoder is post-admission only',
         'scientific_results_changed':False,'full_criterion_met':False})
    validation=json.loads((HERE/'VALIDATION.json').read_text())
    validation['evidence_files'].update(files)
    validation['evidence_files'][str((HERE/'REVIEW_FIX_MANIFEST.json').relative_to(ROOT))]=sha(HERE/'REVIEW_FIX_MANIFEST.json')
    save('VALIDATION.json',validation)
    status=json.loads((HERE/'STATUS.json').read_text())
    status['review_correction_manifest_sha256']=sha(HERE/'REVIEW_FIX_MANIFEST.json')
    save('STATUS.json',status)
    print(json.dumps({'review_fix_manifest_sha256':sha(HERE/'REVIEW_FIX_MANIFEST.json'),
          'guarded_decoder_cases':len(out),'empty_success_rejected':True}))

if __name__=='__main__':main()
