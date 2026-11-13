"""Recompute live theorem navigation without altering finite scientific acceptance.

The exact accepted scientific contents are authenticated by semantic digests.
Only three navigation fields and their explicit custody consumers may change.
The retained dated-authority replayers and their contracts are untouched.
"""
from __future__ import annotations
import argparse,copy,hashlib,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import tools as tool_package
sys.path.insert(0,next(iter(tool_package.__path__)))
from tools import strict_json
from tools.build_model_theorem_map import validate as validate_navigation
from model.completion.finite_attachment_map import validate as validate_map
from model.completion.program_class_v2 import validate as validate_class

MAP='model/completion/finite_attachment_map.json'
REVIEW='model/completion/program_class_v2_review.json'
CLASS='model/completion/program_class_v2.json'
PROGRAM='claims/completion_program.json'
RECORD='model/completion/navigation_custody_current.json'
SELF='tools/renew_finite_attachment_navigation.py'
FILES=(MAP,REVIEW,CLASS,PROGRAM)
FIELDS=('claim_count','lean_module_count','contract_coverage_sha256')
BASELINE_SEMANTIC_SHA256={
  "model/completion/finite_attachment_map.json": "ae0273e71cd22e0d182a16f06c31385413c12a4d2a56aabc6d53d6859a914f52",
  "model/completion/program_class_v2_review.json": "bd339fe9505fe964ad95429c389c49ac94106f43bbf3148aa58254345e62e0a6",
  "model/completion/program_class_v2.json": "c367a61f93d30aa5e1cb91ac4f2e7e17f582a6b31278748d42148b4ac9ae3b9a",
  "claims/completion_program.json": "0f84541bc77d3f5b8bfdbcae9412ffb629a455d3889d2d7c0238655581cc0153"
}

def canonical(x):return json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n'
def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def raw_digest(x):return hashlib.sha256(x.encode()).hexdigest()
def require(ok,message):
    if not ok:raise ValueError(message)

def pointer_update(value,digests):
    if isinstance(value,dict):
        if value.get('path')in digests and 'sha256'in value:value['sha256']=digests[value['path']]
        for child in value.values():pointer_update(child,digests)
    elif isinstance(value,list):
        for child in value:pointer_update(child,digests)

def semantic_views(documents):
    views=copy.deepcopy(documents)
    for k in FIELDS:views[MAP]['existing_navigation_coverage'][k]='<navigation>'
    views[REVIEW]['source_bindings'][MAP]='<map custody>'
    views[REVIEW].pop('navigation_coverage_custody',None)
    views[CLASS]['adoption']['independent_review']['sha256']='<review custody>'
    pointer_update(views[PROGRAM],{p:'<navigation custody>'for p in (MAP,REVIEW,CLASS)})
    return views

def unchanged_science(documents):
    for path,value in semantic_views(documents).items():
        require(digest(value)==BASELINE_SEMANTIC_SHA256[path],'scientific contents changed outside navigation custody: '+path)

def build(root=ROOT):
    root=Path(root).resolve()
    documents={p:strict_json.load(root/p)for p in FILES};unchanged_science(documents)
    register=strict_json.load(root/'tracking/premise_register.json')
    claims=strict_json.load(root/'claims/claim_registry.yaml')
    errors=validate_navigation(register['model_supply'],register,claims,root)
    require(not errors,'invalid complete typed navigation: '+repr(errors))
    coverage=register['model_supply']['contract_coverage']
    navigation={'claim_count':len(coverage['claim_ids']),'lean_module_count':len(coverage['lean_files']),'contract_coverage_sha256':digest(coverage)}
    documents[MAP]['existing_navigation_coverage'].update(navigation)
    acceptance=validate_map(documents[MAP],root)
    require(acceptance['canonical_premises_discharged']==0,'navigation cannot discharge physical premises')
    map_hash=raw_digest(canonical(documents[MAP]))
    record={'schema':'oph.finite_attachment_navigation_custody.v1',
      'role':'CURRENT_NAVIGATION_ONLY_UNCHANGED_SCIENTIFIC_ACCEPTANCE',
      'producer':{'path':SELF,'sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
      'navigation':navigation,'finite_attachment_map_sha256':map_hash,
      'scientific_semantics_sha256':BASELINE_SEMANTIC_SHA256,
      'original_map_validator_result':acceptance,
      'gate_statuses_changed':False,'physical_premises_discharged':[],
      'retained_authority_contracts_modified':False}
    documents[REVIEW]['source_bindings'][MAP]=map_hash
    documents[REVIEW]['navigation_coverage_custody']={
      'path':RECORD,'sha256':raw_digest(canonical(record)),
      'role':'DOCUMENTARY_NAVIGATION_ONLY_NO_ADDITIONAL_THEOREM_OR_PHYSICAL_ACCEPTANCE'}
    documents[CLASS]['adoption']['independent_review']['sha256']=raw_digest(canonical(documents[REVIEW]))
    pointer_update(documents[PROGRAM],{p:raw_digest(canonical(documents[p]))for p in (MAP,REVIEW,CLASS)})
    unchanged_science(documents)
    return {**{p:canonical(d)for p,d in documents.items()},RECORD:canonical(record)}

def execute(root=ROOT,write=False):
    root=Path(root).resolve();outputs=build(root)
    before={p:(root/p).read_bytes()if(root/p).exists()else None for p in outputs}
    changed=[p for p,text in outputs.items()if before[p]!=text.encode()]
    if write:
        try:
            for p,text in outputs.items():(root/p).write_text(text)
            require(build(root)==outputs,'cold navigation replay differs')
            validate_class(strict_json.load(root/CLASS),root)
        except BaseException:
            for p,raw in before.items():
                if raw is None:(root/p).unlink(missing_ok=True)
                else:(root/p).write_bytes(raw)
            raise
    return {'changed_paths':changed,'source_inventory':json.loads(outputs[RECORD])['navigation'],
            'scientific_semantics_exactly_preserved':True,'written':write}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=ROOT)
    m=p.add_mutually_exclusive_group(required=True);m.add_argument('--write',action='store_true');m.add_argument('--check',action='store_true');m.add_argument('--preview',action='store_true')
    a=p.parse_args();result=execute(a.root,a.write)
    if a.check:
        require(not result['changed_paths'],'stale current navigation custody')
        validate_class(strict_json.load(a.root/CLASS),a.root.resolve())
    print(json.dumps(result,indent=2,sort_keys=True))
if __name__=='__main__':main()
