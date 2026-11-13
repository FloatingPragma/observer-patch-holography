"""Recompute live theorem navigation without altering finite scientific acceptance.

The exact accepted scientific contents are authenticated by semantic digests.
Only three navigation fields, one authenticated candidate presentation-source
pin, and their explicit custody consumers may change. The retained
dated-authority replayers and their contracts are untouched.
"""
from __future__ import annotations
import argparse,ast,copy,hashlib,json,sys
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

CANDIDATE='model/particle_joint_reconstruction/candidate.py'
CANDIDATE_OLD_SHA256='3727f55b3e691b535350a9c47fc8cb43c01aeedf081b11e54816896e0c90e10f'
CANDIDATE_CURRENT_SHA256='64d72fd61333566c3ca3acc7e09d2de433910e774b62bb18e1d61c8e9f7f9582'
CANDIDATE_PRESENTATION_ADDITION=(
    b'    # Exact presentation classification only; scientific proofs stay bound.\n'
    b'    paths.difference_update(gap.PRESENTATION_INPUTS)\n')
PRESENTATION_PARENT='model/particle_yang_mills_gap_contract/transfer.py'
PRESENTATION_PARENT_SHA256='f3d4ce01607a6b781b0e9f704573d72f1ffb0dd245b9f192e5b9a15a54a9b8d2'


def candidate_presentation_custody(documents,root=ROOT):
    """Authenticate one exact presentation edit against the retained raw source.

    Removing precisely these two lines recovers the prior cryptographic digest;
    hence every other byte, including every operative expression, is identical.
    The sole added AST expression is in _bindings, and its referenced literal is
    exactly one paper. No arbitrary source, row, old pin or digest is exempted.
    The original semantic digest continues to authenticate every other field.
    """
    root=Path(root)
    program=documents[PROGRAM]
    try:stage=program['stages'][6];row=stage['progress_evidence'][38]
    except (KeyError,IndexError,TypeError)as error:
        raise ValueError('candidate presentation custody row missing')from error
    require(stage.get('id')=='CP-05' and row.get('id')=='conditional_joint_chiral_candidate'
            and row.get('path')==CANDIDATE
            and row.get('sha256')==CANDIDATE_CURRENT_SHA256,
            'candidate presentation custody identity or current hash changed')
    raw=(root/CANDIDATE).read_bytes()
    require(hashlib.sha256(raw).hexdigest()==CANDIDATE_CURRENT_SHA256
            and raw.count(CANDIDATE_PRESENTATION_ADDITION)==1,
            'candidate presentation source changed')
    retained=raw.replace(CANDIDATE_PRESENTATION_ADDITION,b'',1)
    require(hashlib.sha256(retained).hexdigest()==CANDIDATE_OLD_SHA256,
            'candidate presentation inverse does not recover retained source')
    tree=ast.parse(raw);old_tree=ast.parse(retained)
    functions=[node for node in tree.body if isinstance(node,ast.FunctionDef)and node.name=='_bindings']
    require(len(functions)==1,'candidate presentation binding owner changed')
    expected=ast.parse('paths.difference_update(gap.PRESENTATION_INPUTS)').body[0]
    matches=[node for node in functions[0].body
             if ast.dump(node,include_attributes=False)==ast.dump(expected,include_attributes=False)]
    require(len(matches)==1,'candidate presentation AST edit changed')
    functions[0].body.remove(matches[0])
    require(ast.dump(tree,include_attributes=False)==ast.dump(old_tree,include_attributes=False),
            'candidate operative AST changed')
    parent_raw=(root/PRESENTATION_PARENT).read_bytes()
    require(hashlib.sha256(parent_raw).hexdigest()==PRESENTATION_PARENT_SHA256,
            'candidate presentation parent implementation changed')
    parent=ast.parse(parent_raw)
    assignments=[node for node in ast.walk(parent)if isinstance(node,ast.Assign)
                 and any(isinstance(target,ast.Name)and target.id=='PRESENTATION_INPUTS'
                         for target in node.targets)]
    require(len(assignments)==1 and assignments[0]in parent.body
            and ast.literal_eval(assignments[0].value)==('extra/yang_mills_gap_clay_problem.tex',),
            'candidate presentation parent scope changed')


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
    views[PROGRAM]['stages'][6]['progress_evidence'][38]['sha256']=CANDIDATE_OLD_SHA256
    return views

def unchanged_science(documents,root=ROOT):
    candidate_presentation_custody(documents,root)
    for path,value in semantic_views(documents).items():
        require(digest(value)==BASELINE_SEMANTIC_SHA256[path],'scientific contents changed outside navigation custody: '+path)

def build(root=ROOT):
    root=Path(root).resolve()
    documents={p:strict_json.load(root/p)for p in FILES};unchanged_science(documents,root)
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
    unchanged_science(documents,root)
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
