"""Authenticate seven upstream columns and rerun three exact implementations."""
from collections import Counter
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys
import tempfile

HERE=Path(__file__).resolve().parent
CONTRACT=json.loads((HERE/'contract.json').read_text())
MODELS=('M5','M7','M4','M6','M8','M9','M10')
ICONS={'✅':'validated','⚠️':'partial','❌':'negative','🚧':'untested'}

def rows(text):
    return [[v.strip()for v in line.strip().strip('|').split('|')]for line in text.splitlines()if line.startswith('|')]

def parse(text):
    block=text.split('### Summary Status',1)[1].split('### Simplest Test per Criterion',1)[0]
    out=[]
    for cells in rows(block):
        if not cells[0]or cells[0].startswith(('**','Criteria','---')):continue
        if len(cells)!=9 or cells[-1]not in('static','dynamic','both') or any(v not in ICONS for v in cells[1:-1]):
            raise ValueError('upstream matrix shape or status changed')
        out.append({'id':len(out)+1,'criterion':cells[0],'regime':cells[-1],
                    'upstream_status':dict(zip(MODELS,cells[1:-1]))})
    if len(out)!=31 or len({x['criterion']for x in out})!=31:raise ValueError('all31distinct criteria required')
    counts={m:{s:sum(ICONS[r['upstream_status'][m]]==s for r in out)for s in ICONS.values()}for m in MODELS}
    for v in counts.values():v['scored']=v['validated']+v['partial']
    summary=text.split('### Summary Count',1)[1].split('### Summary Status',1)[0]
    seen=set()
    for cells in rows(summary):
        icon=next((i for i in ICONS if cells[0].startswith(i)),None)
        if icon:
            if icon in seen or len(cells)!=8 or list(map(int,cells[1:]))!=[counts[m][ICONS[icon]]for m in MODELS]:
                raise ValueError('summary count differs from complete matrix')
            seen.add(icon)
    if seen!=set(ICONS):raise ValueError('complete tally required')
    key=lambda m:(-counts[m]['scored'],-counts[m]['validated'],counts[m]['negative'])
    if tuple(sorted(MODELS,key=key))!=MODELS:raise ValueError('column order differs from ranking rule')
    return {'criteria':out,'upstream_counts':counts,'upstream_column_order':list(MODELS),
        'oph_admitted':True,'oph_model_id':'M10','oph_official_scored_rows':counts['M10']['scored'],
        'oph_official_rank':1+sum(key(m)<key('M10')for m in MODELS),
        'ranking_rule':'validated+partial descending, validated descending, negative ascending'}

def authenticate(root):
    root=Path(root).resolve()
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip()
    if head!=CONTRACT['upstream_commit']:raise ValueError('upstream commit differs')
    for name,digest in CONTRACT['upstream_sha256'].items():
        if hashlib.sha256((root/name).read_bytes()).hexdigest()!=digest:raise ValueError('upstream source changed: '+name)
    return root

def run(root):
    root=authenticate(root);status=parse((root/'MODELS.md').read_text());records=[]
    base='openwave/xperiments/m10_oph/research/m10_1/returns/'
    with tempfile.TemporaryDirectory(prefix='oph-m10-recompute-')as directory:
        for label,script,output in [('solver','solve.py','results.json'),('audit_a','audit.py','audit_results.json'),('audit_b','audit.py','audit_results.json')]:
            work=Path(directory)/label;work.mkdir();source=root/(base+label+'/'+script)
            result=subprocess.run([sys.executable,str(source)],cwd=work,env={**os.environ,'GIT_OPTIONAL_LOCKS':'0'},capture_output=True,text=True)
            if result.returncode:raise RuntimeError(result.stderr[-3000:])
            actual=json.loads((work/output).read_text());expected=json.loads((root/(base+label+'/'+output)).read_text())
            if actual!=expected:raise ValueError('complete original implementation output differs: '+label)
            records.append({'implementation':label,'command':'python '+base+label+'/'+script+' (in an empty output directory)',
                            'exit_code':result.returncode,'complete_output_matches_upstream':True,
                            'result_sha256':hashlib.sha256((work/output).read_bytes()).hexdigest()})
    authenticate(root)
    adjudication=json.loads((root/'openwave/xperiments/m10_oph/research/data/m10_1_adjudication.json').read_text())
    if len(adjudication['rows'])!=7 or any(x['verdict']!='agree'for x in adjudication['rows']):raise ValueError('certification comparison differs')
    return {'schema':'oph.m10_current_review_gate.v1','verdict':'PASS_CURRENT_M10_STATUS_AND_THREE_EXACT_RECOMPUTATIONS',
        'upstream_commit':CONTRACT['upstream_commit'],'upstream_sha256':CONTRACT['upstream_sha256'],
        'official_status':status,'actual_recomputations':records,
        'upstream_certification':{'task':'M10.1','compared_values':adjudication['rows'],'criterion_cells_awarded':0,
            'scope':'Independent-method exact arithmetic on the declared twelve-direction scalar dispersion through degree six. No physical identification or criterion award.'},
        'prediction_exposure':CONTRACT['prediction_exposure'],
        'native_process_acceptance_gate':CONTRACT['gate'],
        'scope':{'new_upstream_submission':False,'strict_blindness_claimed':False,'new_observational_prediction':False,
            'new_Nature_validation_rows':0,'local_declared_source_checks_automatically_award_cells':False},
        'required_join':'Freeze one twelve-port repair/memory program, apparatus and coarse readout; implement the criterion observable through its actual events and compare the complete results, including failures, to the predeclared physical gate. Source-family labels or equivalent spectra alone do not establish this implementation.'}
