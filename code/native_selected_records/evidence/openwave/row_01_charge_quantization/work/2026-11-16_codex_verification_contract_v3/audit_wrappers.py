"""Exercise real scientific predicates of every distinct historical replay entry.

Only disposable copies are edited. Repeated --k instances share one audited
entry; their unmutated seven-domain execution remains in the 28-job replay.
"""
import ast
import hashlib
import json
from pathlib import Path
import shutil
import sys
import tempfile
from validate_descendants import HERE, SOURCE, invoke

JOBS=json.loads((HERE/'INPUTS.json').read_text())['jobs']
# Explicit predicate fragments select scientific comparisons, never SHA checks.
SITES={
    1:'r.stdout==expected',2:"len(old['source_records'])==30",3:'max(abs(np.array(planted)-np.array(actual)-.1))<1e-6',
    4:'np.max(abs(values-np.array([0.,4.,6.])))<1e-12',
    5:'visibility < 0 if law == "native" and turns == 1 else visibility > 0',
    6:'preparation == F(case["preparation"]["accepted_probability"])',
    8:'max(errors) < 2e-10',9:'v0 == F(case["exact_visibility"])',
    10:"json.dumps(fresh,sort_keys=True,separators=(',',':'))==json.dumps(recorded,sort_keys=True,separators=(',',':'))",
    11:'-1e-10<=eta<=ETA',12:'len(interior) == 2',13:'json.loads(json.dumps(fresh))==stored',
    14:'error<tol',21:"total==Fraction(s['success_probability']) and total+Fraction(s['timeout_probability'])==1",
    22:"result['unique_component']",23:'energy[1]<OLD_PAIR',24:'residual < tolerance',
    25:'stable(actual)==stable(filed)',27:"np.linalg.norm(np.einsum('rxa,rxb->ab',F.conj(),F)-np.eye(4))<1e-9",
}


def mutate_predicate(text, fragment):
    # Locate the actual assert or unconditional require/need call, not a report bit.
    candidates=[]
    for node in ast.walk(ast.parse(text)):
        if isinstance(node,ast.Assert):candidates.append(node.test)
        if isinstance(node,ast.Call) and isinstance(node.func,ast.Name) and node.func.id in ('require','need'):
            if node.args:candidates.append(node.args[0])
    matches=[n for n in candidates if fragment in ast.get_source_segment(text,n)]
    if len(matches)!=1:raise RuntimeError('Ambiguous actual predicate: '+fragment)
    old=ast.get_source_segment(text,matches[0])
    return text.replace(old,'not ('+old+')',1)


def main():
    output=HERE/'wrapper_processes';output.mkdir(exist_ok=False);results=[]
    with tempfile.TemporaryDirectory(prefix='contract-wrapper-mutants-') as temp:
        copied=Path(temp)/'source';shutil.copytree(SOURCE,copied)
        for index,fragment in SITES.items():
            directory,name,args=JOBS[index];folder=copied/'evidence/openwave'/directory;script=folder/name
            original=script.read_bytes();script.write_text(mutate_predicate(original.decode(),fragment))
            try:
                result=invoke([sys.executable,'-B',name,*args],folder,f'job_{index:02}_scientific_predicate',output,False)
                # A custody rejection or timeout is not successful gate coverage.
                result['mutation_fragment']=fragment
                if result.get('timeout') or 'MANIFEST.sha256' in result.get('stderr',''):
                    result['passed']=False
                results.append(result)
            finally:script.write_bytes(original)
        # Reservoir wrapper binds itself prospectively. Update only that test-copy
        # prereg entry so the intended scientific comparison, not SHA rejection, runs.
        folder=copied/'evidence/openwave'/JOBS[0][0];script=folder/'headline.py'
        original=script.read_text();script.write_text(mutate_predicate(original,"receipt['passed']"))
        prereg=folder/'PREREG.json';data=json.loads(prereg.read_text())
        data['files']['headline.py']=hashlib.sha256(script.read_bytes()).hexdigest();prereg.write_text(json.dumps(data))
        results.append(invoke([sys.executable,'-B','headline.py'],folder,'reservoir_primary_gate',output,False))
        # The custom-check menu audit has no assertions. Change its actual required
        # reconstructed transfer operand, retaining real execution and error path.
        directory,name,args=JOBS[26];folder=copied/'evidence/openwave'/directory;script=folder/name
        text=script.read_text();old="near(P,np.eye(3)-eps**2*";new="near(2*P,np.eye(3)-eps**2*"
        if text.count(old)!=1:raise RuntimeError('Unique menu operator check missing')
        script.write_text(text.replace(old,new));results.append(invoke([sys.executable,'-B',name,*args],folder,'menu_transfer_check',output,False))
    summary={'passed':all(r['passed'] for r in results),'distinct_entry_mutations':len(results),
             'repeated_parameter_domains':'independent.py k1..6 and alternate spin receipt share audited code; all execute unmutated in full replay',
             'test_copy_prereg_rebinding':'Reservoir headline.py only; original immutable registration unchanged.',
             'scientific_emergence_credit':False,'results':[{k:v for k,v in r.items() if k not in ('stdout','stderr')} for r in results]}
    (HERE/'WRAPPER_AUDIT.json').write_text(json.dumps(summary,indent=2)+'\n')
    raise SystemExit(0 if summary['passed'] else 1)


if __name__=='__main__':main()
