"""Real-job altered saved-baseline and equal-length multi-record failure tests."""
import json
from pathlib import Path
import sys
import tempfile
from validate_descendants import HERE, SOURCE, invoke, snapshot


def main():
    output=HERE/'validation_supplement';output.mkdir(exist_ok=False);results=[]
    module=HERE/'descendants/verify_source.py'
    with tempfile.TemporaryDirectory(prefix='contract-extra-') as temp:
        folder=Path(temp);contracts=folder/'contracts.json'
        original=json.loads((HERE/'descendants/SOURCE_JOB_CONTRACTS.json').read_text())
        changed=original['jobs'][23]['expected_output'][0]['rows'][0]
        # Select the actual bounded energy interval, never a status flag.
        def change_number(value):
            if isinstance(value,dict):
                for key,item in value.items():
                    if isinstance(item,float):value[key]=item+.1;return True
                    if isinstance(item,(dict,list)) and change_number(item):return True
            elif isinstance(value,list):
                for i,item in enumerate(value):
                    if isinstance(item,float):value[i]=item+.1;return True
                    if isinstance(item,(dict,list)) and change_number(item):return True
            return False
        if not change_number(changed):raise RuntimeError('No saved scientific reading')
        contracts.write_text(json.dumps(original));before=snapshot(folder)
        results.append(invoke([sys.executable,'-B','-O',str(module),'--source',str(SOURCE),
                               '--contracts',str(contracts),'--job','23'],folder,'real_job_changed_saved_baseline',output,False))
        if before!=snapshot(folder):raise RuntimeError('Harness rewrote its saved contract')
        expected=[{'passed':True},{'passed':True}]
        contract={'expected_output':expected,'runtime_paths':[],'required_verdict_paths':[[0,'passed'],[1,'passed']]}
        script=folder/'earlier.py';script.write_text('import runpy\nf=runpy.run_path('+repr(str(module))+')["check_output"]\nf('+repr([{'passed':False},{'passed':True}])+','+repr(contract)+')\n')
        result=invoke([sys.executable,'-B','-O',str(script)],folder,'equal_length_earlier_failure',output,False)
        if 'A current result declares failure' not in result['stderr']:raise RuntimeError('Wrong rejection path')
        results.append(result)
        target=folder/'existing.json';target.write_text('existing record\n');before=target.read_bytes()
        results.append(invoke([sys.executable,'-B',str(module),'--source',str(SOURCE),
                               '--job','23','--record-to',str(target)],folder,'harness_record_overwrite',output,False))
        if target.read_bytes()!=before:raise RuntimeError('Harness overwrote receipt')
    summary={'passed':all(r['passed'] for r in results),'results':[{k:v for k,v in r.items() if k not in ('stdout','stderr')} for r in results]}
    (output/'SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
    raise SystemExit(0 if summary['passed'] else 1)


if __name__=='__main__':main()
