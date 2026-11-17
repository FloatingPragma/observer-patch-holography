"""Real-process rejection and read-only checks for the producer's certificate CLI."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time

HERE=Path(__file__).resolve().parent
PRODUCER=HERE.parents[1]/'work/2026-11-16_codex_interval_source_certificate'
ENV=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')


def source_snapshot():
    names=['certificate.py','interval_math.py','interval_reader.py','HEADLINE.json','CERTIFICATE.json']
    return {name:(hashlib.sha256((PRODUCER/name).read_bytes()).hexdigest(),(PRODUCER/name).stat().st_mtime_ns) for name in names}


def main():
    output=HERE/'processes';output.mkdir(exist_ok=False);before=source_snapshot();results=[]
    with tempfile.TemporaryDirectory(prefix='interval-certificate-cli-') as temp:
        folder=Path(temp)
        commands=[('readonly_headline',[sys.executable,'-B',str(PRODUCER/'certificate.py')],0)]
        original=(PRODUCER/'certificate.py').read_text()
        for label,old,new in [('impurity_gate',"single_off < real('0.002')","single_off < real('-1')"),
                              ('wrong_source_spectrum','for lam in (4, 6):','for lam in (5, 6):')]:
            if original.count(old)!=1:raise AssertionError('Unique scientific mutation site')
            code=original.replace(old,new);script=folder/(label+'.py')
            script.write_text('import sys\nsys.path.insert(0,'+repr(str(PRODUCER))+')\nexec(compile('+repr(code)+','+repr(str(PRODUCER/'certificate.py'))+',"exec"),{"__name__":"__main__","__file__":'+repr(str(PRODUCER/'certificate.py'))+'})\n')
            commands.append((label,[sys.executable,'-B','-O',str(script),'--k','0','--references','stationary','--record-to',str(folder/(label+'.json'))],1))
        changed=folder/'changed.json';saved=json.loads((PRODUCER/'HEADLINE.json').read_text());saved['cases'][0]['single_off_band']['lo_dyadic'][1]='1';changed.write_text(json.dumps(saved))
        commands.append(('changed_saved_certificate',[sys.executable,'-B',str(PRODUCER/'certificate.py'),'--verify',str(changed)],1))
        existing=folder/'existing.json';existing.write_text('retained expectation\n')
        commands.append(('refused_overwrite',[sys.executable,'-B',str(PRODUCER/'certificate.py'),'--record-to',str(existing)],1))
        for label,command,expected in commands:
            start=time.monotonic();run=subprocess.run(command,capture_output=True,text=True,env=ENV,timeout=60)
            result={'label':label,'command':command,'exit_code':run.returncode,'stdout':run.stdout,
                    'stderr':run.stderr,'seconds':time.monotonic()-start,'passed':run.returncode==expected}
            if expected and 'ValueError' not in run.stderr:result['passed']=False
            (output/(label+'.json')).write_text(json.dumps(result,indent=2)+'\n');results.append(result)
            print(json.dumps({k:v for k,v in result.items() if k not in ('stdout','stderr')}),flush=True)
        if existing.read_text()!='retained expectation\n':raise AssertionError('Overwrote saved record')
    summary={'passed':all(x['passed'] for x in results),'producer_saved_data_unchanged':before==source_snapshot(),
             'results':[{k:v for k,v in r.items() if k not in ('stdout','stderr')} for r in results],
             'credit':'Execution/integrity only; spectrum corruption is not a native competing law'}
    (HERE/'PROCESS_AUDIT.json').write_text(json.dumps(summary,indent=2)+'\n')
    if not summary['passed'] or not summary['producer_saved_data_unchanged']:raise AssertionError('Process audit failed')


if __name__=='__main__':main()
