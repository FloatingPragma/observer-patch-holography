"""Verify every retained scientific receipt; any failed check fails the job."""
from pathlib import Path
import hashlib,importlib,json,os,subprocess,sys,time,traceback
ROOT=Path.cwd();sys.path.insert(0,str(ROOT));HERE=Path(__file__).resolve().parent
from model.verify_particle_extensions import PACKAGES
OUT=Path(os.environ['AUDIT_OUTPUT']);OUT.mkdir(parents=True,exist_ok=True)
manifest=json.loads((HERE/'manifest.json').read_text());index=int(sys.argv[1]);names=manifest['replay_packages'][index::8]
assert len(PACKAGES)==48 and list(PACKAGES)==manifest['packages']
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def source_check():
    errors=[p for p,h in manifest['input_sha256'].items()if not(ROOT/p).is_file()or sha(ROOT/p)!=h]
    assert not errors,('immutable scientific inputs changed',errors)
source_check()
summary={'shard':index,'commit':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),'manifest_sha256':sha(HERE/'manifest.json'),'input_files':len(manifest['input_sha256']),'packages':[],'no_writer_invoked':True,'scope':'Only changed or previously failed package colds; unchanged PASS evidence is separately authenticated'}
for name in names:
    out=OUT/name;out.mkdir();record={'package':name,'retained_receipt_sha256':sha(ROOT/'model'/name/'results.json')}
    start=time.monotonic()
    try:
        module=importlib.import_module('model.'+name+'.__main__')
        before=module.pins();expected=manifest['package_pins'][name]
        record['source_inventory_matches_frozen_manifest']=before==expected
        receipt=json.loads((ROOT/'model'/name/'results.json').read_text())
        record['retained_source_inventory_current']=receipt['source_sha256']==before
        assert record['source_inventory_matches_frozen_manifest'] and record['retained_source_inventory_current']
        argv=[sys.executable,'-m','model.'+name,'--check']
        if name=='openwave_boosted_records':argv+=['--openwave-root',os.environ['OPH_OPENWAVE_ROOT'],'--python',os.environ['OPH_TAICHI_PYTHON']]
        record['actual_command']=argv
        with (out/'original-cold.log').open('w')as log:
            try:record['original_cold_exit']=subprocess.run(argv,stdout=log,stderr=subprocess.STDOUT,timeout=2400).returncode
            except subprocess.TimeoutExpired:record['original_cold_exit']=124
        record['source_inventory_unchanged_after_execution']=module.pins()==before
        record['retained_receipt_bytes_unchanged']=sha(ROOT/'model'/name/'results.json')==record['retained_receipt_sha256']
        source_check()
    except Exception as error:
        record['exception']=repr(error);(out/'exception.txt').write_text(traceback.format_exc())
    record['elapsed_seconds']=time.monotonic()-start
    record['passed']=record.get('original_cold_exit')==0 and record.get('source_inventory_unchanged_after_execution',False) and record.get('retained_receipt_bytes_unchanged',False) and 'exception' not in record
    (out/'report.json').write_text(json.dumps(record,indent=2)+'\n');summary['packages'].append(record)
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(record),flush=True)
source_check()
clean=subprocess.run(['git','diff','--exit-code'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
(OUT/'tracked-worktree-diff.txt').write_text(clean.stdout)
summary['all_scheduled_original_checks_executed']=len(summary['packages'])==len(names) and all('original_cold_exit'in r for r in summary['packages'])
summary['all_checks_passed']=all(r['passed']for r in summary['packages'])
summary['all_immutable_source_bytes_unchanged']=True
summary['tracked_worktree_unchanged']=clean.returncode==0
summary['failed_packages']=[r['package']for r in summary['packages']if not r['passed']]
(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary),flush=True)
assert summary['all_scheduled_original_checks_executed'] and summary['all_checks_passed'] and summary['tracked_worktree_unchanged'],summary['failed_packages']
