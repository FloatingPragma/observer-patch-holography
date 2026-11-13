"""Cold execution against immutable source bytes, with no receipt writers."""
from pathlib import Path
import hashlib, importlib, json, os, subprocess, sys, time
ROOT=Path.cwd(); HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
OUT=Path(os.environ['AUDIT_OUTPUT']); OUT.mkdir(parents=True,exist_ok=True)
manifest=json.loads((HERE/'manifest.json').read_text())
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def source_check():
    bad=[p for p,h in manifest['input_sha256'].items() if not (ROOT/p).is_file() or sha(ROOT/p)!=h]
    assert not bad, ('immutable source changed',bad)
source_check()
if sys.argv[1]=='source-check':
    print('PASS_EXACT_FROZEN_SCIENTIFIC_SOURCE_INVENTORY')
    sys.exit(0)
name=manifest['packages'][int(sys.argv[1])]
module=importlib.import_module('model.'+name+'.__main__')
receipt=ROOT/'model'/name/'results.json'; retained=sha(receipt)
before=module.pins()
assert before==manifest['package_pins'][name]
assert before==json.loads(receipt.read_text())['source_sha256']
start=time.monotonic()
record={'package':name,'commit':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
    'manifest_sha256':sha(HERE/'manifest.json'),'retained_receipt_sha256':retained,
    'source_inventory_matches_frozen_manifest':True,'no_writer_invoked':True}
argv=[sys.executable,'-m','model.'+name,'--check']; record['actual_cold_command']=argv
with (OUT/'original-cold.log').open('w') as log:
    record['original_cold_exit']=subprocess.run(argv,stdout=log,stderr=subprocess.STDOUT,timeout=2400).returncode
test=[sys.executable,'-W','error','-m','pytest','-q','model/'+name]
record['actual_test_command']=test
with (OUT/'focused-tests.log').open('w') as log:
    record['focused_tests_exit']=subprocess.run(test,stdout=log,stderr=subprocess.STDOUT,timeout=2400).returncode
source_check()
record['source_inventory_unchanged_after_execution']=module.pins()==before
record['retained_receipt_bytes_unchanged']=sha(receipt)==retained
clean=subprocess.run(['git','diff','--exit-code'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
(OUT/'tracked-worktree-diff.txt').write_text(clean.stdout)
record['tracked_worktree_unchanged']=clean.returncode==0
record['elapsed_seconds']=time.monotonic()-start
record['passed']=all((record['original_cold_exit']==0,record['focused_tests_exit']==0,
    record['source_inventory_unchanged_after_execution'],record['retained_receipt_bytes_unchanged'],record['tracked_worktree_unchanged']))
(OUT/'summary.json').write_text(json.dumps(record,indent=2,sort_keys=True)+'\n')
print(json.dumps(record),flush=True)
assert record['passed'], name
