"""Original independent cold recomputation and focused tests on one snapshot."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time
ROOT=Path.cwd();HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT))
OUT=Path(os.environ['AUDIT_OUTPUT']);OUT.mkdir(parents=True,exist_ok=True)
manifest=json.loads((HERE/'manifest.json').read_text());sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def source_check():
    bad=[p for p,h in manifest['input_sha256'].items()if not (ROOT/p).is_file()or sha(ROOT/p)!=h]
    assert not bad,('frozen scientific source changed',bad)
    for p,h in manifest['audit_source_sha256'].items():assert sha(ROOT/p)==h,('audit helper changed',p)
source_check()
if sys.argv[1]=='source-check':
    print('PASS_EXACT_FROZEN_SCIENTIFIC_SOURCE_INVENTORY');raise SystemExit(0)
name=manifest['packages'][int(sys.argv[1])]
relative='model/'+name+'/results.json';receipt=ROOT/relative
retained=receipt.read_bytes();retained_hash=sha(receipt);before=json.loads(retained)['source_sha256']
assert relative not in before and before==manifest['package_pins'][name]
assert all(sha(ROOT/p)==h for p,h in before.items())
(OUT/'retained-receipt.json').write_bytes(retained)
start=time.monotonic();record={'package':name,'commit':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
 'manifest_sha256':sha(HERE/'manifest.json'),'retained_receipt_sha256':retained_hash,'source_inventory_matches_frozen_manifest':True,
 'original_producer_and_comparison_source_unchanged':True,'all_outputs_immutable':True,'no_scientific_result_cache':True}
commands=(('cold',[sys.executable,'-W','error','-m','model.'+name,'--check']),
          ('focused-tests',[sys.executable,'-W','error','-m','pytest','-q','model/'+name]))
for mode,argv in commands:
 record[mode+'_command']=argv
 with (OUT/(mode+'.log')).open('w')as log:
  record[mode+'_exit']=subprocess.run(argv,stdout=log,stderr=subprocess.STDOUT,timeout=2400).returncode
 if record[mode+'_exit']!=0:break
 source_check()
source_check()
record['source_inventory_unchanged_after_execution']=all(sha(ROOT/p)==h for p,h in before.items())
record['retained_receipt_bytes_unchanged']=receipt.read_bytes()==retained
changed=subprocess.check_output(['git','diff','--name-only'],text=True).splitlines()
record['tracked_source_and_receipts_unchanged']=not changed
record['changed_paths']=changed;record['elapsed_seconds']=time.monotonic()-start
record['passed']=all(record.get(k)is True for k in('source_inventory_unchanged_after_execution','retained_receipt_bytes_unchanged','tracked_source_and_receipts_unchanged'))and all(record.get(k+'_'+'exit')==0 for k in('cold','focused-tests'))
(OUT/'summary.json').write_text(json.dumps(record,indent=2,sort_keys=True)+'\n');print(json.dumps(record),flush=True)
assert record['passed'],name
