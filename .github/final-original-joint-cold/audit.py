from pathlib import Path
import hashlib,json,subprocess,sys,os
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent;OUT=ROOT/'audit-artifacts';OUT.mkdir(exist_ok=True)
manifest=json.loads((HERE/'manifest.json').read_text());sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
for key in ('input_sha256','snapshot_sha256'):
 for p,h in manifest[key].items():assert sha(ROOT/p)==h,(key,p)
assert sha(ROOT/manifest['target_output'])==manifest['retained_receipt_sha256']
report={'phase':sys.argv[1],'commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'all_required_inputs_exact':True,'input_files':len(manifest['input_sha256']),'all_snapshot_files_exact':True,'snapshot_files':len(manifest['snapshot_sha256']),'retained_receipt_sha256':sha(ROOT/manifest['target_output']),'no_writer_invoked':True}
(OUT/(sys.argv[1]+'-report.json')).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
