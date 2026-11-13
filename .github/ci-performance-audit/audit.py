"""Authenticate the public-science baseline and reviewed CI overlay."""
from pathlib import Path
import hashlib, json, os, subprocess, sys

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
manifest_path=HERE/'manifest.json'
manifest=json.loads(manifest_path.read_text())
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
phase=sys.argv[1]
assert phase in ('before','after'), phase
output=Path(os.environ['AUDIT_OUTPUT']).resolve()
assert not output.is_relative_to(ROOT), 'audit outputs must be external'
output.mkdir(parents=True,exist_ok=True)
for path,expected in manifest['input_sha256'].items():
    assert sha(ROOT/path)==expected, ('immutable input changed',path)
receipt='model/particle_identification/completion_results.json'
report={'schema':'oph.ci_performance_public_baseline_phase.v1','phase':phase,
    'git_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
    'base_commit':manifest['base_commit'],'manifest_sha256':sha(manifest_path),
    'input_files':len(manifest['input_sha256']),'all_frozen_input_bytes_exact':True,
    'scientific_receipt_path':receipt,'scientific_receipt_sha256':sha(ROOT/receipt),
    'zero_new_scientific_results':True,'zero_Nature_validation_credit':True}
(output/(phase+'-source-report.json')).write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,sort_keys=True),flush=True)
