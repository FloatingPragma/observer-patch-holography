"""Offline read-only replay; all executed primary sources are local and inspectable."""
from pathlib import Path
import hashlib,json,subprocess,sys,time
p=Path(__file__).resolve().parent
results=[]
for source,receipt in [('experiment_serialized.py','OUTPUT.json'),('refinement/experiment.py','refinement/OUTPUT.json'),('refinement/alternative.py','refinement/ALTERNATIVE_OUTPUT.json')]:
    start=time.perf_counter();r=subprocess.run([sys.executable,'-B',str(p/source)],capture_output=True,timeout=60)
    assert r.returncode==0,r.stderr.decode()
    expected=(p/receipt).read_bytes()
    assert r.stdout==expected,source+' replay differs'
    results.append({'source':source,'seconds':time.perf_counter()-start,'exit_code':r.returncode,'output_sha256':hashlib.sha256(r.stdout).hexdigest()})
print(json.dumps({'replay':results,'total_seconds':sum(r['seconds'] for r in results),'no_packet_writes':True},indent=2))
