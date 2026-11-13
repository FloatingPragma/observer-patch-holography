from pathlib import Path
import hashlib,importlib,json,os,subprocess,sys
r=Path.cwd();sys.path.insert(0,str(r));name=sys.argv[1];out=Path(os.environ['AUDIT_OUTPUT']);out.mkdir(parents=True,exist_ok=True);target=r/'model'/name/'results.json';original=target.read_bytes();module=importlib.import_module('model.'+name+'.__main__');pins=module.pins();assert json.loads(original)['source_sha256']==pins
(out/'before.json').write_bytes(original);(out/'source-pins.json').write_text(json.dumps(pins,indent=2)+'\n')
report={'package':name,'commit':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),'source_pin_count':len(pins),'before_sha256':hashlib.sha256(original).hexdigest()}
def run(flag,label):
 with (out/(label+'.log')).open('w')as log:return subprocess.run([sys.executable,'-m','model.'+name,flag],stdout=log,stderr=subprocess.STDOUT).returncode
try:
 report['retained_cold_exit']=run('--check','retained-cold')
 report['actual_writer_exit']=run('--write','actual-writer');assert report['actual_writer_exit']==0
 after=target.read_bytes();(out/'after.json').write_bytes(after);report['after_sha256']=hashlib.sha256(after).hexdigest()
 assert json.loads(after)['source_sha256']==pins==module.pins()
 report['fresh_cold_exit']=run('--check','fresh-cold');assert report['fresh_cold_exit']==0
finally:target.write_bytes(original)
assert module.pins()==pins
subprocess.run(['git','diff','--exit-code'],check=True)
(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
