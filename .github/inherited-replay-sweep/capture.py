"""Temporary immutable-baseline diagnostics; never changes acceptance rules."""
from pathlib import Path
import hashlib,importlib,json,os,subprocess,sys,time,traceback
ROOT=Path.cwd();sys.path.insert(0,str(ROOT))
from model.verify_particle_extensions import PACKAGES
OUT=Path(os.environ['AUDIT_OUTPUT']);OUT.mkdir(parents=True,exist_ok=True)
index=int(sys.argv[1]);names=PACKAGES[index::8];assert len(PACKAGES)==48
summary={'shard':index,'commit':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),'packages':[]}
def sha(b):return hashlib.sha256(b).hexdigest()
def changes(a,b,path=''):
    if type(a) is not type(b):return [{'path':path,'kind':'type','before':a,'after':b}]
    if isinstance(a,dict):
        out=[]
        for k in sorted(a.keys()|b.keys()):
            p=path+'/'+str(k).replace('~','~0').replace('/','~1')
            if k not in a or k not in b:out.append({'path':p,'kind':'key','before':a.get(k),'after':b.get(k)})
            else:out+=changes(a[k],b[k],p)
        return out
    if isinstance(a,list):
        if len(a)!=len(b):return [{'path':path,'kind':'length','before':len(a),'after':len(b)}]
        return [c for i,(x,y) in enumerate(zip(a,b)) for c in changes(x,y,path+'/'+str(i))]
    if a==b:return []
    d={'path':path,'kind':type(a).__name__,'before':a,'after':b}
    if isinstance(a,float):d['absolute_delta']=abs(a-b)
    return [d]
for name in names:
    out=OUT/name;out.mkdir();target=ROOT/'model'/name/'results.json';original=target.read_bytes()
    tracked=subprocess.check_output(['git','ls-files','-z','--','model/'+name]).decode().split('\0');saved={p:(ROOT/p).read_bytes() for p in tracked if p}
    report={'package':name,'before_sha256':sha(original)}
    (out/'before.json').write_bytes(original)
    def run(flag,label):
        argv=[sys.executable,'-m','model.'+name,flag]
        if name=='openwave_boosted_records':argv+=['--openwave-root',os.environ['OPH_OPENWAVE_ROOT'],'--python',os.environ['OPH_TAICHI_PYTHON']]
        start=time.monotonic()
        with (out/(label+'.log')).open('w')as log:
            try:code=subprocess.run(argv,stdout=log,stderr=subprocess.STDOUT,timeout=2400).returncode
            except subprocess.TimeoutExpired:code=124
        report[label+'_seconds']=time.monotonic()-start
        return code
    try:
        module=importlib.import_module('model.'+name+'.__main__');pins=module.pins()
        report['retained_source_pins_current']=json.loads(original)['source_sha256']==pins
        (out/'source-pins.json').write_text(json.dumps(pins,indent=2)+'\n')
        report['retained_cold_exit']=run('--check','retained-cold')
        if report['retained_cold_exit']!=0:
            report['actual_writer_exit']=run('--write','actual-writer')
            after=target.read_bytes();(out/'after.json').write_bytes(after);report['after_sha256']=sha(after)
            (out/'typed-diff.json').write_text(json.dumps(changes(json.loads(original),json.loads(after)),indent=2,allow_nan=False)+'\n')
            report['post_writer_source_pins_current']=json.loads(after)['source_sha256']==module.pins()
            if report['actual_writer_exit']==0:report['fresh_cold_exit']=run('--check','fresh-cold')
    except Exception as error:
        report['exception']=repr(error);(out/'exception.txt').write_text(traceback.format_exc())
    finally:
        for p,b in saved.items():(ROOT/p).write_bytes(b)
    clean=subprocess.run(['git','diff','--exit-code'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
    report['tracked_snapshot_restored']=clean.returncode==0
    if clean.returncode:(out/'unexpected-worktree-diff.txt').write_text(clean.stdout)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');summary['packages'].append(report)
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(report),flush=True)
summary['all_scheduled_packages_attempted']=len(summary['packages'])==len(names)
summary['retained_failures']=[r['package'] for r in summary['packages'] if r.get('retained_cold_exit')!=0]
summary['all_tracked_snapshots_restored']=all(r['tracked_snapshot_restored']for r in summary['packages'])
(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary),flush=True)
assert summary['all_scheduled_packages_attempted'] and summary['all_tracked_snapshots_restored']
