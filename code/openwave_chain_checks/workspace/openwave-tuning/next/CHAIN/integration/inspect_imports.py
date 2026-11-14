"""Maintainer discovery of local import/read dependencies, without acceptance."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[name]='1'
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,importlib.util
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[3];CHAIN=HERE.parent
seen=set()

def audit(event,args):
    if event=='open' and isinstance(args[0],(str,bytes)):
        p=Path(os.fsdecode(args[0])).absolute()
        if p.is_relative_to(ROOT) and not p.is_relative_to(HERE) and '__pycache__' not in p.parts and p.is_file():seen.add(p)
sys.addaudithook(audit)

def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m

print('Discovering quantum constructor dependencies',flush=True)
p=load('closure_particle',CHAIN/'particle/particle_chain.py');c=p.ParticleChain();c.clock.preparation()
load('closure_dynamic',CHAIN/'dynamical_field/dynamical_field.py')
print('Discovering classical helper dependencies',flush=True)
b=ROOT/'openwave-tuning/auditor/native_field_bridge';f=load('field_feedback',b/'field_feedback.py');f.graph();f.literal_native()
load('closure_reciprocal',b/'reciprocal_feedback.py');load('closure_field',CHAIN/'field_limit/controlled_field.py');load('closure_force',CHAIN/'force/force_measurement.py')
for module in list(sys.modules.values()):
    name=getattr(module,'__file__',None)
    if name:
        path=Path(name).resolve()
        if path.is_relative_to(ROOT) and not path.is_relative_to(HERE) and path.is_file():seen.add(path)
files={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(seen)}
(HERE/'IMPORT_CLOSURE.json').write_text(json.dumps({'schema':'oph.chain.local_runtime_discovery.v1','scope':'Dependency discovery only; no acceptance claim or mutable-root execution used as final evidence.','files':files},indent=2)+'\n')
print(len(files),'local files discovered',flush=True)
