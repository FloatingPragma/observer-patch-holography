"""Read-only strict-verification descendant; v5 science and records untouched.

The original v5 comparator accepted JSON true as integer1. This descendant
requires exact JSON scalar types while replaying the same frozen experiment.
There is deliberately no record-creation or baseline-regeneration option.
"""
import argparse
import hashlib
import json
import math
import time
from pathlib import Path
import confirm

ROOT=Path(__file__).resolve().parent


def require(condition,message):
    if not condition:raise RuntimeError(message)


def compare(a,b,path='$'):
    require(type(a)is type(b),path+' strict JSON type mismatch')
    if isinstance(a,dict):
        require(a.keys()==b.keys(),path+' object keys mismatch')
        for key in a:compare(a[key],b[key],path+'.'+key)
    elif isinstance(a,list):
        require(len(a)==len(b),path+' array length mismatch')
        for i,(x,y)in enumerate(zip(a,b)):compare(x,y,path+'['+str(i)+']')
    elif isinstance(a,float):
        require(math.isfinite(a)and math.isfinite(b),path+' nonfinite number')
        if path.endswith('.selected_clock_probability_given_entry')or path.endswith('.selected_clock_probability_including_calibration_given_entry'):
            require((a==0)==(b==0)and(a==b or abs(a-b)<=1e-9*max(abs(a),abs(b))),path+' rare-cost relative/zero mismatch')
        else:require(abs(a-b)<=2e-10+2e-12*max(abs(a),abs(b)),path+' float roundoff bound exceeded')
    else:require(a==b,path+' exact scalar mismatch')


def main():
    argparse.ArgumentParser(description=__doc__).parse_args()
    binding=json.loads((ROOT/'VERIFIER_V6_BINDING.json').read_text())
    require(hashlib.sha256(Path(__file__).read_bytes()).hexdigest()==binding['verifier_sha256'],'Verifier bytes changed')
    folder=ROOT/'CONFIRMATION'
    require(hashlib.sha256((folder/'INDEX.json').read_bytes()).hexdigest()==binding['source_index_sha256'],'Original experiment index changed')
    confirm.verify_freeze();confirm.verify_custody(folder)
    confirm.compare=compare
    start=time.monotonic();result,artifacts=confirm.run(folder,False)
    compare({'candidate':result['candidate'],'files':artifacts},json.loads((folder/'INDEX.json').read_text()))
    print(json.dumps({'verified':True,'source_candidate':result['candidate'],
        'cases':len(result['rows']),'scientific_gates_pass':result['gates_pass'],
        'verification_revision':'v6 strict JSON scalar types; all tiny analytic path semantics retained',
        'source_index_sha256':binding['source_index_sha256'],'seconds':time.monotonic()-start},indent=2))


if __name__=='__main__':main()
