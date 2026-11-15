"""Writer and cold replay against one authenticated upstream checkout."""
from pathlib import Path
import argparse,json
from hashlib import sha256
from model.openwave_pure_su2_source.custody import ROOT,pin_package,equivalent
from .source import run
HERE=Path(__file__).resolve().parent
def pins():
    parent='model/openwave_admission_packet/results.json'
    paths={'model/openwave_m10_review_gate/contract.json',parent}
    for p,h in json.loads((ROOT/parent).read_text())['source_sha256'].items():
        assert sha256((ROOT/p).read_bytes()).hexdigest()==h,('platform parent source changed',p)
        paths.add(p)
    return pin_package(HERE,paths)
if __name__=='__main__':
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true')
    p.add_argument('--openwave-root',type=Path,required=True);a=p.parse_args()
    before=pins();v=run(a.openwave_root);assert before==pins();v['source_sha256']=before
    path=HERE/'results.json'
    if a.write:path.write_text(json.dumps(v,sort_keys=True,indent=2,ensure_ascii=False)+'\n')
    else:assert equivalent(json.loads(path.read_text()),v),'complete replay differs'
    print(v['verdict'])
