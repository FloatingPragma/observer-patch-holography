from pathlib import Path
import argparse,json
from model.openwave_pure_su2_source.custody import pin_package,equivalent
from .source import run
HERE=Path(__file__).resolve().parent
def pins():return pin_package(HERE,{'model/openwave_m10_review_gate/contract.json'})
if __name__=='__main__':
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true')
    p.add_argument('--openwave-root',type=Path,required=True);a=p.parse_args()
    before=pins();v=run(a.openwave_root);assert pins()==before;v['source_sha256']=before
    out=HERE/'results.json'
    if a.write:out.write_text(json.dumps(v,sort_keys=True,indent=2)+'\n')
    else:assert equivalent(json.loads(out.read_text()),v)
    print(v['verdict'])
