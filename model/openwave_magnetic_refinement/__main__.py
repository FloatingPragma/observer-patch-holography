import argparse,json
from .certificate import HERE,certificate,inputs,equivalent
p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');args=p.parse_args()
path=HERE/'results.json'
if args.check:
    saved=json.loads(path.read_text())
    if saved['input_sha256']!=inputs():raise ValueError('magnetic witness or source changed')
value=certificate()
if args.write:path.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
elif not equivalent(saved,value):raise ValueError('magnetic refinement receipt differs')
print(value['verdict'])
