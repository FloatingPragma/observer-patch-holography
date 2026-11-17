import argparse,json
from .certificate import HERE,certificate,equivalent,inputs
try:
    from tools.pin_scope import bound_pins, scoped
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins, scoped
p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');args=p.parse_args()
path=HERE/'results.json'
if args.check:
    saved=json.loads(path.read_text())
    if bound_pins(saved['input_sha256'])!=scoped(inputs()):raise ValueError('precession source changed')
value=certificate()
if args.write:path.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
elif not equivalent(saved,value):raise ValueError('precession numerical receipt differs')
print(value['verdict'])
