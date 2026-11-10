import argparse,json
from . import model
parser=argparse.ArgumentParser(description=__doc__)
group=parser.add_mutually_exclusive_group(required=True)
group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
args=parser.parse_args(); path=model.HERE/'results.json'
if args.check:
    saved=json.loads(path.read_text());model.filtered.higgs.validate(saved['input_sha256'])
    for name,sha in saved['owned_source_sha256'].items():
        if model.digest(model.HERE/name)!=sha:raise ValueError('owned current source differs: '+name)
actual=model.certificate()
if args.write:path.write_text(json.dumps(actual,indent=2,sort_keys=True)+'\n')
elif actual!=saved:raise ValueError('native Higgs current receipt differs')
print('PASS: same elementary/Higgs excitation, original quantum current and all finite records')
