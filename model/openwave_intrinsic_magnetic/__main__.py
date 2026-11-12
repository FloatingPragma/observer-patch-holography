import argparse
import json
from .certificate import HERE,certificate,validate_pins
parser=argparse.ArgumentParser(description='Recompute the complete finite stationary magnetic source')
mode=parser.add_mutually_exclusive_group(required=True)
mode.add_argument('--write',action='store_true');mode.add_argument('--check',action='store_true')
a=parser.parse_args();path=HERE/'results.json'
if a.check:saved=json.loads(path.read_text());validate_pins(saved)
result=certificate()
if a.write:path.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
elif result!=saved:raise ValueError('fresh stationary magnetic receipt differs')
print(result['verdict'])
