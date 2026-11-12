import argparse
import json
from .certificate import HERE,certificate,validate_pins
parser=argparse.ArgumentParser(description='Compute actual current/Gauss projection and finite Coulomb source choice')
mode=parser.add_mutually_exclusive_group(required=True)
mode.add_argument('--write',action='store_true');mode.add_argument('--check',action='store_true')
a=parser.parse_args();path=HERE/'results.json'
if a.check:saved=json.loads(path.read_text());validate_pins(saved)
value=certificate()
if a.write:path.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
elif value!=saved:raise ValueError('fresh native current/Coulomb receipt differs')
print(value['verdict'])
