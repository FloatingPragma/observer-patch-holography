import argparse
import json
from pathlib import Path
from .certificate import certificate
from model.affine_gauss_charged_clock.certificate import validate_pins

p=argparse.ArgumentParser(description='Replay complete post-clock charged-field records')
g=p.add_mutually_exclusive_group(required=True);g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true')
a=p.parse_args();path=Path(__file__).with_name('results.json')
if a.check:
    saved=json.loads(path.read_text());validate_pins(saved['input_sha256'])
value=certificate()
if a.write:path.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
elif saved!=value:raise ValueError('post-clock record certificate differs')
print('PASS: complete26 post-clock outcomes with preparation, all phase labels and ground-fidelity error')
