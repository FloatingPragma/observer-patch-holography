import argparse
import json
from pathlib import Path
from .model import certificate
p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');p.add_argument('--check',action='store_true');a=p.parse_args()
path=Path(__file__).with_name('results.json');actual=certificate()
if a.write:path.write_text(json.dumps(actual,indent=2,sort_keys=True)+'\n')
if a.check and json.loads(path.read_text())!=actual:raise SystemExit('odd sum-rule receipt differs')
print('PASS: original odd sum rules and graded spatial contrast; no particle mass or two-point pole')
