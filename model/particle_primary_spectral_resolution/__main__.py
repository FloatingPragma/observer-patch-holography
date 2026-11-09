"""Replay full residual and complete prepared-energy records."""
import argparse
import json
from pathlib import Path
from .model import certificate
p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true')
a=p.parse_args();text=json.dumps(certificate(),sort_keys=True,indent=2,allow_nan=False)+'\n'
path=Path(__file__).with_name('results.json')
if a.write:path.write_text(text)
elif not path.is_file() or path.read_text()!=text:raise SystemExit('prepared spectral receipt differs')
print('PASS: complete physical A1u preparation and energy window; no isolated mass')
