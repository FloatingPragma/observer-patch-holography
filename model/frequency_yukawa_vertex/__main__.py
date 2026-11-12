import argparse
import json
from .vertex import HERE,certificate

p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');p.add_argument('--check',action='store_true');a=p.parse_args()
r=certificate();text=json.dumps(r,indent=2,sort_keys=True,allow_nan=False)+'\n';path=HERE/'results.json'
if a.write:path.write_text(text)
if a.check and path.read_text()!=text:raise ArithmeticError('finite-frequency receipt changed')
print('Full-source finite-frequency vertex and external-leg amputation verified.')
