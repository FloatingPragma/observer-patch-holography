"""Recompute the exact finite abundance receipt."""
import argparse,json
from pathlib import Path
from .certificate import certificate

parser=argparse.ArgumentParser(description=__doc__)
mode=parser.add_mutually_exclusive_group(required=True)
mode.add_argument('--write',action='store_true');mode.add_argument('--check',action='store_true')
args=parser.parse_args();result=certificate()
failed=[key for key,value in result['checks'].items() if not value]
if failed:raise SystemExit('failed checks: '+', '.join(failed))
path=Path(__file__).with_name('results.json')
if args.write:path.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
elif json.loads(path.read_text())!=result:raise SystemExit('receipt differs from recomputation')
print(f"PASS: {len(result['checks'])} exact abundance checks; persistent interval {result['persistent_occupied_face_share']}")
