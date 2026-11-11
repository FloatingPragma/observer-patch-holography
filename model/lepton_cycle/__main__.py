import argparse
import json
from pathlib import Path
from tools import strict_json
from .certificate import report


def main():
    parser=argparse.ArgumentParser(description='Complete native lepton-cycle certificate')
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
    args=parser.parse_args();value=report();path=Path(__file__).with_name('results.json')
    if args.check and strict_json.load(path)!=value:
        raise ValueError('complete lepton-cycle receipt differs from replay')
    if args.write:path.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'checks':len(value['checks']),
        'rounds_for_1e_minus_6':value['rounds_for_ideal_trace_distance_1e_minus_6'],
        'ratios':[value['sensitivity']['muon_over_electron'],value['sensitivity']['tau_over_electron']]},indent=2))


if __name__=='__main__':main()
