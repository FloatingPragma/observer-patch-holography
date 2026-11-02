"""Write or independently recompute the current-spectrum proof packet."""
import argparse
import json
from . import packet


def main():
    parser=argparse.ArgumentParser()
    modes=parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--write',action='store_true')
    modes.add_argument('--check',action='store_true')
    args=parser.parse_args()
    if args.write:
        value=packet.build()
        packet.OUTPUT.write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n')
    else:
        from tools.strict_json import loads
        packet.validate(loads(packet.OUTPUT.read_text()))
    print('PASS spin current spectrum')


if __name__=='__main__':main()
