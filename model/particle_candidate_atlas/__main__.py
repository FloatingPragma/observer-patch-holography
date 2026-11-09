"""Reproduce the candidate atlas and its exact comparison audit."""
import argparse
from . import builder

def main():
    p=argparse.ArgumentParser(description=__doc__)
    group=p.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true')
    group.add_argument('--check',action='store_true')
    args=p.parse_args()
    builder.write() if args.write else builder.check()
    print('PASS: particle candidate atlas and exact comparison audit')

if __name__=='__main__':main()
