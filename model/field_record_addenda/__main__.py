"""Replay every additive public package in its declared dependency order."""
import argparse,os,subprocess,sys
from .registry import PACKAGES
from .provenance import ROOT


def main():
    parser=argparse.ArgumentParser();mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write',action='store_true');mode.add_argument('--check',action='store_true')
    args=parser.parse_args();command='--write'if args.write else'--check'
    env=dict(os.environ)
    for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
        env.setdefault(name,'1')
    for package in PACKAGES:
        subprocess.run([sys.executable,'-B','-m','model.field_record_addenda.'+package,command],
                       cwd=ROOT,env=env,check=True)
    print('PASS: all eleven complete additive producers, each in a cold process')


if __name__=='__main__':main()
