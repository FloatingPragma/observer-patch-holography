"""Generate raw loop records; this command performs the Monte Carlo itself."""
import argparse,json
from pathlib import Path
from .coupling import generate


def main():
    p=argparse.ArgumentParser();p.add_argument('--side',type=int,default=6)
    p.add_argument('--beta',type=float,default=6.);p.add_argument('--seed',type=int,required=True)
    p.add_argument('--start',choices=['cold','hot'],default='cold')
    p.add_argument('--samples',type=int,default=400);p.add_argument('--warmup',type=int,default=200)
    args=p.parse_args();path=Path(__file__).parent/'data'/f'N{args.side}_b{args.beta:.1f}_{args.start}_{args.seed}.npz'
    print(json.dumps(generate(path,**vars(args)),indent=2))


if __name__=='__main__':main()
