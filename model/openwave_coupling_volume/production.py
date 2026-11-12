"""Independent larger-box production using the sealed native colour kernel."""
from pathlib import Path
import argparse,json
from model.openwave_scale_coupling.coupling import generate

def main():
 p=argparse.ArgumentParser();p.add_argument('--side',type=int,required=True);p.add_argument('--seed',type=int,required=True);p.add_argument('--start',choices=['cold','hot'],required=True);p.add_argument('--warmup',type=int,default=400);p.add_argument('--samples',type=int,default=400)
 a=p.parse_args();path=Path(__file__).parent/'data'/f'N{a.side}_b6.0_{a.start}_{a.seed}.npz'
 print(json.dumps(generate(path,**vars(a)),indent=2))
if __name__=='__main__':main()
