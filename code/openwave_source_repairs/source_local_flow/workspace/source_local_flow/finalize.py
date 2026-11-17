"""Write a complete result before returning a fail-closed process status."""
import argparse,json
from pathlib import Path

def finish(path,result):
    gates=result.get('gates')
    if not isinstance(gates,dict)or not gates:raise ValueError('nonempty named gates required')
    result['gates']={name:bool(value)for name,value in gates.items()}
    result['all_gates_pass']=all(result['gates'].values())
    Path(path).write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
    return 0 if result['all_gates_pass']else 1

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('input');parser.add_argument('output')
    parser.add_argument('--reject-gate');args=parser.parse_args()
    result=json.loads(Path(args.input).read_text())
    if args.reject_gate:
        if args.reject_gate not in result['gates']:raise ValueError('unknown gate')
        result['gates'][args.reject_gate]=False
        result['review_mutation_only']=args.reject_gate
    raise SystemExit(finish(args.output,result))
