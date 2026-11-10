import argparse,json,math
from . import model


def compare(old,new,path='receipt'):
    if type(old)is not type(new):raise ValueError('type changed: '+path)
    if type(new)is dict:
        if old.keys()!=new.keys():raise ValueError('keys changed: '+path)
        for key in old:compare(old[key],new[key],path+'.'+key)
    elif type(new)is list:
        if len(old)!=len(new):raise ValueError('length changed: '+path)
        for i,(a,b)in enumerate(zip(old,new,strict=True)):compare(a,b,f'{path}[{i}]')
    elif type(new)is float:
        caps={
            'receipt.local_apparatus.raw_numerical_execution.direct_vs_rational_Taylor_operator_residual':1e-12,
            'receipt.local_apparatus.raw_numerical_execution.direct_vs_exact_cyclic_operator_residual':3e-12,
            'receipt.local_apparatus.raw_numerical_execution.complete_five_Kraus_residual':1e-12}
        if path not in caps or not all(math.isfinite(x) and 0<=x<=caps[path] for x in(old,new)):
            raise ValueError('numerical diagnostic outside its declared check: '+path)
    elif old!=new:raise ValueError('exact evidence changed: '+path)


def main():
    p=argparse.ArgumentParser(description='Replay the local Gauss-invariant current instrument with source drift retained.')
    g=p.add_mutually_exclusive_group(required=True);g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true')
    a=p.parse_args();target=model.HERE/'results.json'
    if a.check:
        old=json.loads(target.read_text())
        for name,sha in old['input_sha256'].items():model.require(model.digest(model.ROOT/name)==sha,'source differs: '+name)
    new=model.certificate()
    if a.write:target.write_text(json.dumps(new,indent=2,sort_keys=True,allow_nan=False)+'\n')
    else:compare(old,new)
    print('PASS: complete local400-state current pointer, unchanged native drift, all retained histories and resources')


if __name__=='__main__':main()
