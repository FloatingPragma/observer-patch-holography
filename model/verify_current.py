"""Current replay routing with three explicitly enumerated custody successors.

The original verifier, artifact/source guards, group partition and every other
replay declaration retain their exact source and ordering.
"""
import argparse,hashlib,json,os,subprocess,sys
from pathlib import Path
from . import verify as retained
from tools import strict_json

BASE_VERIFY_SHA256='01a5965d14680094226e90955fa4f770bf8e65b76b5116662fce8bbb6fff3b31'
BASE_REPLAYS_SHA256='87c475d86192b569fe969485ee0bab67bf16d7bfe3251713b9ffb62a3ee1fdfa'
LEGACY=('model.protected_particle_readout',None,'model/protected_particle_readout/results.json')
CURRENT=('model.protected_readout_numerics',None,'model/protected_readout_numerics/results.json')
ROUTES={LEGACY:CURRENT, **{
    ('model.replay_native_fixed',name,'model/'+name+'/results.json'):
    ('model.completion.field_current_measurement',name,'model/'+name+'/results.json')
    for name in ('native_local_energy_clock','native_fair_ground_readout')}}


def digest(value):return hashlib.sha256(value).hexdigest()

def replay_digest(rows):return digest(json.dumps(tuple(rows),separators=(',',':')).encode())

def current_replays(rows=None):
    if digest((retained.ROOT/'model/verify.py').read_bytes())!=BASE_VERIFY_SHA256:
        raise ValueError('retained verifier source changed; explicit route review required')
    rows=retained.REPLAYS if rows is None else tuple(rows)
    if (replay_digest(rows)!=BASE_REPLAYS_SHA256
            or any(rows.count(old)!=1 or new in rows for old,new in ROUTES.items())):
        raise ValueError('unreviewed replay inventory/order or duplicate successor')
    return tuple(ROUTES.get(row,row) for row in rows)


def verify_successor_source_bindings():
    from .protected_readout_numerics import model
    from .protected_particle_readout import packet
    model.specification()
    value=strict_json.load(retained.ROOT/CURRENT[2]);pins=dict(packet.source_bindings())
    owned={p.name:model.digest(p)for p in sorted(model.HERE.iterdir())if p.suffix in('.py','.md')}
    pins.update({str((model.HERE/name).relative_to(retained.ROOT)):sha for name,sha in owned.items()})
    for name in('contract.json','retained.json'):pins[str((model.HERE/name).relative_to(retained.ROOT))]=model.digest(model.HERE/name)
    spec,_,_=model.specification();pins[spec['retained_array']]=spec['retained_array_sha256']
    if value.get('input_sha256')!=pins or value.get('owned_source_sha256')!=owned:
        raise ValueError('numerical successor source scope or custody changed')
    from .completion import field_current_measurement as measurement
    from .native_energy_measurement.packet import source_bindings
    measurement.authenticate(retained.ROOT/measurement.retained.PARENT,
                             measurement.retained.RETAINED_SHA256,source_bindings())
    for name in measurement.SOURCES:
        source=measurement.owner(name)
        receipt=strict_json.load(source.OUTPUT)
        if receipt['source_bindings']!=source.source_bindings():
            raise ValueError('current complete-measurement consumer source changed: '+name)


def verify_all_hashes():
    current_replays()
    count=retained.verify_artifacts()
    verify_successor_source_bindings()
    return count


def selected_groups(index=None,count=None):
    return retained.selected_replay_groups(index,count,replays=current_replays())


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--hashes-only',action='store_true')
    p.add_argument('--shard-index',type=int);p.add_argument('--shard-count',type=int)
    a=p.parse_args()
    try:groups=selected_groups(a.shard_index,a.shard_count)
    except ValueError as error:p.error(str(error))
    count=verify_all_hashes();print(f'Verified {count} retained artifacts and the numerical successor source',flush=True)
    if a.hashes_only:return
    environment=dict(os.environ)
    for name in('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):environment.setdefault(name,'1')
    for group in groups:
        if group==retained.EMPIRICAL_REPLAYS:
            subprocess.run([sys.executable,'-m','model.replay_empirical'],cwd=retained.ROOT,env=environment,check=True)
            continue
        (module,command,receipt),=group
        argv=[sys.executable,'-m',module]
        if command is not None:argv.append(command)
        argv.extend(['--check',receipt]);print(f'Replaying {receipt}',flush=True)
        subprocess.run(argv,cwd=retained.ROOT,env=environment,check=True)
    print(f'All {sum(map(len,groups))} selected current receipts reproduced',flush=True)


if __name__=='__main__':main()
