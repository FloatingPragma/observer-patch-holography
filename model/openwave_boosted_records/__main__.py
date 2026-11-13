"""Export and replay the exact common-source reducing-sector kernel."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from .source import source_packet,boost_checks,theorem_metadata
from model.particle_common_contract.provenance import import_closure
from model.openwave_larmor_dynamics.certificate import equivalent

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
UPSTREAM='ef55175d04737e73c92f33f9f9cf3d17c593da54'


def pins():
    parent='model/openwave_common_electromagnetism/results.json'
    seeds=set(json.loads((ROOT/parent).read_text())['source_sha256'])|{parent}
    seeds|={str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py','.md','.tex') or p.name=='source_packet.json'}
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(import_closure(seeds))}


def run(openwave,python):
    revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=openwave,text=True).strip()
    assert revision==UPSTREAM
    assert not subprocess.check_output(['git','diff','--name-only','HEAD','--'],cwd=openwave).strip()
    assert json.loads((HERE/'source_packet.json').read_text())==source_packet()
    overlay=openwave/'openwave/xperiments/oph_common_boost';overlay.mkdir(exist_ok=True)
    assert all(p.name in ('__init__.py','engine.py','__pycache__') for p in overlay.iterdir())
    for name in ('__init__.py','engine.py'):
        target=overlay/name
        if target.exists() and target.read_bytes()!=(HERE/name).read_bytes():
            previous=json.loads((HERE/'results.json').read_text())['source_sha256']
            assert hashlib.sha256(target.read_bytes()).hexdigest()==previous[str((HERE/name).relative_to(ROOT))], 'unrecognized existing overlay'
        shutil.copy2(HERE/name,target)
    before=pins()
    with tempfile.TemporaryDirectory(prefix='oph-tenth-boost-') as temp:
        out=Path(temp)/'records.json'
        subprocess.run([python,'-m','openwave.xperiments.oph_common_boost.engine',
            '--source',str(HERE/'source_packet.json'),'--output',str(out)],cwd=openwave,
            env={**os.environ,'PYTHONPATH':str(openwave),'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'},check=True)
        records=json.loads(out.read_text())
    assert pins()==before
    return {'schema':'oph.common_source_boosted_records.v1',
        'verdict':'PASS_COMPLETE_REDUCING_SECTOR_RECORDS_AND_CONTROLLED_BOOST_DEFECT',
        'upstream_commit':revision,'upstream_tracked_source_modified':False,
        'official_admission_or_score_award':False,'source_sha256':before,
        'source_domain':source_packet()['source_domain'],
        'theorem':theorem_metadata(),'boosted_energy_and_velocity_checks':boost_checks(),
        'portable_Taichi_execution':records}


def main():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true')
    p.add_argument('--openwave-root',type=Path,default=os.environ.get('OPH_OPENWAVE_ROOT'))
    p.add_argument('--python',default=os.environ.get('OPH_TAICHI_PYTHON'))
    args=p.parse_args()
    if args.openwave_root is None or args.python is None:
        p.error('provide --openwave-root/--python or OPH_OPENWAVE_ROOT/OPH_TAICHI_PYTHON; a genuine pinned OpenWave plus Taichi execution is required')
    packet=HERE/'source_packet.json'
    if args.write:packet.write_text(json.dumps(source_packet(),sort_keys=True,indent=2)+'\n')
    result=run(args.openwave_root.resolve(),args.python);path=HERE/'results.json'
    if args.write:path.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    else:assert equivalent(json.loads(path.read_text()),result)
    print(result['verdict'])


if __name__=='__main__':main()
