"""Regenerate and execute the portable full-native OpenWave solver."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from .export import source_packet
from model.particle_common_contract.provenance import import_closure
from model.common_source.chiral_candidate.action import SOURCE_FILES
from model.openwave_larmor_dynamics.certificate import equivalent

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
UPSTREAM='ef55175d04737e73c92f33f9f9cf3d17c593da54'


def pins():
    seeds=set(SOURCE_FILES)
    seeds.update(str(p.relative_to(ROOT)) for p in HERE.iterdir()
                 if p.suffix in ('.py','.md','.tex') or p.name=='source_packet.json')
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
            for p in sorted(import_closure(seeds))}


def run(openwave, python):
    actual=subprocess.check_output(['git','rev-parse','HEAD'],cwd=openwave,text=True).strip()
    assert actual==UPSTREAM, ('upstream revision',actual)
    assert not subprocess.check_output(['git','diff','--name-only','HEAD','--'],cwd=openwave).strip(), 'modified upstream tracked source'
    assert json.loads((HERE/'source_packet.json').read_text())==source_packet()
    overlay=openwave/'openwave/xperiments/oph_observer'
    overlay.mkdir(exist_ok=True)
    allowed={'__init__.py','engine.py','__pycache__'}
    assert all(p.name in allowed for p in overlay.iterdir())
    for name in ('__init__.py','engine.py'):
        target=overlay/name
        if target.exists():
            assert target.read_bytes()==(HERE/name).read_bytes(), 'different existing overlay'
        else:
            shutil.copy2(HERE/name,target)
    before=pins()
    with tempfile.TemporaryDirectory(prefix='oph-native-platform-') as temp:
        out=Path(temp)/'records.json'
        subprocess.run([python,'-m','openwave.xperiments.oph_observer.engine',
            '--source',str(HERE/'source_packet.json'),'--output',str(out)],cwd=openwave,
            env={**os.environ,'PYTHONPATH':str(openwave)},check=True)
        records=json.loads(out.read_text())
    assert pins()==before
    assert records['OpenWave_version']=='26.12'
    cases=[c for row in records['runs'] for c in row['cases']]
    return {'schema':'oph.openwave_platform_kernel.v1',
        'verdict':'PASS_FULL_NATIVE_HEADLESS_PLATFORM_EVOLUTION',
        'source_sha256':before,'upstream_commit':actual,
        'upstream_modified_tracked_files':False,
        'official_admission_or_award':False,
        'full_native_execution':records,
        'maximum_twenty_step_record_TV':max(c['complete_detector_TV'] for c in cases if c['steps']==20),
        'independent_native_derivative_error':max(r['full_carrier_derivative_error'] for r in records['runs']),
        'negative_control_detected':all(r['kinetic_sign_mutation_difference']>.1 for r in records['runs']),
        'source_scope':'Exact exported native comparison member, all90CAR/180Nambu coordinates; no physical continuum or common interacting-spectrum transfer.'}


def main():
    p=argparse.ArgumentParser()
    g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true')
    p.add_argument('--openwave-root',type=Path,required=True)
    p.add_argument('--python',default=sys.executable)
    a=p.parse_args();result=run(a.openwave_root.resolve(),a.python);path=HERE/'results.json'
    if a.write:path.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    else:assert equivalent(json.loads(path.read_text()),result)
    print(result['verdict'],result['maximum_twenty_step_record_TV'])


if __name__=='__main__':main()
