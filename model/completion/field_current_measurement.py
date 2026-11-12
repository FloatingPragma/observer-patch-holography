"""Current custody adapter for two unchanged complete-measurement consumers.

The two source files and their scientific functions remain immutable. Each is
loaded in a private module namespace; only its named parent-custody service is
replaced. The new service authenticates the actually replayed complete 66-word
receipt and delegates the exact payload/source-scope comparison to the sealed
adapter. No spectral, instrument, clock or numerical calculation is replaced.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
from tools import strict_json
from . import native_measurement_parent_adapter as retained
ROOT=Path(__file__).resolve().parents[2]
REPLAY='model/completion/field_current_measurement_replay.json'
REPLAY_SHA256='42d77ef93b0c4e62b98c1d3f7047db5c6d0ecc778a3a55deecea51fa628a9cfd'
CURRENT_SHA256='09263ef294f95ae60cf912c806708b947bee59d328cb2ed62d8764c2683c645b'
OWNED=('model/completion/field_current_measurement.py','model/completion/test_field_current_measurement.py',REPLAY)
SOURCE_FILES=(*retained.SOURCE_FILES,*OWNED)
SOURCES={'native_local_energy_clock': ('model/native_local_energy_clock/packet.py', '1bb8cda89d1b6f37be71578bed9cdfe7fd82c9ceeb71841f2c553523c4e91aa1'), 'native_fair_ground_readout': ('model/native_fair_ground_readout/packet.py', '147ba199056d7e9b7743a95c78d293adc1123b96006228af294a28427a85cfb1')}


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def authenticate(path,retained_sha256,live_sources,root=ROOT):
    root=Path(root).resolve()
    if Path(path).resolve()!=root/retained.PARENT or retained_sha256!=retained.RETAINED_SHA256:
        raise ValueError('unreviewed complete measurement parent')
    if sha(root/REPLAY)!=REPLAY_SHA256 or sha(root/retained.PARENT)!=CURRENT_SHA256:
        raise ValueError('actual current complete-measurement replay changed')
    if sha(root/retained.ARCHIVE)!=retained.RETAINED_SHA256:
        raise ValueError('retained complete-measurement authority changed')
    replay=strict_json.load(root/REPLAY)
    if (replay['fresh_sha256']!=CURRENT_SHA256 or replay['actual_binary64_replay_byte_exact'] is not True
            or replay['producer']!='model.native_energy_measurement.packet:produce'):
        raise ValueError('complete actual owner execution was not certified')
    retained.compare_payloads(strict_json.loads((root/retained.ARCHIVE).read_text()),
                              strict_json.loads((root/retained.PARENT).read_text()),live_sources)
    for p,h in live_sources.items():
        q=(root/p).resolve()
        if not q.is_relative_to(root) or not q.is_file() or sha(q)!=h:
            raise ValueError('current complete-measurement source changed: '+p)
    return CURRENT_SHA256


def service():
    values={name:getattr(retained,name)for name in dir(retained)if not name.startswith('__')}
    values.update(authenticate=authenticate,SOURCE_FILES=SOURCE_FILES,REPLAYED_SHA256=CURRENT_SHA256)
    return SimpleNamespace(**values)


def owner(name):
    if name not in SOURCES:raise ValueError('unreviewed current measurement consumer')
    path,h=SOURCES[name]
    if sha(ROOT/path)!=h:raise ValueError('unchanged scientific consumer source changed')
    module_name='model.'+name+'._field_current_owner'
    spec=importlib.util.spec_from_file_location(module_name,ROOT/path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.parent_custody=service()
    return module


def produce_local():return owner('native_local_energy_clock').produce()


def produce_ground():return owner('native_fair_ground_readout').produce()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('owner',choices=tuple(SOURCES))
    parser.add_argument('--check',type=Path,required=True)
    args=parser.parse_args();module=owner(args.owner);fresh=module.produce()
    actual=strict_json.load(args.check)
    if retained.canonical(actual)!=retained.canonical(fresh):
        raise ValueError('current consumer differs from complete original scientific replay')
    print('Verified unchanged complete-measurement consumer with explicit current custody: '+args.owner)


if __name__=='__main__':main()
