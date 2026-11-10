"""Replay the existing edge theorem with the computed source premise."""
import argparse
import hashlib
import json
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from .model import ROOT,SOURCE,certificate_data


def certificate():
    seeds={str(p.relative_to(ROOT)) for p in Path(__file__).parent.iterdir() if p.suffix in ('.py','.md')}
    seeds.update((str(SOURCE.relative_to(ROOT)),
                  'model/coherent_defect_wave/chiral_joint/two_cell/READOUT_PROOF.md',
                  'model/spin_shell_readout/THRESHOLD_READOUT.md','model/requirements.txt'))
    source=json.loads(SOURCE.read_text())
    seeds.update(source['input_sha256'])
    hashes={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(import_closure(seeds))}
    return {'schema':'oph.coherent_defect_wave.chiral_joint.spectral_edges.v1',
            **certificate_data(),'input_sha256':hashes}


def main():
    parser=argparse.ArgumentParser();group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
    args=parser.parse_args();path=Path(__file__).with_name('results.json')
    rendered=json.dumps(certificate(),indent=2,sort_keys=True)+'\n'
    if args.write:path.write_text(rendered)
    elif path.read_text()!=rendered:raise SystemExit('FAIL stale charged spectral edge attachment')
    print('PASS original-weight finite source edge and reused non-atomic criterion')


if __name__=='__main__':main()
