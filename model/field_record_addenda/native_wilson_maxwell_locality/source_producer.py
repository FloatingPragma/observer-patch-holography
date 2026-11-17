"""Source-pinned same-cube field/current propagation certificate."""
import argparse,hashlib,json
from pathlib import Path
from . import locality
from ..provenance import import_closure
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[3]))
    from tools.pin_scope import bound_pins
HERE=Path(__file__).resolve().parent
RER=locality.RER


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def certificate():
    value=locality.scientific();parent=locality.admit()
    seeds=set(bound_pins(parent['input_sha256']))|{locality.PARENT,
        'model/native_wilson_maxwell/field.py','model/native_wilson_maxwell/current.py',
        'model/native_wilson_maxwell_refinement/bounds.py',
        'model/particle_common_contract/provenance.py'}
    value['input_sha256']=bound_pins({p:digest(RER/p) for p in sorted(import_closure(seeds))})
    value['owned_source_sha256']=bound_pins({p.name:digest(p) for p in sorted(HERE.iterdir()) if p.suffix in ('.py','.md')})
    return value


def main():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true')
    args=p.parse_args();value=certificate();path=HERE/'results.json'
    if args.write:path.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
    elif value!=json.loads(path.read_text()):raise ValueError('locality certificate differs')
    print('PASS: same-source fixed-cube gauge-invariant propagation and complete records')


if __name__=='__main__':main()
