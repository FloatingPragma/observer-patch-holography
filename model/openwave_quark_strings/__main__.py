import argparse,hashlib,json
from pathlib import Path
from .strings import certificate
from model.particle_common_contract.provenance import import_closure
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]


def produce():
    result=certificate()
    paths=set(HERE.glob('*.py'))|set(HERE.glob('*.md'))
    # The import closure records the exact native algebra used by the producer.
    closure=import_closure([str((HERE/'strings.py').relative_to(ROOT))])
    paths.update(ROOT/p for p in closure)
    result['source_sha256']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}
    return result


def main():
    ap=argparse.ArgumentParser();g=ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true')
    args=ap.parse_args();r=produce();p=HERE/'results.json'
    if args.write:p.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n')
    elif json.loads(p.read_text())!=r:raise SystemExit('native quark string source result differs')
    print(r['verdict'])


if __name__=='__main__':main()
