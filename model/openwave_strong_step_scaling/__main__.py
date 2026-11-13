import argparse
import hashlib
import json
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from model.openwave_larmor_dynamics.certificate import equivalent
from .certificate import calculate

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PARENTS = ('openwave_common_hyperfine_spectrum', 'openwave_string_positive_cone')


def pins():
    files = {str(p.relative_to(ROOT)) for p in HERE.iterdir()
             if p.suffix in ('.py', '.md', '.tex') or p.name == 'inputs.json'}
    for parent in PARENTS:
        receipt = ROOT / 'model' / parent / 'results.json'
        files.add(str(receipt.relative_to(ROOT)))
        bindings = json.loads(receipt.read_text())['source_sha256']
        for path, digest in bindings.items():
            assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == digest, path
        files.update(bindings)
    files = import_closure(files)
    return {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in sorted(files)}


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    before = pins()
    result = calculate()
    assert before == pins()
    result['source_sha256'] = before
    path = HERE / 'results.json'
    if args.write:
        path.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    else:
        assert equivalent(json.loads(path.read_text()), result)
    print(result['verdict'])


if __name__ == '__main__':
    main()
