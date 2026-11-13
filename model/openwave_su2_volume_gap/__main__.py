"""Original exact premise-map producer and independent replay."""
from argparse import ArgumentParser
from hashlib import sha256
import json
from pathlib import Path
from .source import run
import model.openwave_pure_su2_source.source as parent

HERE = Path(__file__).resolve().parent
ROOT = Path(parent.__file__).resolve().parents[2]


def pins():
    paths = ('model/openwave_pure_su2_source/source.py',
             'model/openwave_pure_su2_source/PROOF.md',
             'model/openwave_pure_su2_source/results.json',
             'model/openwave_su2_transfer/transfer.py',
             'model/openwave_su2_transfer/PROOF.md',
             'model/openwave_su2_transfer/results.json')
    result = {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}
    for rel in ('model/openwave_pure_su2_source/results.json', 'model/openwave_su2_transfer/results.json'):
        receipt = json.loads((ROOT/rel).read_text())
        for p, expected in receipt['source_sha256'].items():
            actual = sha256((ROOT/p).read_bytes()).hexdigest()
            assert actual == expected, ('stale parent', p)
            if p in result:
                assert result[p] == actual
            result[p] = actual
    for p in HERE.iterdir():
        if p.suffix in ('.py', '.md', '.tex'):
            result['model/openwave_su2_volume_gap/'+p.name] = sha256(p.read_bytes()).hexdigest()
    return result


def main():
    parser = ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run(); result['source_sha256'] = pins()
    text = json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n'
    output = HERE/'results.json'
    if args.write:
        output.write_text(text)
    else:
        assert output.read_text() == text, 'full exact producer mismatch'
    print(result['verdict'])

if __name__ == '__main__':
    main()
