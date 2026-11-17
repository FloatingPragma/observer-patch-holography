"""Eight falsifiable verifier controls, each retaining the registered fixture."""
import copy
import gzip
import hashlib
import json
from pathlib import Path
import tempfile
from verify_strict import check

HERE = Path(__file__).resolve().parent


def main():
    base_receipt = json.loads((HERE / 'RECEIPT.json').read_text())
    base_data = json.loads(gzip.decompress((HERE / 'EVIDENCE.json.gz').read_bytes()))
    cases = ('branch_coin', 'branch_current', 'missing_history', 'wrong_seed',
             'wrong_summary', 'false_award', 'probe_current', 'timeout_as_success')
    outcomes = {}
    with tempfile.TemporaryDirectory(prefix='row01-mutants-') as temp:
        out = Path(temp)
        for name in cases:
            receipt, data = copy.deepcopy(base_receipt), copy.deepcopy(base_data)
            if name == 'branch_coin':
                data['sectors'][0]['branches'][1][2] = 2
            elif name == 'branch_current':
                data['sectors'][2]['branches'][0][3] += 1
            elif name == 'missing_history':
                data['histories'].pop()
            elif name == 'wrong_seed':
                data['histories'][0]['seed'] += 1
            elif name == 'wrong_summary':
                receipt['histories']['native_balanced_within_horizon'] = 49
            elif name == 'false_award':
                receipt['full_row_solved'] = True
            elif name == 'probe_current':
                data['histories'][0]['probe'][0][3] += 1
            elif name == 'timeout_as_success':
                data['forced_timeout']['probe_outcome'] = 'conclusive'
            raw = gzip.compress((json.dumps(data, sort_keys=True,
                separators=(',', ':')) + '\n').encode(), mtime=0)
            receipt['evidence_sha256'] = hashlib.sha256(raw).hexdigest()
            (out / 'EVIDENCE.json.gz').write_bytes(raw)
            (out / 'RECEIPT.json').write_text(json.dumps(receipt))
            try:
                check(out)
            except (AssertionError, ValueError, IndexError, KeyError):
                outcomes[name] = 'rejected'
            else:
                outcomes[name] = 'MISSED'
    print(json.dumps({'cases': outcomes, 'passed': all(v == 'rejected' for v in
                       outcomes.values())}, indent=2))
    raise SystemExit(0 if all(v == 'rejected' for v in outcomes.values()) else 1)


if __name__ == '__main__':
    main()
