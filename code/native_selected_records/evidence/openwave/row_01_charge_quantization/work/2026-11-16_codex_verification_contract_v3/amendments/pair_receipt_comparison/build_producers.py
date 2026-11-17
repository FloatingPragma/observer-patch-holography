"""Construct repaired direct producers without changing any scientific calculation."""
import hashlib
import json
from pathlib import Path
import shutil
from build_descendants import unconditional, HERE, INTEGRATION

BASE=INTEGRATION/'source_export_v2/source/evidence/openwave/row_01_charge_quantization/work'
DEST=HERE/'descendants/producers'
DEST.mkdir()
changes=[]


def copy_code(old, new):
    text,count=unconditional(old.read_text());new.write_text(text)
    changes.append({'original':str(old.relative_to(INTEGRATION)),'descendant':str(new.relative_to(HERE)),
                    'original_sha256':hashlib.sha256(old.read_bytes()).hexdigest(),
                    'unconditional_checks':count})


old=BASE/'2026-11-15_codex_pr723_vacuum_reservoir';new=DEST/'reservoir';new.mkdir()
for name in ('experiment.py','model.py','native_function.py','reader.py'):
    copy_code(old/name,new/name)
for name in ('TASK.json','seams.json'):shutil.copyfile(old/name,new/name)
(new/'saved').mkdir()
for name in ('RECEIPT.json','EVIDENCE.json.gz'):shutil.copyfile(old/name,new/'saved'/name)
text=(new/'experiment.py').read_text().split("if __name__ == '__main__':")[0]
text+='''if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--verify', type=Path, default=HERE/'saved')
    mode.add_argument('--record-to', type=Path)
    args = parser.parse_args()
    if args.record_to is not None:
        _require(not args.record_to.exists(), 'First-record directory already exists')
    receipt, evidence = calculate()
    payload = gzip.compress(encode(evidence), mtime=0)
    receipt['evidence_sha256'] = hashlib.sha256(payload).hexdigest()
    if args.record_to is not None:
        args.record_to.mkdir()  # Exclusive creation; failed scientific records are retained.
        (args.record_to/'EVIDENCE.json.gz').write_bytes(payload)
        (args.record_to/'RECEIPT.json').write_bytes(encode(receipt))
    else:
        _require(payload == (args.verify/'EVIDENCE.json.gz').read_bytes(), 'Saved full evidence differs')
        _require(encode(receipt) == (args.verify/'RECEIPT.json').read_bytes(), 'Saved receipt differs')
    print(json.dumps(receipt, indent=2))
    raise SystemExit(0 if receipt['passed'] and all(receipt['checks'].values()) else 1)
'''
(new/'experiment.py').write_text(text)

old=BASE/'2026-11-16_codex_common_refinement_family';new=DEST/'pair';new.mkdir()
for name in ('core.py','pair_source.py'):copy_code(old/name,new/name)
(new/'saved').mkdir()
for name in ('PAIR_SOURCE.json.gz','PAIR_RECEIPT.json'):shutil.copyfile(old/name,new/'saved'/name)
text=(new/'pair_source.py').read_text().split("if __name__=='__main__':")[0]
text+='''if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--verify', type=Path, default=Path(__file__).resolve().parent/'saved')
    mode.add_argument('--record-to', type=Path)
    args = parser.parse_args()
    if args.record_to is not None:
        _require(not args.record_to.exists(), 'First-record directory already exists')
    result = run()
    payload = gzip.compress((json.dumps(result,sort_keys=True,separators=(',',':'))+'\\n').encode(),mtime=0)
    if args.record_to is not None:
        args.record_to.mkdir()
        (args.record_to/'PAIR_SOURCE.json.gz').write_bytes(payload)
    else:
        expected = json.loads(gzip.decompress((args.verify/'PAIR_SOURCE.json.gz').read_bytes()))
        actual = json.loads(json.dumps(result))
        for record in (actual, expected):record.pop('seconds')
        _require(actual == expected, 'Saved complete pair evidence differs')
    summary = {k:result[k] for k in ('candidate','passed','literal_native_calls','seconds')}
    summary['evidence_sha256'] = hashlib.sha256(payload).hexdigest()
    if args.record_to is not None:
        (args.record_to/'PAIR_RECEIPT.json').write_text(json.dumps(summary,indent=2)+'\\n')
    print(json.dumps(summary,indent=2))
    raise SystemExit(0 if result['passed'] else 1)
'''
(new/'pair_source.py').write_text(text)
for row in changes:row['descendant_sha256']=hashlib.sha256((HERE/row['descendant']).read_bytes()).hexdigest()
(HERE/'PRODUCER_BUILD.json').write_text(json.dumps(changes,indent=2)+'\n')
