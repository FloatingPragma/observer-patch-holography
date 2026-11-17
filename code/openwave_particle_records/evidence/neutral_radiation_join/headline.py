#!/usr/bin/env python3
"""Replay retained public events and native fine words; do not rerun the source."""
import gzip
import json
from probe import HERE,load_sources,sha,encode,validate


def main():
    for line in (HERE/'PREREG.txt').read_text().splitlines():
        expected,name = line.split('  ',1)
        assert sha((HERE/name).read_bytes()) == expected,name
    manifest = HERE/'MANIFEST.sha256'
    if manifest.exists():
        for line in manifest.read_text().splitlines():
            expected,name = line.split('  ',1)
            assert sha((HERE/name).read_bytes()) == expected,name
    rw,mem = load_sources()
    raw = gzip.decompress((HERE/'RECORDS.json.gz').read_bytes())
    stored = json.loads((HERE/'RECEIPT.json').read_text())
    assert sha(raw) == stored['records']['uncompressed_sha256']
    result = validate(json.loads(raw),rw,mem)
    assert {k:v for k,v in stored.items() if k != 'records'} == result
    assert result['all_checks_pass']
    report = dict(checks=result['checks'],field_cases=result['field']['cases'],
        contact_branches=len(result['field']['flights']),modes=len(result['field']['modes']),
        source_native_calls=result['field']['native_calls'],
        memory_native_calls=result['memory']['record_replay']['attempts'],
        normal_services=result['memory']['normal_services'],
        normal_conclusive=result['memory']['normal_conclusive'],
        forced_timeouts_retained=True,full_physical_rows_solved=[])
    print(json.dumps(report,sort_keys=True,indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
