#!/usr/bin/env python3
"""Canonical replay of byte-preserved source experiment; no new source law."""
from pathlib import Path
import gzip
import hashlib
import json
import sys
import time

sys.set_int_max_str_digits(0)  # Exact finite inputs below are frozen and bounded.
import source_experiment_v1 as source

HERE=Path(__file__).resolve().parent


def verify():
    original=json.loads((HERE/'PREREG.json').read_text())
    for name,digest in original['files'].items():
        preserved='source_experiment_v1.py' if name=='headline.py' else name
        assert hashlib.sha256((HERE/preserved).read_bytes()).hexdigest()==digest,preserved
    amended=json.loads((HERE/'REPLAY_PREREG.json').read_text())
    for name,digest in amended['files'].items():
        assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,name


if __name__=='__main__':
    verify(); started=time.perf_counter()
    result,records=source.execute()
    raw=source.canonical(records).encode()
    result['records_sha256']=hashlib.sha256(raw).hexdigest()
    # Canonical JSON is the retained representation; tuples encode as lists.
    assert source.canonical(result)==(HERE/'RECEIPT.json').read_text()
    assert gzip.decompress((HERE/'RECORDS.json.gz').read_bytes())==raw
    print(json.dumps({'candidate':result['candidate'],'literal_native_calls':source.CALLS,
                      'contexts':len(result['source_tests']),
                      'native_opposite_start_TV':'0','hold_and_swap_TV':'1',
                      'source_algorithm_unchanged':True,
                      'records_reproduced_byte_for_byte':True,
                      'seconds':time.perf_counter()-started},sort_keys=True))
