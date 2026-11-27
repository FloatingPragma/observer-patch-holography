"""Verify the deposited sources, the current acquisition code, and every local-channel output."""
import json
from hashlib import sha256

from model.experiments.record_provenance import DESCRIPTIONS, verify_deposit
from ..signed_effects.sequential import compare

PRODUCER = 'model/experiments/interval_echo_consensus/local_interface/measure.py'
# Sorted compact JSON of source_sha256 in the original b2bd43af4 acquisition.
DEPOSITED_SOURCES = '06da137e7b152ba2e246da1ceabfa1f74a8c73565cb688a876f87c54c18c6780'
# In this revision only the final, post-acquisition comparison delegates to this verifier.
REVIEWED_PRODUCER = 'aeadc2cf6ef7ca2b51bc3d9a168463151e2276ac12c005169cb48806a30a8799'


def compare_replay(saved, actual):
    expected, live = saved['source_sha256'], actual['source_sha256']
    encoded = json.dumps(expected, sort_keys=True, separators=(',', ':')).encode()
    if sha256(encoded).hexdigest() != DEPOSITED_SOURCES:
        raise ValueError('source manifest differs from the deposit for this local-interface acquisition')
    if live.keys() != expected.keys():
        raise ValueError('local-interface replay source set changed')
    verify_deposit({name: digest for name, digest in expected.items() if name in DESCRIPTIONS or name == PRODUCER})
    if live[PRODUCER] not in (expected[PRODUCER], REVIEWED_PRODUCER):
        raise ValueError('unreviewed local-interface acquisition code')
    for name, digest in expected.items():
        if name not in DESCRIPTIONS and name != PRODUCER and live[name] != digest:
            raise ValueError(f'local-interface executable or carrier changed: {name}')
    compare({k: v for k, v in saved.items() if k != 'source_sha256'},
            {k: v for k, v in actual.items() if k != 'source_sha256'})
