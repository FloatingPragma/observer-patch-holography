"""Build a candidate, independently replay it, and only then replace the receipt."""
import json
from . import source, observables, verify


def build():
    sources = source.build()
    packet = dict(source=sources, observables=observables.build(verify.HERE/'measurements.json', sources))
    verify.evidence(packet)
    return dict(schema='oph-paired-gravity-v1', sources=verify.pins(), claim=verify.claim_digest(), evidence=packet)


if __name__ == '__main__':
    row = build()
    (verify.HERE/'receipt.json').write_text(json.dumps(row, indent=2, sort_keys=True)+'\n', encoding='ascii', newline='\n')
    print('Built and independently verified complete paired-gravity receipt')
