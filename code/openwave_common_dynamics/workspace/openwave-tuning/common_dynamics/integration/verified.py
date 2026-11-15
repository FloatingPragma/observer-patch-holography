"""Exact-byte admission for the additive joint native program."""
from pathlib import Path
import hashlib
import json
import joint_program

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
EXPECTED={
    'vector_field':'aa4a18d89c1a8893ee86eb5e9151bc75e8a0b3a527640cee66aa25382ca55e76',
    'particles':'678983fa40a2a97b8eb399d6b5ad4ab6ace77a5f441baf448cb0c00bc8a6f9a9',
    'spin':'65075e0bc14f883c5cdb877171650b08ee57c8d1bdffce281e68937bb07eb401',
}


def require_digest(payload,expected):
    if hashlib.sha256(payload).hexdigest()!=expected:
        raise ValueError('source bytes differ from the reviewed frozen member')


def verify_sources():
    pins={}
    for directory,expected in EXPECTED.items():
        path=HERE.parent/directory/'MANIFEST.json'
        payload=path.read_bytes();require_digest(payload,expected)
        pins[str(path.relative_to(ROOT))]=expected
        data=json.loads(payload)
        for key in ('files','source_pins','source_sha256'):
            for relative,wanted in data.get(key,{}).items():
                require_digest((ROOT/relative).read_bytes(),wanted)
                if relative in pins and pins[relative]!=wanted:
                    raise ValueError('incompatible frozen source pins')
                pins[relative]=wanted
    return pins


def complete_program(**parameters):
    pins=verify_sources()
    yield from joint_program.complete_program(**parameters)
    if verify_sources()!=pins:raise ValueError('source changed during program generation')
    yield {'operation':'terminal_decode','position':'all128 address-parity outputs; flag32 nonphysical addresses',
           'field':'all288 four-bit residues, without success normalization',
           'clock':'configured X/Y native fine C/E records',
           'complements':'all native reset/export/commit failures retain original records and last completed macro',
           'unmeasured':'event/padding quantum degrees remain in the final state',
           'source_verified':True}


def error_contract(**parameters):
    pins=verify_sources()
    result=joint_program.error_contract(**parameters)
    if verify_sources()!=pins:raise ValueError('source changed during error calculation')
    return result


def run():
    pins=verify_sources()
    source=HERE.parent/'particles/massive_coupled.py'
    payload=source.read_bytes();expected=hashlib.sha256(payload).hexdigest()
    require_digest(payload,expected)
    mutant=bytes([payload[0]^1])+payload[1:]
    try:require_digest(mutant,expected)
    except ValueError:rejected=True
    else:raise AssertionError('actual source digest gate accepted mutation')
    parameters=joint_program.admit()
    assert parameters['clock_wire'] not in parameters['flags']
    receipt={'schema':'oph.common.integration.source_admission.v1',
             'checks':{'exact_frozen_source_bytes':True,'one_byte_source_mutation_rejected':rejected,
                       'verified_callable_has_pre_and_post_custody':True,'disjoint_default_allocation':True},
             'interfaces':['verified.complete_program','verified.error_contract'],
             'source_sha256':pins,'official_rows_changed':[]}
    receipt['source_sha256'][str(Path(__file__).relative_to(ROOT))]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    receipt['source_sha256'][str((HERE/'joint_program.py').relative_to(ROOT))]=hashlib.sha256((HERE/'joint_program.py').read_bytes()).hexdigest()
    (HERE/'SOURCE_ADMISSION.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'checks':receipt['checks'],'pins':len(pins)}))


if __name__=='__main__':run()
