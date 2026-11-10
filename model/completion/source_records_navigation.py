"""Typed source-record evidence projection from the sealed field checkpoint.

Only four enumerated prose fields, primary evidence, two architecture premises
and their existing binding unions change. Physical statuses and every unlisted
semantic value are conserved. The predecessor validators remain immutable.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from tools import strict_json
from . import field_navigation as retained

ROOT=Path(__file__).resolve().parents[2]
AUTHORITY='model/completion/source_records_navigation_authority.json'
AUTHORITY_SHA256='e52da9fe4c7dd9e612351af36f1bf61e0168a46f1d054fc25e3f0e9fac2a8be0'
REGISTRY=retained.REGISTRY
TRACKING=retained.TRACKING
DICTIONARY=retained.DICTIONARY
JSON_FILES=retained.JSON_FILES
canonical=retained.canonical
digest=retained.digest
bytes_digest=retained.bytes_digest
require=retained.require
at=retained.at
normalize=retained.normalize
project=retained.project


def authority(root=ROOT):
    root=Path(root)
    raw=(root/AUTHORITY).read_bytes()
    require(bytes_digest(raw)==AUTHORITY_SHA256,'source-record navigation authority changed')
    approved=strict_json.loads(raw.decode())
    require(approved['immutable_checkpoint']=='c7b09261853f5a3c5a67e21294c3d9511ed59fc7',
            'retained field checkpoint changed')
    for key in ('predecessor_custody_authority','predecessor_navigation_authority'):
        row=approved[key]
        require(bytes_digest((root/row['path']).read_bytes())==row['sha256'],
                'retained field authority changed')
    for name,row in approved['baseline_files'].items():
        require(bytes_digest((root/row['retained_path']).read_bytes())==row['sha256'],
                'retained field ledger bytes changed: '+name)
    for name,expected in approved['new_evidence_sha256'].items():
        path=(root/name).resolve()
        require(path.is_relative_to(root.resolve()) and path.is_file()
                and bytes_digest(path.read_bytes())==expected,
                'reviewed source-record evidence changed: '+name)
    return approved


def validate(root=ROOT,approved=None):
    root=Path(root);approved=authority(root) if approved is None else approved
    documents={name:strict_json.load(root/name) for name in JSON_FILES}
    normalize(documents,approved)
    raw=(root/DICTIONARY).read_bytes();suffix=approved['dictionary_append'].encode()
    require(bool(suffix) and raw.endswith(suffix),'reviewed architecture definitions changed')
    require(bytes_digest(raw[:-len(suffix)])==approved['baseline_files'][DICTIONARY]['sha256'],
            'retained architecture dictionary changed')
    return {'schema':'oph.source_records_navigation_validation.v1',
            'retained_checkpoint':approved['immutable_checkpoint'],
            'prose_fields_replaced':len(approved['prose_replacements']),
            'added_evidence_files':len(approved['new_evidence_sha256']),
            'added_assumptions':list(approved['new_assumption_definitions']),
            'claim_count':len(documents[REGISTRY]['claims']),
            'physical_statuses_and_nonlisted_semantics_exactly_preserved':True,
            'scientific_producer_replay_claimed_by_this_projection':False,
            'current_files':{name:bytes_digest((root/name).read_bytes())
                             for name in (*JSON_FILES,DICTIONARY)}}


def apply(root=ROOT):
    root=Path(root);approved=authority(root)
    raw={name:(root/name).read_bytes() for name in (*JSON_FILES,DICTIONARY)}
    if all(bytes_digest(value)==approved['baseline_files'][name]['sha256'] for name,value in raw.items()):
        result=project({name:strict_json.loads(raw[name].decode()) for name in JSON_FILES},approved)
        outputs={name:(json.dumps(result[name],indent=2,ensure_ascii=False)+'\n').encode() for name in JSON_FILES}
        outputs[DICTIONARY]=raw[DICTIONARY]+approved['dictionary_append'].encode()
        for name in raw:require((root/name).read_bytes()==raw[name],'ledger changed during projection')
        for name,value in outputs.items():(root/name).write_bytes(value)
    return validate(root,approved)


def main():
    parser=argparse.ArgumentParser(description=__doc__);group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--apply',action='store_true');group.add_argument('--check',action='store_true')
    args=parser.parse_args();print(json.dumps(apply() if args.apply else validate(),indent=2,sort_keys=True))


if __name__=='__main__':main()
